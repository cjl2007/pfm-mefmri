"""Replace per-task grayplots with the current PFM network-ordered reports."""
from pathlib import Path
import fcntl,hashlib,html,json,os,re,shutil,subprocess,sys,tempfile

def source_stats(subject,task,qc_name):
    records=[];seen=set()
    for base in [subject/'func'/task,subject/'func/unprocessed'/task,subject/'physio/unprocessed'/task]:
        if not base.exists():continue
        for directory,dirs,files in os.walk(base,followlinks=True):
            real=str(Path(directory).resolve())
            if real in seen:dirs[:]=[];continue
            seen.add(real);dirs[:]=sorted(d for d in dirs if d not in [qc_name,'QC'] and not d.startswith('.MotionGrayPlots'))
            for name in sorted(files):
                p=Path(directory)/name
                if p.is_file():s=p.stat();records.append((str(p),s.st_size,s.st_mtime_ns))
    return records

def publish(stage,target):
    """Replace the generated report only after rendering succeeds; restore on failure."""
    old=target.with_name('.'+target.name+'.previous')
    if target.is_symlink():raise ValueError('Refusing to replace a symlinked report directory: '+str(target))
    if old.exists():
        if not target.exists():old.replace(target)
        else:shutil.rmtree(old)
    if target.exists():target.replace(old)
    try:stage.replace(target)
    except BaseException:
        if old.exists():old.replace(target)
        raise
    if old.exists():shutil.rmtree(old)

def retire_legacy(root,task):
    """Remove recognized generated hash directories after their replacement is ready."""
    folder=root/'QC/NetworkGrayPlots'/task;retired=[]
    if (root/'QC/NetworkGrayPlots').is_symlink():return retired
    if folder.is_symlink():return retired
    if folder.is_dir():
        for p in folder.iterdir():
            if p.is_symlink() or not p.is_dir() or not re.fullmatch('[a-f0-9]{16}',p.name):continue
            marker=p/'complete.json'
            if not marker.is_file() or not (p/'motion_grayplots_qc.html').is_file():continue
            try:record=json.loads(marker.read_text())
            except (ValueError,OSError):continue
            if not all(k in record for k in ['label_sha256','files','code']):continue
            shutil.rmtree(p);retired.append(p.name)
        if not any(folder.iterdir()):folder.rmdir()
    return retired

def run(medir,subject,start,env=None,runner=subprocess.run):
    env=dict(os.environ if env is None else env);medir=Path(medir);subject=Path(subject);root=Path(env['PFM_DATASET_ROOT'])
    label=Path(env.get('PFM_NETWORK_GRAYPLOTS_LABEL') or root/'RidgePFM'/(env.get('PFM_RF_OUTFILE','RidgePFM_VTX')+'.dlabel.nii'))
    if not env.get('PFM_NETWORK_GRAYPLOTS_LABEL') and env.get('PFM_STRATEGY','ridgepfm').lower()=='infomap':label=root/'Infomap'/(env.get('PFM_INFOMAP_LABEL_OUTFILE','InfomapNetworkLabels')+'_ModeConsensus.dlabel.nii')
    if not label.is_file():raise ValueError('Network grayplots require a fitted network dlabel: '+str(label))
    renderer=medir/'lib/mefmri_qc_motion_grayplots.py';specs=json.loads(env['PFM_DATASET_TASKS_JSON']);qc_name=env.get('QC_OUTPUT_DIRNAME','QC');outputs=[]
    for spec in specs:
        task=spec['task'];prefix=spec['prefix']
        if not task or Path(task).name!=task or task in ['.','..']:raise ValueError('Invalid task directory')
        parent=subject/'func'/task/qc_name;parent.mkdir(parents=True,exist_ok=True);out=parent/'MotionGrayPlots'
        codes=[renderer,medir/'lib/mefmri_network_grayplot.py',medir/'lib/mefmri_motor_task_strip.py',medir/'lib/mefmri_motion.py',medir/'lib/mefmri_figure_style.py',Path(__file__)]
        provenance=dict(label=str(label.resolve()),label_sha256=hashlib.sha256(label.read_bytes()).hexdigest(),files=source_stats(subject,task,qc_name),code={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in codes if p.is_file()},max_rows=env.get('MOTION_GRAYPLOTS_QC_MAX_ROWS','1200'),fd_threshold=env.get('MOTION_GRAYPLOTS_QC_FD_THRESHOLD','0.3'),qcfc_sample=env.get('MOTION_GRAYPLOTS_QC_QCFC_SAMPLE_SIZE','300'),start_session=str(start))
        fingerprint=hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest()
        with (parent/'.MotionGrayPlots.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            marker=out/'complete.json';current={}
            if marker.exists():
                try:current=json.loads(marker.read_text())
                except (ValueError,OSError):pass
            if current.get('fingerprint')==fingerprint and (out/'motion_grayplots_qc.html').exists() and (out/'network_order.json').exists():
                print('[network-grayplots] reuse '+str(out),flush=True)
            else:
                with tempfile.TemporaryDirectory(prefix='.MotionGrayPlots-',dir=parent) as temp:
                    stage=Path(temp)/'report';child=env.copy();child['MEFMRI_GRAYPLOT_LOCK_HELD']='1'
                    runner([sys.executable,str(renderer),'--subject-dir',str(subject),'--func-name',task,'--func-prefix',prefix,'--start-session',str(start),'--out-dir',str(stage),'--network-labels',str(label),'--max-grayplot-rows',provenance['max_rows'],'--fd-threshold',provenance['fd_threshold'],'--qcfc-sample-size',provenance['qcfc_sample']],env=child,check=True)
                    if not (stage/'motion_grayplots_qc.html').is_file() or not (stage/'network_order.json').is_file():raise ValueError('Incomplete network grayplot rendering')
                    if source_stats(subject,task,qc_name)!=provenance['files'] or hashlib.sha256(label.read_bytes()).hexdigest()!=provenance['label_sha256']:raise ValueError('Grayplot inputs changed while rendering; previous report retained')
                    (stage/'complete.json').write_text(json.dumps(dict(status='complete',fingerprint=fingerprint,provenance=provenance),indent=2)+'\n')
                    publish(stage,out)
                print('[network-grayplots] updated '+str(out),flush=True)
            retired=retire_legacy(root,task)
        outputs.append(dict(task=task,report=str(out/'motion_grayplots_qc.html'),fingerprint=fingerprint,retired_legacy_reports=retired))
    qc=root/'QC';qc.mkdir(parents=True,exist_ok=True)
    from mefmri_figure_style import report_css
    links=''.join('<li><a href="'+html.escape(os.path.relpath(v['report'],qc),quote=True)+'">'+html.escape(v['task'])+'</a></li>' for v in outputs)
    (qc/'MotionGrayPlots.html').write_text('<!doctype html><meta charset="utf-8"><title>Motion and physiology</title><style>'+report_css()+'</style><main><h1>Motion and physiology</h1><ul>'+links+'</ul></main>')
    (qc/'network_grayplots.json').write_text(json.dumps(dict(label=str(label),outputs=outputs),indent=2)+'\n')
    legacy=qc/'NetworkGrayPlots';index=legacy/'index.html'
    if index.is_file() and 'PFM network-ordered grayplots' in index.read_text():index.unlink()
    if legacy.is_dir() and not any(legacy.iterdir()):legacy.rmdir()
    return outputs

if __name__=='__main__':
    try:run(sys.argv[1],sys.argv[2],sys.argv[3])
    except (ValueError,OSError,subprocess.CalledProcessError) as e:raise SystemExit('[network-grayplots] '+str(e))
