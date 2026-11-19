#!/usr/bin/env python3
"""Run dataset-bound NSI or seed-FC QC, with input/settings/output receipts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from mefmri_pfm_dataset import atomic_json, signature


def run(stage, subject, study, medir, start, task, prefix):
    root = Path(os.environ['PFM_DATASET_ROOT']).resolve()
    manifest = root / 'dataset.json'
    dataset = json.loads(manifest.read_text())
    source = Path(os.environ['PFM_INPUT_CIFTI']).resolve()
    if not source.is_file():
        raise ValueError('Missing selected PFM input: ' + str(source))
    qc = root / 'QC'
    out = qc / ('NSI' if stage == 'nsi' else 'CrawlingSeedFC')
    module = Path(os.environ.get('FUNC_NSI_MODULE' if stage == 'nsi' else 'FUNC_CRAWLING_SEED_FC_MODULE',
                      str(Path(medir) / 'modules' / ('mefmri_func_' + stage + '.sh'))))
    env = os.environ.copy()
    env.update(NSI_INPUT_CIFTI=str(source), NSI_OUTDIR=str(qc/'NSI'),
               NSI_SCAN_INFO_JSON=str(Path(env['PFM_CONCAT_DIR'])/'ScanInfo.json'),
               CRAWLING_SEED_FC_INPUT_CIFTI=str(source), CRAWLING_SEED_FC_OUTDIR=str(qc/'CrawlingSeedFC'))
    setting_prefix = 'NSI_' if stage == 'nsi' else 'CRAWLING_SEED_FC_'
    settings = {k:v for k,v in env.items() if k.startswith(setting_prefix)}
    for key in ('PFM_QC_FORCE', 'CRAWLING_SEED_FC_FORCE', 'CRAWLING_SEED_FC_FORCE_DCONN'):
        settings.pop(key, None)
    helpers = ['mefmri_pfm_qc.py', 'mefmri_pfm_dataset.py']
    if stage != 'nsi': helpers.append('crawling_seed_fc_movie.py')
    code = {n: hashlib.sha256((Path(medir)/'lib'/n).read_bytes()).hexdigest() for n in helpers}
    resources = []
    if stage == 'nsi':
        backend = Path(env.get('NSI_EXTERNAL_ROOT', ''))
        if env.get('NSI_EXTERNAL_ROOT') and backend.is_dir():
            for pattern in ('pfm_nsi/**/*.py', '**/*model*.json.gz'):
                resources.extend(signature(p) for p in sorted(backend.glob(pattern)) if p.is_file())
    else:
        for key, fallback in [('CRAWLING_SEED_FC_SCENE_TEMPLATE','res0urces/CrawlingSeedFC/FlatMaps+Inflated.scene'),
                              ('CRAWLING_SEED_FC_VERTICES','res0urces/CrawlingSeedFC/VerticesToSample.txt')]:
            p = Path(env.get(key) or str(Path(medir)/fallback))
            if p.is_file(): resources.append(signature(p))
        folders = [Path(study)/subject/'anat/T1w/fsaverage_LR32k',
                   Path(study)/subject/'anat/MNINonLinear/fsaverage_LR32k']
        folders += [Path(env[k]) for k in ('CRAWLING_SEED_FC_SURFACE_RESOURCE_DIR', 'CRAWLING_SEED_FC_FLAT_SURFACE_RESOURCE_DIR') if env.get(k)]
        for folder in dict.fromkeys(folders):
            resources.extend(signature(p) for p in sorted(folder.glob('*.surf.gii')))
    binding = dict(dataset_fingerprint=dataset['input_fingerprint'], input=signature(source),
                   settings=settings, module_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),
                   code_sha256=code, resources=resources)
    fingerprint = hashlib.sha256(json.dumps(binding, sort_keys=True).encode()).hexdigest()
    receipt = qc / (stage + '_dataset.json')
    old = json.loads(receipt.read_text()) if receipt.is_file() else {}
    force = env.get('PFM_QC_FORCE','0') == '1' or (stage != 'nsi' and
            any(env.get(k,'0') == '1' for k in ('CRAWLING_SEED_FC_FORCE','CRAWLING_SEED_FC_FORCE_DCONN')))
    intact = old.get('outputs') and all(Path(x['path']).is_file() and signature(x['path']) == x for x in old['outputs'])
    if not force and old.get('fingerprint') == fingerprint and old.get('status') == 'complete' and intact:
        print('[pfm-qc] reusing matching ' + stage + ': ' + str(out), flush=True)
        return
    # A changed dataset must not reuse a previous movie/dconn cached in this directory.
    if stage != 'nsi':
        env.update(CRAWLING_SEED_FC_FORCE='1', CRAWLING_SEED_FC_FORCE_DCONN='1')
    record = dict(schema_version=1, stage=stage, status='running', fingerprint=fingerprint,
                  binding=binding, tasks=dataset.get('tasks', []), output_dir=str(out))
    atomic_json(receipt, record)
    before = {str(p.resolve()): signature(p) for p in out.rglob('*') if p.is_file()}
    try:
        subprocess.run(['bash', str(module), subject, study, medir, str(start), task, prefix], env=env, check=True)
    except Exception:
        atomic_json(receipt, {**record, 'status':'failed'})
        raise
    outputs = [signature(p) for p in sorted(out.rglob('*')) if p.is_file()
               and before.get(str(p.resolve())) != signature(p)]
    atomic_json(receipt, {**record, 'status':'complete' if outputs else 'skipped', 'outputs':outputs})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['nsi','crawling_seed_fc'])
    parser.add_argument('subject'); parser.add_argument('study'); parser.add_argument('medir')
    parser.add_argument('start', type=int); parser.add_argument('task'); parser.add_argument('prefix')
    args = parser.parse_args()
    run(args.stage,args.subject,args.study,args.medir,args.start,args.task,args.prefix)
