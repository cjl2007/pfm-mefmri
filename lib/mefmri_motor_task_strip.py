"""Motor-only display annotation; shares the recorded acquisition clock with task GLMs."""
import csv,json
from pathlib import Path
import numpy as np
COLORS={'rest':'#b8bec7','index':'#377eb8','thumb':'#e69f00'}
def load_motor_strip(subject, task, run, n_frames, tr):
    if task != 'motor':return None
    raw=Path(subject)/'func/unprocessed/motor'/f'session_{run.session}'/f'run_{run.run}'
    events=raw/f'Motor_S{run.session}_R{run.run}_events.csv'
    if not events.is_file():return None
    result={'segments':[],'notes':[]}
    try:
        meta=json.loads((raw/f'Motor_S{run.session}_R{run.run}_taskmodel.json').read_text())
        if meta.get('alignment_status') not in ('confirmed','assumed_authorized'):
            raise ValueError('timing alignment pending')
        offset=float(meta['onset_offset_seconds'])
        removed=int((run.run_dir/'rmVols.txt').read_text())if(run.run_dir/'rmVols.txt').exists()else 0
        times=np.asarray(json.loads((raw/f'Motor_S{run.session}_R{run.run}_E1.json').read_text())['SliceTiming'])
        shifts=np.loadtxt(run.run_dir/'SliceTiming.txt',ndmin=1)
        if times.shape!=shifts.shape:raise ValueError('slice timing mismatch')
        reference=float(np.median(times+shifts*tr))
        rows=list(csv.DictReader(events.open()))
        end=-np.inf
        for row in rows:
            onset=float(row['onset']);duration=float(row['duration']);condition=row['trial_type']
            if condition not in COLORS or not np.isfinite([onset,duration,offset]).all() or duration<=0 or onset<end-2e-6:
                raise ValueError('invalid event intervals')
            end=onset+duration
            start=(onset+offset-reference)/tr-removed
            stop=start+duration/tr
            a=max(-.5,start);b=min(n_frames-.5,stop)
            if b>a:result['segments'].append((a,b-a,condition))
        result['notes']=[meta['alignment_status'].replace('_',' ')+'; unlogged time left blank']
    except (OSError,ValueError,KeyError,TypeError) as e:
        result={'segments':[],'notes':['Task strip unavailable: '+str(e)]}
    return result

def draw_motor_strip(ax, strip):
    from matplotlib.patches import Patch
    ax.set_facecolor('#fafafa')
    for onset,duration,condition in strip['segments']:
        ax.broken_barh([(onset,duration)],(0,1),facecolors=COLORS[condition],edgecolors='none')
    ax.set_ylim(0,1);ax.set_yticks([]);ax.set_ylabel('task')
    ax.legend(handles=[Patch(color=c,label=k)for k,c in COLORS.items()],loc='upper left',bbox_to_anchor=(0,1.8),ncol=3,fontsize=7,frameon=False)
    ax.text(1,1.12,'; '.join(dict.fromkeys(strip['notes'])),transform=ax.transAxes,ha='right',fontsize=6)
