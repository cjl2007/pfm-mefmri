"""Redraw saved coregistration measurements in the documentation style."""
from pathlib import Path
import os,tempfile,json
os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='pfm-doc-plots-'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import font_manager
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];A=ROOT/'docs/assets/coreg-stability'
try:font_manager.findfont('Arial',fallback_to_default=False);family='Arial'
except ValueError:family='Liberation Sans'
plt.rcParams.update({'font.family':family,'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'normal','pdf.fonttype':42,'svg.fonttype':'none'})
data=pd.read_csv(A/'scan_metrics.csv');sharp=pd.read_csv(A/'refinement_sharpness.csv');colors=['#84919e','#168b80']
legend=[Patch(facecolor=c,label=l) for c,l in zip(colors,['Refinement off','Refinement on'])]
def save(fig,name):
    for ext in ['png','pdf','svg']:fig.savefig(A/(name+'.'+ext),dpi=180,bbox_inches='tight',facecolor='white')
    plt.close(fig)
fig,axes=plt.subplots(4,3,figsize=(15,12),sharex=False,constrained_layout=True)
for i in range(6):
    person=f'Individual {i+1}'
    for j,(method,color) in enumerate(zip(['baseline','template'],colors)):
        group=data[(data.individual==person)&(data.method==method)].sort_values('scan')
        for k,key in enumerate(['adjacent_rms_percent','median_rms_to_other_scans_percent']):
            ax=axes[2*(i//3)+k,i%3];ax.plot(group.scan,group[key],color=color,lw=1.25)
            if j==1:
                for scan in group.loc[group.between_sessions.fillna(False).astype(bool),'scan']:ax.axvline(scan,color='#c7cdd1',lw=.6,zorder=0)
                ax.grid(axis='y',color='#ebedef');ax.set_xlabel('Scan')
                if k==0:ax.set_title(f'{person} · {len(group)} scans',loc='left',fontsize=12)
                if i%3==0:ax.set_ylabel('Consecutive scans (RMS %)' if k==0 else 'Other scans (median RMS %)')
for k in range(2):
    group=[axes[2*(i//3)+k,i%3] for i in range(6)];lo=min(a.get_ylim()[0] for a in group);hi=max(a.get_ylim()[1] for a in group)
    for ax in group:ax.set_ylim(lo,hi)
axes[0,0].legend(handles=legend,frameon=False,fontsize=9)
save(fig,'six_subject_stability')
def boxes(ax,key):
    for i in range(4):
        pair=[]
        for j,method in enumerate(['baseline','template']):
            if key=='sharpness':values=sharp[(sharp.individual==f'Individual {i+1}')&(sharp.refinement==['Off','On'][j])].sort_values('scan').normalized_gradient_per_mm.to_numpy()
            else:values=data[(data.individual==f'Individual {i+1}')&(data.method==method)].sort_values('scan').median_rms_to_other_scans_percent.to_numpy()
            pair.append(values)
            box=ax.boxplot(values,positions=[3*i+.9*j],widths=.6,patch_artist=True,showfliers=False,manage_ticks=False,medianprops=dict(color='#26333c'))
            box['boxes'][0].set(facecolor=colors[j],edgecolor=colors[j],alpha=.65)
        for a,b in zip(*pair):ax.plot([3*i,3*i+.9],[a,b],color='#bbc3c9',lw=.5,alpha=.18,zorder=0)
    ax.set_xticks([3*i+.45 for i in range(4)]);ax.set_xticklabels([f'Individual {i+1}\n{len(data[(data.individual==f"Individual {i+1}")&(data.method=="baseline")])} scans' for i in range(4)])
    ax.grid(axis='y',color='#ebedef');ax.set_ylabel('Spatial gradient (1/mm)' if key=='sharpness' else 'Other scans (median RMS %)')
fig,axes=plt.subplots(1,2,figsize=(13,4.4),constrained_layout=True);boxes(axes[0],'stability');boxes(axes[1],'sharpness');axes[0].legend(handles=legend,frameon=False,ncol=2,fontsize=10);save(fig,'refinement_stability_sharpness')
fig,ax=plt.subplots(figsize=(9,4.2),constrained_layout=True);boxes(ax,'stability');ax.legend(handles=legend,frameon=False,ncol=2);save(fig,'refinement_spread')
# Check the frozen scan values against the previously published summaries.
summary=pd.read_csv(A/'subject_summaries.csv')
for row in summary[summary.endpoint=='median_rms_to_other_scans_percent'].itertuples():
    values=data[(data.individual==row.subject)&(data.method==row.method)].median_rms_to_other_scans_percent
    np.testing.assert_allclose(np.median(values),row.median,rtol=0,atol=1e-12)
(A/'style_render.json').write_text(json.dumps(dict(font=family,records=len(data),original_summary_medians_verified=True,registration_recomputed=False),indent=2)+'\n')
print('Three coregistration figures redrawn; saved summary medians verified')
