"""Redraw the frozen cohort summaries; no FC calculation or map fitting."""
from pathlib import Path
import os,json
ROOT=Path(__file__).resolve().parents[1]
import tempfile
os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='ridgepfm-plots-'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
try:font_manager.findfont("Arial",fallback_to_default=False);doc_font="Arial"
except ValueError:doc_font="Liberation Sans"
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
D=json.loads((ROOT/'docs/assets/ridgepfm/ridgepfm-cv-grid-validation.json').read_text())
rows=D['summary'];colors=D['colors'];weights=D['weights'];lambdas=D['lambdas'];markers=['o','s','D','^','v','P','X','*']
assets=ROOT/'docs/assets/ridgepfm';assets.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':doc_font,'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
specs=[('ridgepfm_cv_tradeoff','fc','dice','Held-out FC profile homogeneity\n(above spatial null)','Network map reliability (Dice)'),
       ('ridgepfm_individual_information','fc','specificity','Held-out FC profile homogeneity\n(above spatial null)','Individual specificity (own − other Dice)'),
       ('ridgepfm_reliability_specificity','specificity_N20','dice_N20','Individual specificity (own − other Dice)','Network map reliability (Dice)'),
       ('ridgepfm_parameter_ranges','fc','dice','Held-out FC profile homogeneity\n(above spatial null)','Network map reliability (Dice)')]
for name,xf,yf,xlabel,ylabel in specs:
    fig,axes=plt.subplots(1,2,figsize=(12,6),sharex=True,sharey=True)
    fig.subplots_adjust(left=.085,right=.985,bottom=.32,top=.90,wspace=.23)
    xs=[r[xf] for r in rows];ys=[r[yf] for r in rows];dx=max(xs)-min(xs);dy=max(ys)-min(ys)
    for ax,regime,title in zip(axes,['clinical','dense'],['~30 minutes','Multiple hours']):
        ax.grid(color='#e9ecef',lw=.65,zorder=0)
        for j,lam in enumerate(lambdas):
            group=sorted([r for r in rows if r['regime']==regime and r['lam']==lam],key=lambda r:r['spatial'])
            ax.plot([r[xf] for r in group],[r[yf] for r in group],color='#a5adb3',lw=.8,zorder=1)
            for r in group:ax.scatter(r[xf],r[yf],s=55 if j==7 else 36,marker=markers[j],color=colors[weights.index(r['spatial'])],edgecolor='white',lw=.5,zorder=3)
        ax.set_xlim(min(xs)-.08*dx,max(xs)+.08*dx);ax.set_ylim(min(ys)-.08*dy,max(ys)+.08*dy)
        if name=='ridgepfm_parameter_ranges':
            zone=[r for r in rows if r['regime']==regime and r['lam']>=(1 if regime=='clinical' else 3) and r['spatial']>=(.025 if regime=='clinical' else .01) and r['spatial']<=.1]
            xlo=min(r[xf] for r in zone)-.003;xhi=max(r[xf] for r in zone)+.003
            ylo=min(r[yf] for r in zone)-.01;yhi=max(r[yf] for r in zone)+.01
            ax.add_patch(Rectangle((xlo,ylo),xhi-xlo,yhi-ylo,facecolor='#e7edf1',alpha=.65,edgecolor='#425d70',linewidth=1.3,linestyle='--',zorder=2))
            ax.text(xlo,yhi+.009,'Suggested starting point',fontsize=9,color='#425d70')
        ax.set_title(title,loc='left',fontsize=12,fontweight='normal',pad=12)
        ax.set_xlabel(xlabel,labelpad=9);ax.set_ylabel(ylabel,labelpad=8);ax.tick_params(labelleft=True)
    spatial=[Line2D([0],[0],marker='o',color='none',markerfacecolor=c,markeredgecolor='white',markersize=8,label=f'{w:g}') for c,w in zip(colors,weights)]
    ridge=[Line2D([0],[0],marker=m,color='none',markerfacecolor='#6f7980',markeredgecolor='white',markersize=8,label=f'{l:g}') for m,l in zip(markers,lambdas)]
    fig.legend(handles=spatial,title='Spatial weight',loc='lower center',bbox_to_anchor=(.53,.13),ncol=6,frameon=False,columnspacing=1.9,handletextpad=.4)
    fig.legend(handles=ridge,title='λ',loc='lower center',bbox_to_anchor=(.53,.015),ncol=8,frameon=False,columnspacing=1.6,handletextpad=.4)
    for ext in ['png','svg','pdf']:fig.savefig(assets/(name+'.'+ext),dpi=200,bbox_inches='tight',facecolor='white')
    plt.close(fig)
(assets/'ridgepfm-figure-render.json').write_text(json.dumps(dict(source=str(ROOT/'docs/assets/ridgepfm/ridgepfm-cv-grid-validation.json'),points_per_panel=48,FC_participants=19,map_only_participants=20,statistics_unchanged=True,axes_shared=True),indent=2)+'\n')
print('Four cohort figures redrawn from frozen summaries')
