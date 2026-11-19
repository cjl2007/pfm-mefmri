"""Display-only network ordering. Colors and names come from the PFM dlabel table."""
from pathlib import Path
import html,json,re
import nibabel as nib
import numpy as np

def network_order(image, sampled, label_path):
    labels=nib.load(str(label_path))
    if labels.header.get_axis(1)!=image.header.get_axis(1):raise ValueError('PFM labels and grayplot CIFTI geometry differ')
    if labels.shape[0]!=1:raise ValueError('Select a single-map network dlabel')
    axis=labels.header.get_axis(0)
    if not isinstance(axis,nib.cifti2.LabelAxis):raise ValueError('Network file must contain a CIFTI LabelAxis')
    values=np.asarray(labels.dataobj)[0,sampled];table=axis.label[0]
    if not np.isfinite(values).all()or not np.equal(values,values.astype(int)).all():raise ValueError('Invalid network IDs')
    ids=values.astype(int);keys=sorted(set(ids)-{0})+([0]if 0 in ids else [])
    if any(k not in table for k in keys if k!=0):raise ValueError('Label IDs missing from network color table')
    # Keep cortical vertices and volume grayordinates separate; use the same
    # network IDs/colors within each block. Zero stays unassigned in each block.
    brain=image.header.get_axis(1)
    cortex=np.isin(brain.name[sampled],['CIFTI_STRUCTURE_CORTEX_LEFT','CIFTI_STRUCTURE_CORTEX_RIGHT'])
    ordered=[];sections=[];cursor=0
    for name,mask in [('Cortex',cortex),('Subcortex',~cortex)]:
        if not mask.any():continue
        block=np.concatenate([np.flatnonzero(mask & (ids==k))for k in keys])
        ordered.append(block);sections.append(dict(name=name,start=cursor,stop=cursor+len(block),assigned_rows=int(np.count_nonzero(ids[mask])),sampled_rows=int(mask.sum())))
        cursor+=len(block)
    order=np.concatenate(ordered);colors=[];legend=[]
    for k in keys:
        name,rgba=table.get(k,('Unassigned',(0.7,0.7,0.7,1)))
        if k==0:name='Unassigned';rgba=(.7,.7,.7,1)
        color=[float(v)for v in rgba];color[3]=1.0
        count=int((ids==k).sum());legend.append({'id':int(k),'name':name,'rgba':color,'sampled_rows':count,'cortical_rows':int(((ids==k)&cortex).sum()),'subcortical_rows':int(((ids==k)&~cortex).sum())})
    color_table={item['id']:item['rgba'] for item in legend}
    colors=np.asarray([color_table[int(ids[i])] for i in order])
    return {'order':order,'colors':colors,'legend':legend,'compartments':sections,'label_file':str(Path(label_path).resolve())}

def draw_network_strip(ax, info):
    strip=ax.inset_axes([-.013,0,.008,1])
    strip.imshow(info['colors'][:,None,:],aspect='auto',interpolation='nearest',origin='upper')
    sections=info.get('compartments',[])
    if len(sections)>1:
        for section in sections[1:]:
            y=section['start']-.5
            # A display-only gap leaves the sampled signals and metrics intact.
            ax.axhline(y,color='white',lw=2,zorder=5)
            strip.axhline(y,color='white',lw=2,zorder=5)
        for section in sections:
            center=(section['start']+section['stop'])/2/len(info['order'])
            strip.text(-.9,1-center,section['name'],transform=strip.transAxes,ha='right',va='center',rotation=90,fontsize=6,color='#59636c',clip_on=False)
    strip.set_xticks([]);strip.set_yticks([])
    for spine in strip.spines.values():spine.set_visible(False)

def annotate_report(folder,info):
    meta={k:v for k,v in info.items()if k not in ['order','colors']}
    (folder/'network_order.json').write_text(json.dumps(meta,indent=2))
    blocks=[]
    for item in info['legend']:
        color=','.join(str(round(255*x))for x in item['rgba'][:3])
        blocks.append('<span style="display:inline-block;margin:5px"><i style="display:inline-block;width:12px;height:12px;background:rgb('+color+')"></i> '+html.escape(item['name'])+'</span>')
    path=folder/'motion_grayplots_qc.html';text=path.read_text();note='<section id="network-order"><h2>Network order</h2><p>Cortex and subcortex are grouped separately by network, with a thin divider. Colors come from the fitted map; unassigned rows are gray. The order is shared across runs and stages.</p>'+''.join(blocks)+'<p><a href="network_order.json">Network labels and sampled row counts</a></p></section>'
    text=re.sub(r'<section id="network-order">.*?</section>','',text,flags=re.S)
    text=text.replace('<body>','<body>'+note,1)if '<body>'in text else text.replace('</body>',note+'</body>')
    path.write_text(text)
