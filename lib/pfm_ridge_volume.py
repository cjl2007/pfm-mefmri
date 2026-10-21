#!/usr/bin/env python3
"""Ribbon-constrained cortical export plus CIFTI volume components on a shared grid."""
from __future__ import annotations
import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to


def run(cmd):
    env=os.environ.copy();env.setdefault('OMP_NUM_THREADS','4')
    subprocess.run([str(x) for x in cmd],check=True,env=env)


def _label_table_file(path, table):
    lines=[]
    for key,(name,rgba) in table.items():
        if int(key)==0: continue
        color=[int(round(255*float(value))) for value in rgba]
        lines.extend([str(name),f'{int(key)} {color[0]} {color[1]} {color[2]} {color[3]}'])
    path.write_text('\n'.join(lines)+'\n')


def export(cifti, reference, surfaces, output, target_reference=None, acpc_to_target_xfm=None):
    image=nib.load(cifti);axis=image.header.get_axis(1);ref=nib.load(reference)
    label=isinstance(image.header.get_axis(0),nib.cifti2.LabelAxis)
    if not isinstance(image.header.get_axis(0),(nib.cifti2.LabelAxis,nib.cifti2.ScalarAxis,nib.cifti2.SeriesAxis)):
        raise ValueError('Expected dense label or scalar maps')
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='ridge_volume_') as folder:
        tmp=Path(folder);parts=[];supports=[]
        for hemi,structure in [('L','CORTEX_LEFT'),('R','CORTEX_RIGHT')]:
            gift=tmp/(hemi+('.label.gii' if label else '.func.gii'));vol=tmp/(hemi+'.nii.gz')
            run(['wb_command','-cifti-separate',cifti,'COLUMN','-label' if label else '-metric',structure,gift])
            mid,white,pial=surfaces[hemi]
            run(['wb_command','-label-to-volume-mapping' if label else '-metric-to-volume-mapping',gift,mid,reference,vol,
                 '-ribbon-constrained',white,pial,'-voxel-subdiv','3'])
            values=nib.load(vol).get_fdata(dtype=np.float32)
            if values.ndim==3:values=values[...,None]
            if label:
                # Workbench can renumber labels when the source's zero label
                # differs from its default unassigned label. Restore original
                # identities by the exported label names before combining.
                for column, original in enumerate(image.header.get_axis(0).label):
                    names={entry[0]:key for key,entry in original.items()}
                    if len(names)!=len(original):
                        raise ValueError('Volume label export requires unique names within each label table')
                    table_path=tmp/(hemi+'_labels_%d.txt'%column)
                    run(['wb_command','-volume-label-export-table',vol,str(column+1),table_path])
                    lines=table_path.read_text().splitlines()
                    corrected=np.zeros(values.shape[:3],np.float32)
                    for line in range(0,len(lines),2):
                        name=lines[line];key=int(lines[line+1].split()[0])
                        if name in names:
                            corrected[values[...,column]==key]=names[name]
                        elif key!=0 and np.any(values[...,column]==key):
                            raise ValueError('Workbench introduced an unknown label: '+name)
                    values[...,column]=corrected
            # Map a valid-surface mask through the exact same ribbon weights.
            bm=next(b for name,sl,b in axis.iter_structures() if name=='CIFTI_STRUCTURE_'+structure)
            mask=np.zeros(axis.nvertices['CIFTI_STRUCTURE_'+structure],np.float32);mask[bm.vertex]=1
            mg=tmp/(hemi+'_support.func.gii');mv=tmp/(hemi+'_support.nii.gz')
            support_image=nib.gifti.GiftiImage(darrays=[nib.gifti.GiftiDataArray(mask,intent='NIFTI_INTENT_SHAPE')])
            support_image.meta['AnatomicalStructurePrimary']='CortexLeft' if hemi=='L' else 'CortexRight'
            nib.save(support_image,mg)
            run(['wb_command','-set-structure',mg,structure])
            run(['wb_command','-metric-to-volume-mapping',mg,mid,reference,mv,'-ribbon-constrained',white,pial,'-voxel-subdiv','3'])
            support=nib.load(mv).get_fdata(dtype=np.float32)
            if support.ndim==4:support=support[...,0]
            parts.append(values);supports.append(support)
        shape=ref.shape[:3];nmap=image.shape[0]
        combined=np.zeros(shape+(nmap,),np.float32);coverage=np.zeros(shape,np.uint8)
        if label:
            for values,support in zip(parts,supports):
                keep=np.any(values!=0,axis=-1)&(support>0)&(coverage==0)
                combined[keep]=values[keep];coverage[keep]=1
        else:
            weight=sum(supports)
            np.divide(sum(parts),weight[...,None],out=combined,where=weight[...,None]>1e-8)
            coverage[weight>1e-8]=1
        # Preserve original volumetric values. This is resampling in the SAME
        # physical space; no registration between anatomical spaces is inferred.
        ids=np.flatnonzero(axis.volume_mask)
        if len(ids):
            native=np.zeros(axis.volume_shape+(nmap,),np.float32);mask=np.zeros(axis.volume_shape,np.uint8)
            vox=axis.voxel[ids];loc=tuple(vox.T);native[loc]=np.asarray(image.dataobj)[:,ids].T;mask[loc]=1
            target=(shape,ref.affine)
            vm=resample_from_to(nib.Nifti1Image(mask,axis.affine),target,order=0).get_fdata()>0
            values=np.stack([resample_from_to(nib.Nifti1Image(native[...,i],axis.affine),target,order=0 if label else 1).get_fdata(dtype=np.float32) for i in range(nmap)],axis=-1)
            if not label:
                volume_weight=resample_from_to(nib.Nifti1Image(mask.astype(np.float32),axis.affine),target,order=1).get_fdata(dtype=np.float32)
                np.divide(values,volume_weight[...,None],out=values,where=volume_weight[...,None]>1e-8)
            combined[vm]=values[vm];coverage[vm]=2
        if label:combined=np.rint(combined).astype(np.int32)
        acpc=tmp/'mapped_acpc.nii.gz';acpc_coverage=tmp/'coverage_acpc.nii.gz'
        nib.save(nib.Nifti1Image(combined,ref.affine),acpc)
        nib.save(nib.Nifti1Image(coverage,ref.affine),acpc_coverage)
        target_reference=Path(target_reference) if target_reference else Path(reference)
        target_coverage=output.with_name(output.name.replace('.nii.gz','_Coverage.nii.gz'))
        if acpc_to_target_xfm:
            mapped=tmp/'mapped_target.nii.gz'
            run(['flirt','-in',acpc,'-ref',target_reference,'-applyxfm','-init',acpc_to_target_xfm,
                 '-interp','nearestneighbour' if label else 'trilinear','-out',mapped])
            run(['flirt','-in',acpc_coverage,'-ref',target_reference,'-applyxfm','-init',acpc_to_target_xfm,
                 '-interp','nearestneighbour','-out',target_coverage])
            if not label:
                # Linear interpolation against zero-filled background attenuates
                # continuous values at the transformed tissue boundary. Map a
                # binary valid-support image with the same interpolation and
                # divide it out, then honor the nearest-neighbor coverage mask.
                acpc_support=tmp/'support_acpc.nii.gz';target_support=tmp/'support_target.nii.gz'
                nib.save(nib.Nifti1Image((coverage>0).astype(np.float32),ref.affine),acpc_support)
                run(['flirt','-in',acpc_support,'-ref',target_reference,'-applyxfm','-init',acpc_to_target_xfm,
                     '-interp','trilinear','-out',target_support])
                values=nib.load(mapped).get_fdata(dtype=np.float32)
                support=nib.load(target_support).get_fdata(dtype=np.float32)
                target_valid=nib.load(target_coverage).get_fdata(dtype=np.float32)>0
                if values.ndim==3: values=values[...,None]
                np.divide(values,support[...,None],out=values,where=support[...,None]>1e-8)
                values[~target_valid]=0
                nib.save(nib.Nifti1Image(values,nib.load(target_reference).affine),output)
        else:
            mapped=acpc if label else output
            if not label: shutil.copyfile(acpc,output)
            shutil.copyfile(acpc_coverage,target_coverage)
        if label:
            tables=image.header.get_axis(0).label
            if len(tables)!=1:
                raise ValueError('Label-volume export currently requires a single-map dense label file')
            table_path=tmp/'output_labels.txt';_label_table_file(table_path,tables[0])
            run(['wb_command','-volume-label-import',mapped,table_path,output,'-discard-others'])
        meta=dict(source=str(Path(cifti).resolve()),reference=str(target_reference.resolve()),
                  ribbon_mapping_reference=str(Path(reference).resolve()),
                  method='Workbench ribbon-constrained; cortical scalar values normalized by mapped valid-surface coverage; volume component takes precedence at overlap; transformed scalar targets are normalized by transformed valid-tissue support.',
                  coverage={'0':'outside mapped tissue','1':'cortical ribbon','2':'CIFTI volume'},
                  maps=[str(n) for n in getattr(image.header.get_axis(0),'name',range(nmap))],
                  space='Ribbon mapping is performed in ACPC anatomical space, then optionally transformed to the configured anatomical target.',
                  acpc_to_target_transform=str(Path(acpc_to_target_xfm).resolve()) if acpc_to_target_xfm else None)
        if label:meta['label_tables']=[{str(k):[v[0],list(v[1])] for k,v in table.items()} for table in image.header.get_axis(0).label]
        output.with_name(output.name.replace('.nii.gz','.json')).write_text(json.dumps(meta,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cifti',required=True);p.add_argument('--reference',required=True)
    p.add_argument('--surface-dir',required=True);p.add_argument('--subject',required=True);p.add_argument('--output',required=True)
    p.add_argument('--target-reference');p.add_argument('--acpc-to-target-xfm')
    a=p.parse_args();surfaces={h:[Path(a.surface_dir)/(a.subject+'.'+h+'.'+kind+'.32k_fs_LR.surf.gii') for kind in ['midthickness','white','pial']] for h in ['L','R']}
    export(a.cifti,a.reference,surfaces,a.output,a.target_reference,a.acpc_to_target_xfm)
if __name__=='__main__':main()
