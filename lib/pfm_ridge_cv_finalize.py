#!/usr/bin/env python3
"""Verify that ordinary RidgePFM output uses the CV-selected cortical labels."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import nibabel as nib
from pfm_ridge_cv import write,digest

def finalize(cv_root,map_path,mapping_path):
    root=Path(cv_root);selected=json.loads((root/'selected_parameters.json').read_text());run=Path(selected['run_directory'])
    actual=nib.load(map_path);expected=nib.load(run/'RidgePFM_CV_Selected_Cortex.dlabel.nii');axis=actual.header.get_axis(1);ci=np.flatnonzero(axis.surface_mask)
    if axis[ci]!=expected.header.get_axis(1):raise ValueError('Ordinary mapper changed cortical axis')
    np.testing.assert_array_equal(np.asarray(actual.dataobj)[0,ci],np.asarray(expected.dataobj)[0])
    mapping=json.loads(Path(mapping_path).read_text());parameters=mapping['parameters']
    if parameters['lam']!=selected['lambda_value'] or parameters['spatial_weight']!=selected['spatial_weight']:raise ValueError('Ordinary mapper used different parameters')
    cv=dict(schema=selected['schema'],signature=selected['signature'],criterion=selected['criterion'],max_fc_loss=selected['max_fc_loss'],selected_parameters=str(root/'selected_parameters.json'),report=str(root/'report.html'),cortical_labels_verified=True)
    mapping['cross_validation']=cv;write(mapping_path,mapping)
    selected.update(selected_map=str(Path(map_path).resolve()),selected_map_sha256=digest(map_path),selected_refit_verified=True)
    write(root/'selected_parameters.json',selected);write(run/'selected_parameters.json',selected);write(root/'selected_refit_verification.json',cv)
    print('[ridge-cv] ordinary mapper cortical labels exactly match selected all-data grid map')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cv-root',required=True);p.add_argument('--map',required=True);p.add_argument('--mapping',required=True);a=p.parse_args();finalize(a.cv_root,a.map,a.mapping)
