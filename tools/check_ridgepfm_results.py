#!/usr/bin/env python3
"""Check saved RidgePFM confidence products independently from their writer."""
import argparse,json
from pathlib import Path
import nibabel as nib
import numpy as np


def check(folder,prefix='RidgePFM_VTX'):
    folder=Path(folder);conf=folder/'Confidence';metadata=json.loads((conf/'confidence.json').read_text())
    point=nib.load(folder/(prefix+'.dlabel.nii'));axis=point.header.get_axis(1)
    def read(suffix):
        img=nib.load(conf/(prefix+suffix));assert img.header.get_axis(1)==axis
        return np.asarray(img.dataobj)
    sessions=read('_SessionAssignments.dlabel.nii');frequency=read('_SessionAssignmentFrequency.dscalar.nii')
    valid=(sessions>0).sum(0)
    for k,row in enumerate(frequency,1):
        expected=np.divide((sessions==k).sum(0),valid,out=np.full(len(axis),np.nan),where=valid>0)
        np.testing.assert_allclose(row,expected,atol=1e-7,equal_nan=True)
    np.testing.assert_allclose(frequency[:,valid>0].sum(0),1,atol=1e-6)
    split=read('_SplitAssignments.dlabel.nii');r=read('_FCReproducibility.dscalar.nii')
    ids=np.asarray(metadata['sessions']);plan=metadata['signature']['splits']
    minimum=np.min([(sessions[np.isin(ids,group)]>0).sum(0) for pair in plan for group in pair],axis=0)
    eligible=(split>0).all(0)&np.isfinite(r).all(0)&(minimum>=2)
    same=np.all(split==split[0],axis=0);high=np.min(r,axis=0)>=metadata['fc_threshold']
    expected=np.where(eligible,1+(~same).astype(int)+2*(~high).astype(int),0)
    state=read('_AssignmentStates.dlabel.nii')[0];np.testing.assert_array_equal(state,expected)
    membership=read('_SplitLabelMembership.dscalar.nii')
    for k,row in enumerate(membership,1):np.testing.assert_array_equal(row,(split==k).any(0))
    stripe=read('_NetworkStripes.dlabel.nii')[0].astype(int);wta=np.asarray(point.dataobj)[0]
    np.testing.assert_array_equal(stripe[state!=2],wta[state!=2])
    for g in np.flatnonzero(state==2):assert stripe[g] in split[:,g]
    return dict(status='passed',grayordinates=len(axis),sessions=len(ids),networks=len(frequency),minimum_valid_sessions=int(minimum.min()),
        state_counts={str(i):int(np.sum(state==i)) for i in range(5)},
        frequencies_recomputed_from_every_session=True,stripes_restricted_to_state_2=True,all_split_memberships_checked=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder');p.add_argument('--prefix',default='RidgePFM_VTX');a=p.parse_args()
    result=check(a.folder,a.prefix);(Path(a.folder)/'confidence_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
