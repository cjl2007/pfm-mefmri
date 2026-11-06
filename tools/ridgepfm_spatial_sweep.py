#!/usr/bin/env python3
"""Render a spatial-weight sweep from saved ridge coefficients on a checked axis."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import nibabel as nib
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from pfm_ridge_fusion import clamp_prob_smooth, inject_subcortical_spatial_priors, spatial_prior_term
from pfm_ridge_scene import build_scene


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coefficients', required=True, help='NPZ with TAXONOMY_beta arrays (grayordinates x networks)')
    parser.add_argument('--reference', required=True, help='Matching WTA at ridge=1, FC=1, spatial=.05; checked exactly')
    parser.add_argument('--taxonomy', default='Lynch2024')
    parser.add_argument('--subcortical-priors', default='')
    parser.add_argument('--subject-dir', required=True)
    parser.add_argument('--work-dir', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    work = Path(args.work_dir); work.mkdir(parents=True, exist_ok=True)
    subject = Path(args.subject_dir)
    ref = nib.load(args.reference); axis = ref.header.get_axis(1)
    beta = np.load(args.coefficients)[args.taxonomy+'_beta'].astype(float)
    ngray, nnet = beta.shape
    assert ngray == len(axis) and np.isfinite(beta).all()
    prior = loadmat(root/'res0urces/PriorsPFM'/args.taxonomy/'PriorsCortical.mat', squeeze_me=True, struct_as_record=False)['Priors']
    spatial = np.full((ngray, nnet), .5, np.float32)
    spatial[axis.surface_mask] = prior.Spatial[:, :nnet]
    inject_subcortical_spatial_priors(axis, spatial, args.subcortical_priors, nnet)
    spatial = spatial_prior_term(clamp_prob_smooth(spatial), 'logit')
    fc = (beta-beta.mean(1, keepdims=True))/np.maximum(beta.std(1, keepdims=True), 1e-8)
    labels = lambda weight: (np.argmax(fc+weight*spatial, axis=1)+1).astype(np.int32)
    np.testing.assert_array_equal(labels(.05), np.asarray(ref.dataobj)[0])
    weights = [0., .05, .1, .25, .5]
    frames = []; records = []
    env = os.environ.copy(); env.update(OMP_NUM_THREADS='4', QT_QPA_PLATFORM='offscreen')
    fontpath = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
    baseline = labels(0)
    for weight in weights:
        name = 'spatial_'+str(weight).replace('.', 'p')
        data = labels(weight)
        cifti = work/(name+'.dlabel.nii')
        nib.save(nib.Cifti2Image(data[None], ref.header, nifti_header=ref.nifti_header), cifti)
        scene = work/(name+'.scene'); png = work/(name+'.png')
        borders=[]
        for hemi in ['L','R']:
            border=work/(name+'.'+hemi+'.border');borders.append(border)
            surface=subject/'anat/T1w/fsaverage_LR32k'/(subject.name+'.'+hemi+'.midthickness.32k_fs_LR.surf.gii')
            subprocess.run(['wb_command','-cifti-label-to-border',str(cifti),'-border',str(surface),str(border)],check=True,env=env)
        build_scene(root/'res0urces/CrawlingSeedFC/FlatMaps+Inflated.scene', cifti,
                    subject/'anat/T1w/fsaverage_LR32k', subject/'anat/MNINonLinear/fsaverage_LR32k', subject.name, scene, 0, border_paths=borders, compact=True)
        with (work/(name+'.log')).open('w') as log:
            subprocess.run(['wb_command', '-scene-capture-image', str(scene), '1', str(png),
                            '-size-width-height', '1200', '720'], env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        brain = Image.open(png).convert('RGB')
        assert np.asarray(brain).std() > 5, 'Blank scene'
        frame = Image.new('RGB', (1200, 880), 'white'); frame.paste(brain, (0, 100))
        draw = ImageDraw.Draw(frame)
        draw.text((28, 16), f'RidgePFM · {args.taxonomy}', font=ImageFont.truetype(fontpath, 28), fill='#142d45')
        draw.text((28, 58), f'Spatial prior = {weight:g}'+(' · FC only' if weight == 0 else '')+'     |     FC = 1 · ridge = 1',
                  font=ImageFont.truetype(fontpath, 23), fill='#334b61')
        changed = float(np.mean(data[axis.surface_mask] != baseline[axis.surface_mask]))
        draw.text((28, 837), f'{changed:.1%} of cortical assignments differ from FC only · same fitted FC evidence throughout',
                  font=ImageFont.truetype(fontpath, 19), fill='#334b61')
        frame.save(work/(name+'_annotated.png')); frames.append(frame)
        records.append(dict(spatial_weight=weight, cortical_fraction_changed_from_fc_only=changed, source=str(cifti)))
        print('Spatial weight', weight, 'rendered', flush=True)
    output = Path(args.output); output.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(output, save_all=True, append_images=frames[1:], duration=2400, loop=0, disposal=2)
    (work/'sweep.json').write_text(json.dumps(dict(taxonomy=args.taxonomy, ridge=1, fc_weight=1,
        coefficients=str(Path(args.coefficients).resolve()), coefficients_sha256=hashlib.sha256(Path(args.coefficients).read_bytes()).hexdigest(),
        default_labels_match_reference_exactly=True, frames=records,
        interpretation='Increasing the spatial weight strengthens group-prior influence. This sensitivity analysis does not establish an optimal weight.'), indent=2)+'\n')


if __name__ == '__main__':
    main()
