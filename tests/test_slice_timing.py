import importlib.util
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'lib'))
from find_epi_params import build_slice_timing
from validate_slice_timing import nifti_shape, validate_metadata, validate_file


def header(path, n=60, endian='<', version=1):
    import gzip
    h = bytearray(348 if version == 1 else 540)
    struct.pack_into(endian+'i', h, 0, len(h))
    struct.pack_into(endian+('8h' if version == 1 else '8q'), h, 40 if version == 1 else 16, 4, 2, 2, n, 80, 1, 1, 1)
    if str(path).endswith('.gz'):
        with gzip.open(path, 'wb') as f: f.write(h)
    else: path.write_bytes(h)


def old_timing(times, mb):
    base = times[:len(set(times))]
    order = sorted(range(len(base)), key=lambda i: base[i])
    ref = math.floor(len(base)/2+.5)-1
    return [(ref-order.index(i))/len(base) for i in range(len(base))]*mb


class SliceTimingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        self.times = [0,.8,.1,.9,.2,1,.3,1.1,.4,1.2,.5,1.3,.6,1.4,.7]*4

    def tearDown(self): self.tmp.cleanup()

    def test_missing_mb_preserves_sixty_slices_and_legacy_offsets(self):
        self.assertEqual(build_slice_timing(self.times), old_timing(self.times,4))
        self.assertEqual(len(build_slice_timing(self.times)),60)

    def test_mb_does_not_change_vector_length_or_values(self):
        for mb in [None,1,4,8]:
            self.assertEqual(build_slice_timing(self.times,mb),build_slice_timing(self.times))

    def test_nonrepeated_spatial_grouping_preserved(self):
        times=[0,0,.4,.4,.2,.2]
        self.assertEqual(build_slice_timing(times),[1/3,1/3,-1/3,-1/3,0,0])

    def test_single_band_legacy_unchanged(self):
        for times in [[0,.2,.4,.6],[.6,.4,.2,0],[0,.4,.2,.6]]:
            self.assertEqual(build_slice_timing(times),old_timing(times,1))

    def test_all_invalid_metadata_rejected(self):
        for times in [None,[],[0]*15,[float('nan')]*60,[float('inf')]*60,[-.1]*60,[1.6]*60,['0']*60,[True]*60]:
            with self.subTest(times=str(times)[:30]),self.assertRaises(ValueError):
                validate_metadata(dict(SliceTiming=times,RepetitionTime=1.6),60)

    def test_missing_mb_valid_and_unsupported_axis_fails(self):
        d=dict(SliceTiming=self.times,RepetitionTime=1.6)
        self.assertEqual(validate_metadata(d,60),self.times)
        for direction in ['i','j','k-']:
            with self.assertRaisesRegex(ValueError,'remapping'):
                validate_metadata(dict(d,SliceEncodingDirection=direction),60)

    def test_nifti_headers_and_invalid_headers(self):
        for version in [1,2]:
            for endian in ['<','>']:
                p=self.d/f'{version}{endian=="<"}.nii.gz';header(p,version=version,endian=endian)
                self.assertEqual(nifti_shape(p),(2,2,60,80))
        p=self.d/'bad.nii';p.write_bytes(b'bad')
        with self.assertRaises(ValueError):nifti_shape(p)

    def test_runtime_guard_rejects_short_nan_and_multicolumn(self):
        p=self.d/'x.nii.gz';header(p)
        t=self.d/'SliceTiming.txt'
        for content in ['0\n'*15,'nan\n'*60,'0 0\n'*60,'1\n'*60]:
            t.write_text(content)
            with self.assertRaises(ValueError):validate_file(t,p)
        t.write_text('0\n'*60)
        self.assertEqual(validate_file(t,p),60)

    def prepare_subject(self, invalid=False):
        run=self.d/'s/func/unprocessed/rest/session_1/run_1';run.mkdir(parents=True)
        for e in [1,2]:
            stem=f'Rest_S1_R1_E{e}'
            header(run/(stem+'.nii.gz'))
            (run/(stem+'.json')).write_text(json.dumps(dict(EchoTime=.015*e,RepetitionTime=1.6,SliceTiming=self.times if not invalid or e==1 else self.times[:15])))
        bin_dir=self.d/'bin';bin_dir.mkdir();stub=bin_dir/'fslnvols'
        stub.write_text('#!/bin/sh\necho 80\n');stub.chmod(0o755)
        env=dict(os.environ,PATH=str(bin_dir)+':'+os.environ['PATH'])
        return env

    def run_helper(self,env):
        return subprocess.run([sys.executable,str(ROOT/'lib/find_epi_params.py'),'--subdir',str(self.d/'s'),'--no-fieldmap-mode','--validation-qc','0'],env=env,capture_output=True,text=True)

    def test_end_to_end_metadata_with_missing_mb(self):
        result=self.run_helper(self.prepare_subject())
        self.assertEqual(result.returncode,0,result.stderr)
        out=self.d/'s/func/rest/session_1/run_1'
        self.assertEqual(len((out/'SliceTiming.txt').read_text().splitlines()),60)
        self.assertFalse(json.loads((out/'SliceTiming_provenance.json').read_text())['multiband_factor_used_for_timing'])

    def test_bad_metadata_fails_before_prior_outputs_removed(self):
        env=self.prepare_subject(invalid=True)
        out=self.d/'s/func/rest/session_1/run_1';out.mkdir(parents=True)
        (out/'keep.txt').write_text('keep')
        result=self.run_helper(env)
        self.assertNotEqual(result.returncode,0)
        self.assertEqual((out/'keep.txt').read_text(),'keep')

    def test_both_headmotion_paths_guard_each_slicetimer_call(self):
        for name in ['mefmri_func_headmotion.sh','mefmri_func_headmotion_medic.sh']:
            lines=(ROOT/'modules'/name).read_text().splitlines()
            for i,line in enumerate(lines):
                if line.strip().startswith('slicetimer '):
                    self.assertIn('validate_slice_timing.py',lines[i-1])
                    self.assertIn('|| return 1',lines[i-1])
            joined='\n'.join(lines)
            self.assertLess(joined.index('validate_slice_timing.py'),joined.index('rm -rf "$ProcDir'))

    @unittest.skipUnless(Path('/usr/local/fsl/bin/slicetimer').exists(),'FSL required')
    def test_actual_fsl_corrects_known_synthetic_slice_phase(self):
        import numpy as np
        import nibabel as nib
        ranks=np.array([sorted(set(self.times)).index(t) for t in self.times])
        samples=np.arange(80)
        data=np.broadcast_to(100+np.sin(2*np.pi*.05*(samples[None,:]+ranks[:,None]/15)),(2,2,60,80)).astype('float32').copy()
        img=nib.Nifti1Image(data,np.diag([-1.,1.,1.,1.]));img.header.set_zooms((1,1,1,1.6))
        source=self.d/'synthetic.nii.gz';nib.save(img,source)
        timing=self.d/'timing.txt';timing.write_text('\n'.join(str(x) for x in build_slice_timing(self.times))+'\n')
        validate_file(timing,source)
        output=self.d/'corrected.nii.gz'
        subprocess.run(['/usr/local/fsl/bin/slicetimer','-i',str(source),'-o',str(output),'-r','1.6','--tcustom='+str(timing)],check=True,capture_output=True,env=dict(os.environ,FSLOUTPUTTYPE='NIFTI_GZ'))
        target=100+np.sin(2*np.pi*.05*(samples+7/15))
        before=np.sqrt(np.mean((data[0,0,:,10:-10]-target[None,10:-10])**2))
        after=np.sqrt(np.mean((nib.load(output).get_fdata()[0,0,:,10:-10]-target[None,10:-10])**2))
        self.assertLess(after,before*.1)


if __name__ == '__main__':unittest.main()
