"""Run the real coregistration dispatch/worker with synthetic files and tool stubs.

GNU parallel and filesystem run discovery are real; no MRI registration is run.
The extracted blocks deliberately exclude anatomy generation and cleanup.
"""
import ast
from dataclasses import dataclass
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(os.environ.get('MEFMRI_TEST_ROOT', Path(__file__).resolve().parents[1]))

STUB = r'''
import os
from pathlib import Path
import sys

root = Path(os.environ['COREG_TEST_ROOT']).resolve()
flags = {'--out', '--fslregout', '--reg', '--o', '-out', '-omat'}
args = sys.argv[1:]
for index, arg in enumerate(args):
    name, sep, value = arg.partition('=')
    if name not in flags:
        continue
    value = value if sep else args[index + 1]
    path = Path(value).resolve()
    if root not in path.parents:
        raise SystemExit('Stub output outside its temporary fixture: ' + str(path))
    outputs = [path]
    if Path(sys.argv[0]).name == 'epi_reg_dof':
        outputs += [Path(value + suffix) for suffix in ('.nii.gz', '_warp.nii.gz', '.mat')]
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.touch()
'''


@unittest.skipUnless(shutil.which('parallel'), 'GNU parallel is required')
class CoregRunCoverageTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='coreg-coverage-')
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        source = (ROOT / 'modules/mefmri_func_coreg.sh').read_text()
        start = source.index('# Also coregister individual SBrefs to the target anatomical image.')
        end = source.index('# Write run-level pointer files needed by later stages.', start)
        pointer_start = source.index('# Write run-level intermediate target and warp pointers')
        pointer_end = source.index('# Remove the temporary FreeSurfer folder.', pointer_start)
        self.blocks = source[start:end] + '\n' + source[pointer_start:pointer_end]

    def run_fixture(self, layout=None, fieldmaps=(), start_session=1,
                    scan_specific=1, spaced_paths=False):
        layout = layout or {s: (1, 2) for s in range(1, 14)}
        suffix = ' with spaces' if spaced_paths else ''
        subject = self.base / ('study' + suffix) / ('SUB001' + suffix)
        code = self.base / ('code' + suffix)
        func = subject / 'func/rest'
        qa = func / 'QC/Coreg'
        xfms = subject / 'func/xfms/rest'
        fm = func / 'field_maps/AllFMs'
        stub_bin = self.base / 'stub-bin'
        for folder in (qa, xfms, fm, stub_bin, code / 'res0urces'):
            folder.mkdir(parents=True)
        for name in ('AvgSBref.nii.gz', 'AvgSBref2acpc_EpiReg+BBR_warp.nii.gz',
                     'AvgSBref2nonlin_EpiReg+BBR_warp.nii.gz'):
            (xfms / name).touch()
        for session, runs in layout.items():
            for run in runs:
                path = func / f'session_{session}/run_{run}'
                path.mkdir(parents=True)
                (path / 'SBref.nii.gz').touch()
                (path / 'PE.txt').write_text('-y\n')
                (path / 'EffectiveEchoSpacing.txt').write_text('0.0005\n')
                if (session, run) in fieldmaps:
                    for prefix in ('FM_rads_acpc', 'FM_mag_acpc', 'FM_mag_acpc_brain'):
                        (fm / f'{prefix}_S{session}_R{run}.nii.gz').touch()
        executable = self.base / 'stub-tool'
        executable.write_text('#!' + sys.executable + '\n' + STUB)
        executable.chmod(0o755)
        for tool in ('flirt', 'applywarp', 'bbregister', 'tkregister2', 'convertwarp'):
            (stub_bin / tool).symlink_to(executable)
        (code / 'res0urces/epi_reg_dof').symlink_to(executable)
        variables = dict(
            MEDIR=code, Subdir=subject, Subject=subject.name, FuncDirName='rest',
            FuncXfmsDir='rest', StartSession=start_session, NTHREADS=3, DOF=6,
            AtlasTemplate=self.base / ('atlas' + suffix + '.nii.gz'),
            AtlasSpace='T1w', AvgPEDIR='-y', FM_DIR=fm.parent,
            COREG_QA_DIR=qa, WMSEG_NII=self.base / 'wmseg.nii.gz',
            SCAN_SPECIFIC_FM=scan_specific, FUNC_NOFIELDMAP_MODE=0,
        )
        script = 'set -euo pipefail\n' + '\n'.join(
            f'{name}={shlex.quote(str(value))}' for name, value in variables.items())
        script += '\nparallel() { command parallel --plain --will-cite "$@"; }\n'
        script += self.blocks
        environment = os.environ.copy()
        for key in ('BASH_ENV', 'ENV', 'PARALLEL', 'COREG_POINTER_POLICY'):
            environment.pop(key, None)
        environment.update(
            PATH=str(stub_bin) + os.pathsep + environment['PATH'],
            SHELL='/bin/bash', PARALLEL_HOME=str(self.base / 'parallel-state'),
            COREG_TEST_ROOT=str(self.base),
        )
        result = subprocess.run(['bash', '-c', script], env=environment,
                                capture_output=True, text=True, timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout + '\n' + result.stderr)
        expected = {(s, r) for s, runs in layout.items() if s >= start_session for r in runs}
        return subject, qa, xfms, expected

    def output_runs(self, folder, pattern):
        return {tuple(map(int, re.search(r'_S(\d+)_R(\d+)', path.name).groups()))
                for path in folder.glob(pattern)}

    def assert_average_coverage(self, qa, expected):
        for pattern in ('SBref2AvgSBref_S*_R*.mat',
                        'SBref2acpc_EpiReg+BBR_AvgFM_S*_R*.nii.gz',
                        'SBref2nonlin_EpiReg+BBR_AvgFM_S*_R*.nii.gz'):
            self.assertEqual(self.output_runs(qa, pattern), expected, pattern)

    def assert_pointers(self, subject, expected, scan_specific):
        for session, run in expected:
            run_dir = subject / f'func/rest/session_{session}/run_{run}'
            target = Path((run_dir / 'IntermediateCoregTarget.txt').read_text().strip())
            warp = Path((run_dir / 'Intermediate2ACPCWarp.txt').read_text().strip())
            self.assertTrue(target.is_file(), str(target))
            self.assertTrue(warp.is_file(), str(warp))
            if (session, run) in scan_specific:
                self.assertEqual(target, run_dir / 'SBref.nii.gz')
                self.assertEqual(warp.name, f'SBref2acpc_EpiReg+BBR_S{session}_R{run}_warp.nii.gz')
            else:
                self.assertEqual(target.name, 'AvgSBref.nii.gz')
                self.assertEqual(warp.name, 'AvgSBref2acpc_EpiReg+BBR_warp.nii.gz')

    def test_all_26_runs_across_13_sessions_have_qc_inputs(self):
        subject, qa, _, expected = self.run_fixture()
        self.assertEqual(len(expected), 26)
        self.assert_average_coverage(qa, expected)
        self.assert_pointers(subject, expected, set())

    def test_scan_specific_fieldmaps_in_later_sessions_and_runs_are_used(self):
        fieldmaps = {(2, 2), (10, 1), (13, 2)}
        subject, qa, xfms, expected = self.run_fixture(fieldmaps=fieldmaps)
        self.assert_average_coverage(qa, expected)
        self.assertEqual(self.output_runs(qa, 'SBref2acpc_EpiReg+BBR_ScanSpecificFM_S*_R*.nii.gz'), fieldmaps)
        self.assertEqual(self.output_runs(xfms, 'SBref2acpc_EpiReg+BBR_S*_R*_warp.nii.gz'), fieldmaps)
        self.assert_pointers(subject, expected, fieldmaps)

    def test_average_pointer_policy_stays_average(self):
        subject, qa, _, expected = self.run_fixture(
            layout={1: (1, 2), 2: (1, 2)}, fieldmaps={(1, 2), (2, 1)}, scan_specific=0)
        self.assert_average_coverage(qa, expected)
        self.assert_pointers(subject, expected, set())

    def test_start_session_limits_all_later_sessions_not_just_the_first(self):
        subject, qa, _, expected = self.run_fixture(start_session=10)
        self.assertEqual(expected, {(s, r) for s in range(10, 14) for r in (1, 2)})
        self.assert_average_coverage(qa, expected)
        self.assert_pointers(subject, expected, set())
        self.assertFalse((subject / 'func/rest/session_1/run_1/IntermediateCoregTarget.txt').exists())

    def test_nonconsecutive_sessions_and_runs(self):
        subject, qa, _, expected = self.run_fixture(layout={2: (1, 3), 5: (2,), 12: (1, 4)})
        self.assert_average_coverage(qa, expected)
        self.assert_pointers(subject, expected, set())

    def test_paths_with_spaces_keep_parallel_arguments_intact(self):
        subject, qa, _, expected = self.run_fixture(
            layout={1: (1, 2), 10: (1, 3)}, fieldmaps={(10, 3)}, spaced_paths=True)
        self.assert_average_coverage(qa, expected)
        self.assert_pointers(subject, expected, {(10, 3)})

    def test_both_reports_can_discover_every_generated_run(self):
        subject, qa, _, expected = self.run_fixture()
        # Compile the real filesystem discovery without importing plotting libraries.
        report = ROOT / 'lib/mefmri_qc_coreg_sdc_report.py'
        tree = ast.parse(report.read_text())
        selected = [node for node in tree.body
                    if isinstance(node, (ast.ClassDef, ast.FunctionDef))
                    and node.name in {'RunInfo', '_read_text', 'discover_runs'}]
        namespace = dict(Path=Path, os=os, re=re, dataclass=dataclass)
        exec(compile(ast.Module(body=selected, type_ignores=[]), str(report), 'exec'), namespace)
        old_qc = os.environ.pop('QC_OUTPUT_DIRNAME', None)
        try:
            discovered = namespace['discover_runs'](subject, 'rest', 'rest', 'topup')
        finally:
            if old_qc is not None:
                os.environ['QC_OUTPUT_DIRNAME'] = old_qc
        self.assertEqual({(run.session, run.run) for run in discovered}, expected)
        self.assert_average_coverage(qa, expected)


if __name__ == '__main__':
    unittest.main()
