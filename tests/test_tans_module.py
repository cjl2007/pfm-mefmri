#!/usr/bin/env python3

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class TansModuleTests(unittest.TestCase):
    def test_robustness_default_and_explicit_disabled(self):
        for value in ('', '0', '1'):
            with self.subTest(value=value):
                result = self.run_module(TANS_ROBUSTNESS_ENABLE=value)
                self.assertEqual(result.returncode, 0, result.stderr)
                args = self.arguments()
                self.assertEqual(args[args.index('--robustness-enabled') + 1], value or '1')
        result = self.run_module(TANS_ROBUSTNESS_ENABLE='false')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('TANS_ROBUSTNESS_ENABLE', result.stderr)

    def setUp(self):
        self.module = Path(os.environ.get(
            "TEST_TANS_MODULE", Path(__file__).parents[1] /
            "modules" / "mefmri_func_tans.sh"))
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.study = self.root / "study with spaces"
        for subject in ("CIR002", "BioSelect-23002"):
            subject_dir = self.study / subject
            for modality in ("T1w", "T2w"):
                folder = subject_dir / "anat" / modality
                folder.mkdir(parents=True)
                (folder / f"{modality}.nii.gz").write_bytes(b"fixture")
            (subject_dir / "func" / "rest" / "PFM" / "RidgePFM").mkdir(
                parents=True)
        code = self.root / "code"
        code.mkdir()
        (code / "run_pipeline.py").write_text(
            "import json,sys\n"
            "with open(sys.argv[1] + '.args', 'w') as handle:\n"
            "    json.dump(sys.argv[1:], handle)\n")
        self.config = self.root / "tans.yaml"
        self.config.write_text("subject: '{subject}'\n")
        self.capture = self.config.with_suffix(".yaml.args")
        self.environment = {
            key: value for key, value in os.environ.items()
            if not key.startswith("TANS_")}
        self.environment.update({
            "TANS_CODE_ROOT": str(code),
            "TANS_CONFIG_FILE": str(self.config),
            "TANS_PYTHON": "python3",
            "TANS_WORKERS": "4",
            "TANS_TARGET_NETWORK": "Salience",
            "TANS_COARSE_GRID_DENSITY": "0.90",
            "TANS_TARGET_SEARCH_SPACE": str(self.root / "search-space.dscalar.nii"),
            "TANS_INCLUDE_SM_HAND_CALIBRATION": "0",
            "TANS_SM_HAND_COIL_MODEL": str(self.root / "C-B60.ccd"),
            "TANS_TREATMENT_COIL_MODEL": str(self.root / "Cool-B65.ccd"),
            "TANS_CALIBRATION_REGISTRY": str(self.root / "registry.json"),
            "TANS_SM_HAND_DIDT_CALIBRATION": "c60-calibration",
            "TANS_TREATMENT_DIDT_CALIBRATION": "b65-ap-calibration",
        })
        for name in ("search-space.dscalar.nii", "C-B60.ccd", "Cool-B65.ccd"):
            (self.root / name).write_bytes(b"fixture")
        (self.root / "registry.json").write_text("{}")

    def run_module(self, subject="CIR002", **overrides):
        if self.capture.exists():
            self.capture.unlink()
        return subprocess.run(
            ["bash", str(self.module), subject, str(self.study), "/pipeline", "1",
             "rest", "Rest"], capture_output=True, text=True,
            env={**self.environment, **overrides})

    def arguments(self):
        return json.loads(self.capture.read_text())

    def assert_output(self, expected):
        arguments = self.arguments()
        self.assertEqual(arguments[arguments.index("--output") + 1], str(expected))

    def test_module_passes_native_subject_and_pfm_paths(self):
        completed = self.run_module()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        arguments = self.arguments()
        subject_dir = self.study / "CIR002"
        self.assertIn(str(subject_dir), arguments)
        self.assertIn(str(subject_dir / "func/rest/PFM/RidgePFM"), arguments)
        self.assert_output(subject_dir / "tans")
        self.assertIn(f"[tans] output={subject_dir / 'tans'}", completed.stdout)
        for value in ("Salience", "0.90", "--skip-sm-hand-calibration",
                      "c60-calibration", "b65-ap-calibration"):
            self.assertIn(value, arguments)
        for name in ("search-space.dscalar.nii", "C-B60.ccd", "Cool-B65.ccd",
                     "registry.json"):
            self.assertIn(str(self.root / name), arguments)
        self.assertNotIn("anat/charm", completed.stdout)

    def test_output_default_follows_each_subject(self):
        for subject in ("CIR002", "BioSelect-23002"):
            with self.subTest(subject=subject):
                completed = self.run_module(subject)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assert_output(self.study / subject / "tans")

    def test_legacy_folder_fallback_and_new_folder_precedence(self):
        pfm = self.study / 'CIR002/func/rest/PFM'
        (pfm/'RidgePFM').rename(pfm/'RidgeFusion')
        completed = self.run_module()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        arguments = self.arguments()
        self.assertEqual(arguments[arguments.index('--pfm-dir')+1], str(pfm/'RidgeFusion'))
        (pfm/'RidgePFM').mkdir()
        completed = self.run_module()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        arguments = self.arguments()
        self.assertEqual(arguments[arguments.index('--pfm-dir')+1], str(pfm/'RidgePFM'))
        completed = self.run_module(TANS_PFM_DIR=str(pfm/'RidgeFusion'))
        self.assertEqual(completed.returncode, 0, completed.stderr)
        arguments = self.arguments()
        self.assertEqual(arguments[arguments.index('--pfm-dir')+1], str(pfm/'RidgeFusion'))

    def test_unspecified_hardware_preserves_study_yaml_choices(self):
        completed = self.run_module(
            TANS_SM_HAND_COIL_MODEL="", TANS_TREATMENT_COIL_MODEL="",
            TANS_SM_HAND_DIDT_CALIBRATION="", TANS_TREATMENT_DIDT_CALIBRATION="")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        args = self.arguments()
        for option in ("--sm-hand-coil", "--treatment-coil",
                       "--sm-hand-didt-calibration", "--treatment-didt-calibration"):
            self.assertNotIn(option, args)
        self.assertIn("study YAML", completed.stdout)

    def test_empty_folder_name_defaults_to_tans(self):
        completed = self.run_module(TANS_OUTPUT_DIRNAME="")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assert_output(self.study / "CIR002" / "tans")

    def test_custom_folder_name_follows_each_subject(self):
        for subject in ("CIR002", "BioSelect-23002"):
            for folder in ("tans_review-v2", "tans review"):
                with self.subTest(subject=subject, folder=folder):
                    completed = self.run_module(subject, TANS_OUTPUT_DIRNAME=folder)
                    self.assertEqual(completed.returncode, 0, completed.stderr)
                    self.assert_output(self.study / subject / folder)

    def test_folder_name_cannot_be_a_path(self):
        for folder in ("/tmp/tans", "../other", "nested/tans", ".", "..",
                       "tans/", "nested\\tans", "tans\nother", "tans\rother"):
            with self.subTest(folder=folder):
                completed = self.run_module(TANS_OUTPUT_DIRNAME=folder)
                self.assertEqual(completed.returncode, 2)
                self.assertIn("single folder name", completed.stderr)
                self.assertFalse(self.capture.exists())

    def test_legacy_output_path_requires_explicit_migration(self):
        completed = self.run_module(TANS_OUTPUT_DIR="/tmp/old-tans",
                                    TANS_OUTPUT_DIRNAME="tans")
        self.assertEqual(completed.returncode, 2)
        self.assertIn("TANS_OUTPUT_DIR is no longer supported", completed.stderr)
        self.assertFalse(self.capture.exists())

    def test_empty_legacy_setting_does_not_override_folder_name(self):
        completed = self.run_module(TANS_OUTPUT_DIR="", TANS_OUTPUT_DIRNAME="review")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assert_output(self.study / "CIR002" / "review")

    def test_bundled_configs_use_subject_relative_output(self):
        config_root = Path(os.environ.get(
            "TEST_TANS_CONFIG_DIR", Path(__file__).parents[1] / "config"))
        for name in ("mefmri_wrapper_config.sh", "mefmri_wrapper_config_BD2-SAL.sh",
                     "mefmri_wrapper_config_BioSelect-AMN.sh"):
            if name != "mefmri_wrapper_config.sh" and not (config_root / name).exists():
                continue  # Study-specific configs are not distributed in the release.
            with self.subTest(config=name):
                text = (config_root / name).read_text()
                self.assertRegex(text, r'(?m)^TANS_OUTPUT_DIRNAME="tans"\s')
                self.assertNotRegex(text, r"(?m)^TANS_OUTPUT_DIR=")


if __name__ == "__main__":
    unittest.main()
