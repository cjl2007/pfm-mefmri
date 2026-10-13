#!/usr/bin/env python3
"""Regression tests for per-echo tedana PCA/TR input validation."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATE_MODULE = REPO_ROOT / "modules" / "mefmri_validate_inputs.sh"


class ValidateTedanaPcaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.subject_dir = self.root / "ME01"
        (self.subject_dir / "anat" / "unprocessed" / "T1w").mkdir(parents=True)
        (self.subject_dir / "anat" / "unprocessed" / "T1w" / "T1w.nii.gz").write_bytes(b"anat")

        self.run_dir = (
            self.subject_dir / "func" / "unprocessed" / "rest" / "session_1" / "run_1"
        )
        self.run_dir.mkdir(parents=True)
        for echo in range(1, 5):
            stem = f"Rest_S1_R1_E{echo}"
            # Decreasing placeholder sizes avoid unrelated echo-size warnings.
            (self.run_dir / f"{stem}.nii.gz").write_bytes(b"x" * (5 - echo))
            (self.run_dir / f"{stem}.json").write_text(
                json.dumps(
                    {
                        "EchoTime": echo * 0.01,
                        "RepetitionTime": 1.0,
                        "SliceTiming": [0.0],
                    }
                )
            )

        self.bin_dir = self.root / "bin"
        self.bin_dir.mkdir()
        fslnvols = self.bin_dir / "fslnvols"
        fslnvols.write_text(
            """#!/usr/bin/env bash
case "${1##*/}" in
  *_E1.nii.gz) echo "${TEST_E1_TRS}" ;;
  *_E2.nii.gz) echo "${TEST_E2_TRS}" ;;
  *_E3.nii.gz) echo "${TEST_E3_TRS}" ;;
  *_E4.nii.gz) echo "${TEST_E4_TRS}" ;;
  *) exit 2 ;;
esac
"""
        )
        fslnvols.chmod(0o755)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def run_validation(self, trs: tuple[int, int, int, int], mepca: str = "350") -> subprocess.CompletedProcess[str]:
        env = os.environ.copy()
        env.update(
            {
                "PATH": f"{self.bin_dir}:{env['PATH']}",
                "FUNC_NOFIELDMAP_MODE": "1",
                "MEPCA": mepca,
                "PIPELINE_EFFECTIVE_DENOISE_METHOD": "meica",
                "TEST_E1_TRS": str(trs[0]),
                "TEST_E2_TRS": str(trs[1]),
                "TEST_E3_TRS": str(trs[2]),
                "TEST_E4_TRS": str(trs[3]),
            }
        )
        return subprocess.run(
            ["bash", str(VALIDATE_MODULE), str(self.subject_dir), "rest", "Rest", "1"],
            env=env,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_accepts_exact_fifty_percent_margin_for_every_echo(self) -> None:
        result = self.run_validation((525, 525, 525, 525))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Input validation passed", result.stdout)

    def test_checks_nonfinal_echo_for_numeric_pca_limit(self) -> None:
        result = self.run_validation((600, 524, 600, 600))
        self.assertEqual(result.returncode, 2)
        self.assertIn("Tedana PCA/TR requirement failed for Rest_S1_R1_E2.nii.gz", result.stdout)
        self.assertIn("MEPCA=350 requires at least 525 TRs", result.stdout)

    def test_checks_final_echo_for_numeric_pca_limit(self) -> None:
        result = self.run_validation((600, 600, 600, 500))
        self.assertEqual(result.returncode, 2)
        self.assertIn("Tedana PCA/TR requirement failed for Rest_S1_R1_E4.nii.gz", result.stdout)

    def test_named_pca_method_has_no_fixed_component_limit(self) -> None:
        result = self.run_validation((100, 100, 100, 100), mepca="kundu")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
