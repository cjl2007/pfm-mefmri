#!/usr/bin/env python3

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
HELPER = REPO_ROOT / "lib" / "find_epi_params.py"


class FindEpiParamsNoFieldmapTests(unittest.TestCase):
    def run_helper(self, subject_dir: Path, *extra_args: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            [
                sys.executable,
                str(HELPER),
                "--subdir",
                str(subject_dir),
                "--func-name",
                "rest",
                "--func-prefix",
                "Rest",
                *extra_args,
            ],
            capture_output=True,
            text=True,
        )

    def make_subject(self, root: Path) -> Path:
        subject_dir = root / "SUB001"
        run_dir = subject_dir / "func" / "unprocessed" / "rest" / "session_1" / "run_1"
        run_dir.mkdir(parents=True)
        metadata = {
            "EchoTime": 0.012,
            "RepetitionTime": 1.5,
            "SliceTiming": [0.0, 0.75],
            "PhaseEncodingDirection": "j-",
        }
        (run_dir / "Rest_S1_R1_E1.json").write_text(json.dumps(metadata), encoding="utf-8")
        return subject_dir

    def test_no_fieldmap_mode_accepts_missing_effective_echo_spacing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            subject_dir = self.make_subject(Path(tmpdir))
            result = self.run_helper(subject_dir, "--no-fieldmap-mode")
            self.assertEqual(result.returncode, 0, result.stderr)
            run_dir = subject_dir / "func" / "rest" / "session_1" / "run_1"
            self.assertEqual((run_dir / "TE.txt").read_text(), "12\n")
            self.assertFalse((run_dir / "EffectiveEchoSpacing.txt").exists())

    def test_missing_effective_echo_spacing_fails_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            subject_dir = self.make_subject(Path(tmpdir))
            result = self.run_helper(subject_dir)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("EffectiveEchoSpacing", result.stderr)


if __name__ == "__main__":
    unittest.main()
