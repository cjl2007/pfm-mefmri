#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import nibabel as nib
import numpy as np
from nibabel.cifti2 import BrainModelAxis, Cifti2Header, LabelAxis, SeriesAxis


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
from pfm_ridge_repeat_support import (  # noqa: E402
    censored_run_slices,
    load_raw_run_lengths,
    main,
    summarize_run_labels,
)


class RepeatSupportTests(unittest.TestCase):
    def test_censoring_reconstructs_concatenated_run_slices(self) -> None:
        fd = np.asarray([0.1, 0.3, 0.0, 0.2, 0.21, 0.1, 0.0])
        slices, counts = censored_run_slices(fd, [3, 4], 0.2)
        self.assertEqual(counts, [2, 2])
        self.assertEqual([(x.start, x.stop) for x in slices], [(0, 2), (2, 4)])

    def test_scan_info_requires_positive_lengths(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "ScanInfo.json"
            path.write_text(json.dumps([{"n_timepoints": 4}, {"n_timepoints": 6}]))
            lengths, records = load_raw_run_lengths(path)
            self.assertEqual(lengths, [4, 6])
            self.assertEqual(len(records), 2)

    def test_assignment_frequencies_keep_unavailable_outcome_separate(self) -> None:
        labels = np.asarray([[1, 1, 0, 2], [1, 2, 0, 2], [2, 2, 1, 0]], np.uint16)
        point = np.asarray([1, 2, 1, 2], np.uint16)
        result = summarize_run_labels(labels, point, 2)
        self.assertTrue(np.allclose(result["concatenated_label_repeat_support"], [2 / 3, 2 / 3, 1, 1]))
        self.assertTrue(np.allclose(result["unavailable_run_frequency"], [0, 0, 2 / 3, 1 / 3]))
        self.assertTrue(np.allclose(result["conditional_network_frequencies"].sum(0), 1))
        self.assertTrue(np.array_equal(result["available_run_count"], [3, 3, 1, 2]))

    def test_two_run_end_to_end_writes_repeat_support(self) -> None:
        with tempfile.TemporaryDirectory() as folder_name:
            folder = Path(folder_name)
            brain = BrainModelAxis.from_mask(np.ones(3, bool), name="CortexLeft")
            series = SeriesAxis(start=0.0, step=1.0, size=8)
            values = np.vstack(
                [np.zeros((4, 3), np.float32), np.full((4, 3), 10, np.float32)]
            )
            data_path = folder / "input.dtseries.nii"
            nib.save(nib.Cifti2Image(values, Cifti2Header.from_axes((series, brain))), data_path)
            table = {
                0: ("Unavailable", (0.0, 0.0, 0.0, 0.0)),
                1: ("A", (1.0, 0.0, 0.0, 1.0)),
                2: ("B", (0.0, 1.0, 0.0, 1.0)),
            }
            point_path = folder / "point.dlabel.nii"
            nib.save(
                nib.Cifti2Image(
                    np.asarray([[1, 1, 2]], np.uint16),
                    Cifti2Header.from_axes((LabelAxis(["full"], [table]), brain)),
                ),
                point_path,
            )
            fd_path = folder / "FD.txt"
            np.savetxt(fd_path, np.zeros(8))
            scan_path = folder / "ScanInfo.json"
            scan_path.write_text(json.dumps([{"n_timepoints": 4}, {"n_timepoints": 4}]))
            mapper = folder / "fake_mapper.py"
            mapper.write_text(
                "import argparse,pathlib,nibabel as nib,numpy as np\n"
                "from nibabel.cifti2 import Cifti2Header,LabelAxis\n"
                "p=argparse.ArgumentParser();p.add_argument('--in-cifti');p.add_argument('--outdir');p.add_argument('--outfile');a,_=p.parse_known_args()\n"
                "x=nib.load(a.in_cifti); first=float(np.asarray(x.dataobj)[0,0]); labels=np.array([[1,1,2] if first<5 else [1,2,2]],np.uint16)\n"
                "t={0:('Unavailable',(0.,0.,0.,0.)),1:('A',(1.,0.,0.,1.)),2:('B',(0.,1.,0.,1.))};o=pathlib.Path(a.outdir);o.mkdir(exist_ok=True)\n"
                "nib.save(nib.Cifti2Image(labels,Cifti2Header.from_axes((LabelAxis(['run'],[t]),x.header.get_axis(1)))),o/(a.outfile+'.dlabel.nii'))\n"
            )
            out = folder / "repeat"
            argv = [
                "pfm_ridge_repeat_support.py", "--in-cifti", str(data_path),
                "--full-map-dlabel", str(point_path), "--fd-txt", str(fd_path),
                "--scan-info-json", str(scan_path), "--fd-threshold", "0.2",
                "--minimum-frames", "3", "--mapper", str(mapper),
                "--distance-npy", "unused.npy", "--priors-mat", "unused.mat",
                "--left-surf", "unused.L.surf.gii", "--right-surf", "unused.R.surf.gii",
                "--outdir", str(out), "--outfile-prefix", "Synthetic",
            ]
            with mock.patch.object(sys, "argv", argv):
                self.assertEqual(main(), 0)
            support = np.asarray(
                nib.load(out / "Synthetic_ScanRepeatSupport.dscalar.nii").dataobj
            )
            self.assertTrue(np.allclose(support[0], [1.0, 0.5, 1.0]))
            frequencies = np.asarray(
                nib.load(out / "Synthetic_ScanRepeatNetworkFrequencies.dscalar.nii").dataobj
            )
            self.assertTrue(np.allclose(frequencies.sum(0), 1.0))
            status = json.loads((out / "status.json").read_text())
            self.assertEqual(status["state"], "succeeded")
            self.assertEqual(status["eligible_runs"], 2)
            self.assertFalse((out / "working").exists())


if __name__ == "__main__":
    unittest.main()
