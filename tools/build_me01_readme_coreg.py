#!/usr/bin/env python3
"""Build a longitudinal coregistration GIF for a supplied example subject."""

from __future__ import annotations

import re
import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from mefmri_qc_coreg_sdc_report import image_from_three_slice, load_3d, make_mask


OUT = ROOT / "docs/assets/qc-coreg-over-time.gif"
RUN_RE = re.compile(r"_(?:AvgFM|ScanSpecificFM)_S(\d+)_R(\d+)\.nii\.gz$")


def selected_runs(coreg_dir: Path) -> list[tuple[int, int, Path]]:
    chosen: dict[tuple[int, int], Path] = {}
    for path in coreg_dir.glob("SBref2acpc_EpiReg+BBR_*_S*_R*.nii.gz"):
        match = RUN_RE.search(path.name)
        if not match:
            continue
        key = (int(match.group(1)), int(match.group(2)))
        if key not in chosen or "ScanSpecificFM" in path.name:
            chosen[key] = path
    return [(session, run, chosen[(session, run)]) for session, run in sorted(chosen)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    SUBJECT = args.subject_dir
    output = args.output
    t1 = load_3d(SUBJECT / "func/xfms/rest/T1w_acpc_brain_func.nii.gz")
    contour = make_mask(t1)
    slices = tuple(int(np.median(np.where(contour)[axis])) for axis in range(3))
    runs = selected_runs(SUBJECT / "func/qa/CoregQA")
    if not runs:
        raise SystemExit("No ME01 coregistration QA volumes found")

    images = [
        image_from_three_slice(
            t1,
            "Coregistration over time",
            contour,
            slices,
            "T1w anatomical anchor",
        )
    ]
    for scan_index, (session, run, path) in enumerate(runs, start=1):
        images.append(
            image_from_three_slice(
                load_3d(path),
                "Coregistration over time",
                contour,
                slices,
                f"Scan {scan_index:02d} / {len(runs):02d}",
            )
        )

    width = 1100
    height = round(images[0].height * width / images[0].width)
    images = [im.resize((width, height), Image.Resampling.LANCZOS) for im in images]
    frames: list[Image.Image] = []
    durations: list[int] = []
    for index, current in enumerate(images):
        following = images[(index + 1) % len(images)]
        frames.extend((current, Image.blend(current, following, 0.5)))
        durations.extend((900 if index == 0 else 520, 100))

    output.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        output,
        save_all=True,
        append_images=frames[1:],
        duration=durations,
        loop=0,
        optimize=True,
    )
    print(f"Wrote {output} with {len(runs)} runs across {len({s for s, _, _ in runs})} sessions")


if __name__ == "__main__":
    main()
