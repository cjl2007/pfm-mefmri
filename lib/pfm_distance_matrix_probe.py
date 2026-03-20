#!/usr/bin/env python3
"""Probe/validate DistanceMatrix behavior without loading full v7.3 MAT in Python.

Computes selected cortical rows using the same rules as make_distance_matrix.m:
  - within-hemisphere cortical: geodesic distance on midthickness surface
  - cross-hemisphere cortical: constant (default 255) or Euclidean (optional)
  - cortical <-> subcortical: Euclidean distance in mm

Optionally compares those rows against an existing MATLAB DistanceMatrix.mat by
calling MATLAB to extract only the probed rows.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import nibabel as nib
import numpy as np
from nibabel.affines import apply_affine


@dataclass
class Grayordinates:
    n_gray: int
    left_idx: np.ndarray
    right_idx: np.ndarray
    left_vertices: np.ndarray
    right_vertices: np.ndarray
    coords_mm: np.ndarray


def _as_uint8_like_matlab(x: np.ndarray) -> np.ndarray:
    # MATLAB uint8 conversion rounds to nearest integer and saturates to [0,255].
    return np.clip(np.rint(x), 0, 255).astype(np.uint8, copy=False)


def load_grayordinates(ref_cifti: Path, left_surf: Path, right_surf: Path) -> Grayordinates:
    img = nib.load(str(ref_cifti))
    bm = img.header.get_axis(1)
    if not hasattr(bm, "iter_structures"):
        raise RuntimeError("Expected BrainModelAxis as CIFTI axis 1")

    l_surf = nib.load(str(left_surf))
    r_surf = nib.load(str(right_surf))
    l_coords = np.asarray(l_surf.darrays[0].data, dtype=np.float64)
    r_coords = np.asarray(r_surf.darrays[0].data, dtype=np.float64)

    n_gray = bm.size
    coords_mm = np.full((n_gray, 3), np.nan, dtype=np.float64)

    left_idx: List[int] = []
    right_idx: List[int] = []
    left_vertices: List[int] = []
    right_vertices: List[int] = []

    for name, slc, sub in bm.iter_structures():
        stop = slc.stop if slc.stop is not None else n_gray
        gidx = np.arange(slc.start, stop, dtype=np.int32)
        if name == "CIFTI_STRUCTURE_CORTEX_LEFT":
            v = np.asarray(sub.vertex, dtype=np.int32)
            coords_mm[gidx, :] = l_coords[v, :]
            left_idx = gidx.tolist()
            left_vertices = v.tolist()
        elif name == "CIFTI_STRUCTURE_CORTEX_RIGHT":
            v = np.asarray(sub.vertex, dtype=np.int32)
            coords_mm[gidx, :] = r_coords[v, :]
            right_idx = gidx.tolist()
            right_vertices = v.tolist()
        else:
            vox = np.asarray(sub.voxel, dtype=np.float64)
            coords_mm[gidx, :] = apply_affine(bm.affine, vox)

    if len(left_idx) == 0 or len(right_idx) == 0:
        raise RuntimeError("Failed to locate both cortical hemispheres in CIFTI axis")
    if np.isnan(coords_mm).any():
        raise RuntimeError("Some grayordinates are missing coordinates")

    return Grayordinates(
        n_gray=n_gray,
        left_idx=np.asarray(left_idx, dtype=np.int32),
        right_idx=np.asarray(right_idx, dtype=np.int32),
        left_vertices=np.asarray(left_vertices, dtype=np.int32),
        right_vertices=np.asarray(right_vertices, dtype=np.int32),
        coords_mm=coords_mm,
    )


def geodesic_row_to_cortex(
    surf_path: Path, source_vertex: int, cortex_vertices: np.ndarray
) -> np.ndarray:
    with tempfile.TemporaryDirectory(prefix="pfm_geo_") as tdir:
        out_shape = Path(tdir) / "dist.shape.gii"
        cmd = [
            "wb_command",
            "-surface-geodesic-distance",
            str(surf_path),
            str(int(source_vertex)),
            str(out_shape),
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        dist_all = np.asarray(nib.load(str(out_shape)).darrays[0].data, dtype=np.float64)
    return dist_all[cortex_vertices]


def build_probe_rows(
    g: Grayordinates,
    left_surf: Path,
    right_surf: Path,
    source_rows: np.ndarray,
    interhemi_mode: str,
    interhemi_constant: int,
) -> np.ndarray:
    left_set = set(g.left_idx.tolist())
    right_set = set(g.right_idx.tolist())
    sub_idx = np.setdiff1d(np.arange(g.n_gray, dtype=np.int32), np.r_[g.left_idx, g.right_idx], assume_unique=True)
    out = np.zeros((source_rows.size, g.n_gray), dtype=np.uint8)

    left_pos = {int(gi): i for i, gi in enumerate(g.left_idx)}
    right_pos = {int(gi): i for i, gi in enumerate(g.right_idx)}

    for ri, s in enumerate(source_rows.tolist()):
        row = np.zeros((g.n_gray,), dtype=np.float64)
        src_xyz = g.coords_mm[s]

        if s in left_set:
            src_vertex = int(g.left_vertices[left_pos[s]])
            row[g.left_idx] = geodesic_row_to_cortex(left_surf, src_vertex, g.left_vertices)
            if interhemi_mode == "constant":
                row[g.right_idx] = float(interhemi_constant)
            else:
                d = np.linalg.norm(g.coords_mm[g.right_idx] - src_xyz[None, :], axis=1)
                row[g.right_idx] = d
        elif s in right_set:
            src_vertex = int(g.right_vertices[right_pos[s]])
            row[g.right_idx] = geodesic_row_to_cortex(right_surf, src_vertex, g.right_vertices)
            if interhemi_mode == "constant":
                row[g.left_idx] = float(interhemi_constant)
            else:
                d = np.linalg.norm(g.coords_mm[g.left_idx] - src_xyz[None, :], axis=1)
                row[g.left_idx] = d
        else:
            # This tool is intentionally cortex-row-focused.
            raise RuntimeError(f"Source row {s} is not cortical")

        dsub = np.linalg.norm(g.coords_mm[sub_idx] - src_xyz[None, :], axis=1)
        row[sub_idx] = dsub
        out[ri, :] = _as_uint8_like_matlab(row)

    return out


def extract_matlab_rows(distance_mat: Path, rows_1based: np.ndarray, matlab_cmd: str) -> np.ndarray:
    with tempfile.TemporaryDirectory(prefix="pfm_matlab_rows_") as tdir:
        rows_txt = Path(tdir) / "rows.txt"
        out_csv = Path(tdir) / "rows.csv"
        script_m = Path(tdir) / "extract_rows.m"
        np.savetxt(rows_txt, rows_1based, fmt="%d")
        rows_txt_esc = str(rows_txt).replace("'", "''")
        dist_mat_esc = str(distance_mat).replace("'", "''")
        out_csv_esc = str(out_csv).replace("'", "''")
        script_m_esc = str(script_m).replace("'", "''")

        script_m.write_text(
            "\n".join(
                [
                    f"rows = readmatrix('{rows_txt_esc}');",
                    f"S = load('{dist_mat_esc}');",
                    "D = S.D;",
                    "R = D(rows, :);",
                    f"writematrix(double(R), '{out_csv_esc}');",
                    "exit(0);",
                ]
            )
        )
        subprocess.run(
            [matlab_cmd, "-batch", f"run('{script_m_esc}')"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        m = np.loadtxt(out_csv, delimiter=",")
        if m.ndim == 1:
            m = m[None, :]
        return m.astype(np.uint8, copy=False)


def main() -> int:
    ap = argparse.ArgumentParser(description="Probe and validate DistanceMatrix rows")
    ap.add_argument("--ref-cifti", required=True)
    ap.add_argument("--left-surf", required=True)
    ap.add_argument("--right-surf", required=True)
    ap.add_argument("--distance-mat", default="", help="Optional MATLAB DistanceMatrix.mat for row-wise comparison")
    ap.add_argument("--matlab-cmd", default="matlab")
    ap.add_argument("--n-left", type=int, default=8)
    ap.add_argument("--n-right", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--interhemi-mode", choices=["constant", "euclidean"], default="constant")
    ap.add_argument("--interhemi-constant", type=int, default=255)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    g = load_grayordinates(Path(args.ref_cifti), Path(args.left_surf), Path(args.right_surf))
    rng = np.random.default_rng(args.seed)
    left_pick = np.sort(rng.choice(g.left_idx, size=min(args.n_left, g.left_idx.size), replace=False))
    right_pick = np.sort(rng.choice(g.right_idx, size=min(args.n_right, g.right_idx.size), replace=False))
    source_rows = np.r_[left_pick, right_pick].astype(np.int32)

    py_rows = build_probe_rows(
        g,
        Path(args.left_surf),
        Path(args.right_surf),
        source_rows=source_rows,
        interhemi_mode=args.interhemi_mode,
        interhemi_constant=args.interhemi_constant,
    )

    summary: Dict[str, object] = {
        "n_gray": int(g.n_gray),
        "source_rows_1based": (source_rows + 1).tolist(),
        "interhemi_mode": args.interhemi_mode,
        "interhemi_constant": int(args.interhemi_constant),
    }

    if args.distance_mat:
        mat_rows = extract_matlab_rows(Path(args.distance_mat), source_rows + 1, args.matlab_cmd)
        diff = py_rows.astype(np.int16) - mat_rows.astype(np.int16)
        abs_diff = np.abs(diff)
        summary.update(
            {
                "compared_to_matlab": True,
                "max_abs_diff": int(abs_diff.max()),
                "mean_abs_diff": float(abs_diff.mean()),
                "fraction_exact": float((abs_diff == 0).mean()),
                "fraction_absdiff_le_1": float((abs_diff <= 1).mean()),
            }
        )
    else:
        summary["compared_to_matlab"] = False

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
