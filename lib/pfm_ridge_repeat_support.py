#!/usr/bin/env python3
"""Estimate RidgePFM assignment repeatability from independently mapped runs.

The full concatenated RidgePFM map remains the point estimate. Each eligible run
is mapped independently with the same estimator, then categorical assignments are
averaged exactly. This is temporal repeat support, not P(true network), fractional
membership, or a calibrated forecast for a future scan.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import nibabel as nib
import numpy as np
from nibabel.cifti2 import Cifti2Header, LabelAxis, ScalarAxis, SeriesAxis


def bool_int(value: str) -> int:
    parsed = int(value)
    if parsed not in (0, 1):
        raise argparse.ArgumentTypeError("expected 0 or 1")
    return parsed


def load_raw_run_lengths(path: Path) -> tuple[list[int], list[dict]]:
    records = json.loads(path.read_text())
    if not isinstance(records, list) or not records:
        raise ValueError("ScanInfo.json must contain a nonempty run list")
    try:
        lengths = [int(record["n_timepoints"]) for record in records]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Every ScanInfo run needs a positive n_timepoints") from exc
    if any(length < 1 for length in lengths):
        raise ValueError("Every ScanInfo run needs a positive n_timepoints")
    return lengths, records


def censored_run_slices(
    fd: np.ndarray, raw_lengths: list[int], threshold: float
) -> tuple[list[slice], list[int]]:
    fd = np.asarray(fd, dtype=np.float64).reshape(-1)
    if sum(raw_lengths) != fd.size:
        raise ValueError(f"run lengths sum to {sum(raw_lengths)}, but FD has {fd.size} rows")
    slices = []
    kept_counts = []
    raw_start = 0
    kept_start = 0
    for length in raw_lengths:
        kept = int(np.count_nonzero(fd[raw_start : raw_start + length] < threshold))
        slices.append(slice(kept_start, kept_start + kept))
        kept_counts.append(kept)
        raw_start += length
        kept_start += kept
    return slices, kept_counts


def summarize_run_labels(
    run_labels: np.ndarray, point_labels: np.ndarray, network_count: int
) -> dict[str, np.ndarray]:
    labels = np.asarray(run_labels)
    point = np.asarray(point_labels)
    if labels.ndim != 2 or point.shape != (labels.shape[1],):
        raise ValueError("run labels and concatenated point labels do not align")
    if labels.dtype.kind not in "iu" or labels.min() < 0 or labels.max() > network_count:
        raise ValueError("run labels are outside the admitted taxonomy")
    available = labels > 0
    available_count = available.sum(0)
    unavailable_frequency = 1.0 - available.mean(0)
    counts = np.zeros((network_count + 1, labels.shape[1]), dtype=np.uint16)
    columns = np.arange(labels.shape[1])
    for row in labels:
        counts[row, columns] += 1
    conditional = np.full((network_count, labels.shape[1]), np.nan, dtype=np.float32)
    np.divide(
        counts[1:].astype(np.float32),
        available_count[np.newaxis, :],
        out=conditional,
        where=available_count[np.newaxis, :] > 0,
    )
    point_support = np.full(labels.shape[1], np.nan, dtype=np.float32)
    valid_point = (point > 0) & (point <= network_count) & (available_count > 0)
    point_count = np.zeros(labels.shape[1], dtype=np.uint16)
    point_count[valid_point] = counts[point[valid_point], columns[valid_point]]
    point_support[valid_point] = point_count[valid_point] / available_count[valid_point]
    modal_label = counts[1:].argmax(0).astype(np.uint16) + 1
    modal_label[available_count == 0] = 0
    modal_support = np.full(labels.shape[1], np.nan, dtype=np.float32)
    valid_modal = available_count > 0
    modal_support[valid_modal] = counts[modal_label[valid_modal], columns[valid_modal]] / available_count[valid_modal]
    entropy = np.full(labels.shape[1], np.nan, dtype=np.float32)
    if network_count > 1:
        safe = np.where(np.isfinite(conditional) & (conditional > 0), conditional, 1.0)
        value = -np.sum(np.where(conditional > 0, conditional * np.log(safe), 0.0), axis=0) / np.log(network_count)
        entropy[valid_modal] = value[valid_modal]
    return {
        "conditional_network_frequencies": conditional,
        "concatenated_label_repeat_support": point_support,
        "competing_label_frequency": 1.0 - point_support,
        "modal_label": modal_label,
        "modal_label_support": modal_support,
        "unavailable_run_frequency": unavailable_frequency.astype(np.float32),
        "available_run_count": available_count.astype(np.float32),
        "normalized_label_entropy": entropy,
        "counts": counts,
    }


def save_scalar(reference, names: list[str], values: np.ndarray, path: Path) -> None:
    axis = reference.header.get_axis(1)
    image = nib.Cifti2Image(
        np.asarray(values, np.float32),
        Cifti2Header.from_axes((ScalarAxis(names), axis)),
        nifti_header=reference.nifti_header,
    )
    image.set_data_dtype(np.float32)
    nib.save(image, path)


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-cifti", required=True)
    parser.add_argument("--full-map-dlabel", required=True)
    parser.add_argument("--fd-txt", required=True)
    parser.add_argument("--scan-info-json", required=True)
    parser.add_argument("--fd-threshold", type=float, required=True)
    parser.add_argument("--minimum-frames", type=int, default=30)
    parser.add_argument("--mapper", required=True)
    parser.add_argument("--distance-npy", required=True)
    parser.add_argument("--priors-mat", required=True)
    parser.add_argument("--subcort-priors-nii", default="")
    parser.add_argument("--fc-weight", type=float, default=1.0)
    parser.add_argument("--fc-demean", type=bool_int, default=0)
    parser.add_argument("--positive-fc-only", type=bool_int, default=0)
    parser.add_argument("--spatial-weight", type=float, default=0.05)
    parser.add_argument("--spatial-prior-transform", choices=("logit", "log"), default="logit")
    parser.add_argument("--lambda", dest="lam", type=float, default=1.0)
    parser.add_argument("--local-exclusion-mm", type=float, default=10.0)
    parser.add_argument("--brain-structures-csv", default="")
    parser.add_argument("--left-surf", required=True)
    parser.add_argument("--right-surf", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--outfile-prefix", default="RidgePFM")
    parser.add_argument("--keep-run-artifacts", type=bool_int, default=0)
    args = parser.parse_args()

    if args.minimum_frames < 3:
        raise ValueError("minimum-frames must be at least 3")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    status = out / "status.json"
    source = nib.load(args.in_cifti)
    point_image = nib.load(args.full_map_dlabel)
    if source.header.get_axis(1) != point_image.header.get_axis(1):
        raise ValueError("concatenated data and full RidgePFM map use different grayordinate axes")
    point_labels = np.asarray(point_image.dataobj)[0].astype(np.uint16)
    label_axis = point_image.header.get_axis(0)
    if not isinstance(label_axis, LabelAxis):
        raise ValueError("full-map-dlabel must use a CIFTI label axis")
    table = label_axis.label[0]
    network_ids = sorted(int(key) for key in table if int(key) > 0)
    if network_ids != list(range(1, len(network_ids) + 1)):
        raise ValueError("network label keys must be consecutive starting at 1")
    network_count = len(network_ids)
    network_names = [str(table[index][0]) for index in network_ids]

    raw_lengths, scan_records = load_raw_run_lengths(Path(args.scan_info_json))
    fd = np.loadtxt(args.fd_txt, dtype=np.float64)
    slices, kept_counts = censored_run_slices(fd, raw_lengths, args.fd_threshold)
    if sum(kept_counts) != source.shape[0]:
        raise ValueError(
            f"censor-derived retained frames {sum(kept_counts)} != input CIFTI rows {source.shape[0]}"
        )
    eligible = [index for index, count in enumerate(kept_counts) if count >= args.minimum_frames]
    metadata = {
        "schema": "ridgepfm_scan_repeat_support_v1",
        "eligible_run_indices_zero_based": eligible,
        "retained_frames_per_run": kept_counts,
        "minimum_frames": args.minimum_frames,
        "fd_threshold": args.fd_threshold,
        "scan_info": scan_records,
        "interpretation": (
            "Independent per-run WTA assignment frequencies around the full concatenated point map. "
            "Not biological membership probabilities, simultaneous memberships, or calibrated future-repeat forecasts."
        ),
    }
    write_json(out / "run_inventory.json", metadata)
    if len(eligible) < 2:
        write_json(status, {"state": "not_estimable", "reason": "fewer_than_two_eligible_runs", **metadata})
        print("[ridge-repeat] fewer than two eligible runs; repeat support is not estimable")
        return 0

    run_labels = []
    run_names = []
    work = out / "working"
    # Keep one compact categorical checkpoint per completed run. A long
    # multi-run job can then resume after interruption without recomputing
    # finished scans. Checkpoints are removed only after the combined CIFTIs
    # and count archive have been written successfully.
    work.mkdir(exist_ok=True)
    try:
        for position, run_index in enumerate(eligible, start=1):
            run_slice = slices[run_index]
            run_name = f"run_{run_index + 1:03d}"
            run_names.append(run_name)
            checkpoint = work / f"{run_name}_labels.npy"
            if checkpoint.exists():
                labels = np.load(checkpoint, allow_pickle=False)
                if labels.shape != point_labels.shape or labels.dtype.kind not in "iu":
                    raise ValueError(f"invalid saved checkpoint: {checkpoint}")
                labels = labels.astype(np.uint16, copy=False)
                if labels.min() < 0 or labels.max() > network_count:
                    raise ValueError(f"saved checkpoint labels outside taxonomy: {checkpoint}")
                run_labels.append(labels)
                write_json(
                    status,
                    {
                        "state": "running",
                        "completed_runs": position,
                        "eligible_runs": len(eligible),
                        "resumed_checkpoint": run_name,
                    },
                )
                continue
            run_dir = work / run_name
            if run_dir.exists():
                shutil.rmtree(run_dir)
            run_dir.mkdir()
            values = np.asarray(source.dataobj[run_slice, :], dtype=np.float32)
            run_series = SeriesAxis(start=0.0, step=source.header.get_axis(0).step, size=values.shape[0])
            run_input = run_dir / "input.dtseries.nii"
            nib.save(
                nib.Cifti2Image(
                    values,
                    Cifti2Header.from_axes((run_series, source.header.get_axis(1))),
                    nifti_header=source.nifti_header,
                ),
                run_input,
            )
            command = [
                sys.executable, args.mapper,
                "--in-cifti", str(run_input),
                "--distance-npy", args.distance_npy,
                "--priors-mat", args.priors_mat,
                "--outdir", str(run_dir),
                "--outfile", "RunMap",
                "--fc-weight", str(args.fc_weight),
                "--fc-demean", str(args.fc_demean),
                "--positive-fc-only", str(args.positive_fc_only),
                "--spatial-weight", str(args.spatial_weight),
                "--spatial-prior-transform", args.spatial_prior_transform,
                "--lambda", str(args.lam),
                "--local-exclusion-mm", str(args.local_exclusion_mm),
                "--fit-r2-mode", "exact" if (args.fc_demean or args.positive_fc_only) else "skip",
                "--brain-structures-csv", args.brain_structures_csv,
                "--left-surf", args.left_surf,
                "--right-surf", args.right_surf,
            ]
            if args.subcort_priors_nii:
                command.extend(["--subcort-priors-nii", args.subcort_priors_nii])
            subprocess.run(command, check=True)
            mapped = nib.load(run_dir / "RunMap.dlabel.nii")
            if mapped.header.get_axis(1) != source.header.get_axis(1):
                raise ValueError(f"{run_name} output changed the native axis")
            labels = np.asarray(mapped.dataobj)[0].astype(np.uint16)
            temporary_checkpoint = checkpoint.with_suffix(".npy.tmp")
            with temporary_checkpoint.open("wb") as stream:
                np.save(stream, labels, allow_pickle=False)
            temporary_checkpoint.replace(checkpoint)
            run_labels.append(labels)
            write_json(
                status,
                {"state": "running", "completed_runs": position, "eligible_runs": len(eligible)},
            )
            if not args.keep_run_artifacts:
                shutil.rmtree(run_dir)

        labels = np.asarray(run_labels, dtype=np.uint16)
        summary = summarize_run_labels(labels, point_labels, network_count)
        per_run = nib.Cifti2Image(
            labels,
            Cifti2Header.from_axes((LabelAxis(run_names, [table] * len(run_names)), source.header.get_axis(1))),
            nifti_header=source.nifti_header,
        )
        per_run.set_data_dtype(np.uint16)
        nib.save(per_run, out / f"{args.outfile_prefix}_PerScan_WTA.dlabel.nii")
        save_scalar(
            source,
            [f"{name} conditional scan-assignment frequency" for name in network_names],
            summary["conditional_network_frequencies"],
            out / f"{args.outfile_prefix}_ScanRepeatNetworkFrequencies.dscalar.nii",
        )
        diagnostic_names = [
            "concatenated-label repeat support given available",
            "competing-label frequency given available",
            "modal-label support given available",
            "unavailable-run frequency",
            "available-run count",
            "normalized label entropy",
        ]
        diagnostic_values = np.vstack(
            [
                summary["concatenated_label_repeat_support"],
                summary["competing_label_frequency"],
                summary["modal_label_support"],
                summary["unavailable_run_frequency"],
                summary["available_run_count"],
                summary["normalized_label_entropy"],
            ]
        )
        save_scalar(
            source,
            diagnostic_names,
            diagnostic_values,
            out / f"{args.outfile_prefix}_ScanRepeatSupport.dscalar.nii",
        )
        np.savez_compressed(
            out / f"{args.outfile_prefix}_ScanRepeatCounts.npz",
            counts=summary["counts"],
            point_labels=point_labels,
            eligible_run_indices=np.asarray(eligible, np.int32),
            retained_frames=np.asarray([kept_counts[index] for index in eligible], np.int32),
        )
        shutil.rmtree(work)
        write_json(
            status,
            {
                "state": "succeeded",
                "eligible_runs": len(eligible),
                "native_targets": int(labels.shape[1]),
                "network_count": network_count,
                **metadata,
            },
        )
        print(f"[ridge-repeat] wrote repeat-support maps from {len(eligible)} eligible runs")
        return 0
    except BaseException as exc:
        write_json(status, {"state": "failed_preserve_working_outputs", "error": repr(exc)})
        raise


if __name__ == "__main__":
    raise SystemExit(main())
