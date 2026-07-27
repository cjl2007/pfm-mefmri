#!/usr/bin/env python3
"""Evaluate one-, two-, and three-network parcel descriptions across time splits.

Parcel geometry is held fixed. Candidate network mixtures are selected from two
random temporal thirds and scored against the held-out third. Contiguous blocks
are assigned together to reduce train/test leakage from autocorrelation. A
parcelwise candidate-identity null and cortex/subcortex FDR calibrate the
default descriptive selection. A repeat-level exact sign-flip max-statistic
null is retained as a stricter whole-brain high-confidence tier.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import nibabel as nib
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
from nibabel.cifti2.cifti2_axes import LabelAxis

from pfm_parcel_uncertainty_stripes import (
    parcel_network_from_label,
    save_dscalar,
)


def parse_int_csv(value: str) -> list[int]:
    out = [int(x.strip()) for x in value.split(",") if x.strip()]
    if not out or any(x <= 0 for x in out):
        raise argparse.ArgumentTypeError("expected comma-separated positive integers")
    return out


def load_raw_run_lengths(scan_info_json: Path) -> list[int]:
    with scan_info_json.open() as f:
        scan_info = json.load(f)
    if not isinstance(scan_info, list) or not scan_info:
        raise ValueError(f"Expected a non-empty run list in {scan_info_json}")
    try:
        lengths = [int(run["n_timepoints"]) for run in scan_info]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid n_timepoints entries in {scan_info_json}") from exc
    if any(length <= 0 for length in lengths):
        raise ValueError(f"Run lengths in {scan_info_json} must be positive")
    return lengths


def cortical_grayordinates(brain_axis) -> np.ndarray:
    indices = []
    for name, slc, _ in brain_axis.iter_structures():
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            stop = slc.stop if slc.stop is not None else brain_axis.size
            indices.extend(range(slc.start, stop))
    return np.asarray(indices, dtype=np.int64)


def load_network_priors(priors_mat: Path, cortical_count: int) -> tuple[np.ndarray, list[str]]:
    pri = sio.loadmat(str(priors_mat), squeeze_me=True, struct_as_record=False)["Priors"]
    fc = np.asarray(pri.FC, dtype=np.float64)
    if fc.shape[0] != cortical_count:
        raise ValueError(f"Priors.FC rows ({fc.shape[0]}) != cortical grayordinates ({cortical_count})")
    names = [str(x).strip() for x in np.asarray(pri.NetworkLabels).ravel()[: fc.shape[1]]]
    return fc, names


def censored_run_slices(fd: np.ndarray, raw_run_lengths: list[int], threshold: float) -> list[slice]:
    if int(np.sum(raw_run_lengths)) != fd.size:
        raise ValueError(f"Raw run lengths sum to {sum(raw_run_lengths)}, but FD has {fd.size} values")
    slices = []
    raw_start = 0
    censored_start = 0
    for length in raw_run_lengths:
        kept = int(np.sum(fd[raw_start : raw_start + length] < threshold))
        slices.append(slice(censored_start, censored_start + kept))
        raw_start += length
        censored_start += kept
    return slices


def parcel_membership(parcel_ids: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
    parcel_values = np.unique(parcel_ids[parcel_ids > 0]).astype(np.int32)
    lookup = {int(parcel_id): col for col, parcel_id in enumerate(parcel_values)}
    gray = np.where(parcel_ids > 0)[0]
    cols = np.asarray([lookup[int(parcel_ids[g])] for g in gray], dtype=np.int32)
    counts = np.bincount(cols, minlength=parcel_values.size).astype(np.float64)
    weights = 1.0 / counts[cols]
    membership = sp.csr_matrix(
        (weights, (gray, cols)),
        shape=(parcel_ids.size, parcel_values.size),
        dtype=np.float64,
    )
    return membership, parcel_values


def l2_normalize_time_columns(x: np.ndarray) -> np.ndarray:
    out = np.asarray(x, dtype=np.float64)
    out = out - out.mean(axis=0, keepdims=True)
    norm = np.sqrt((out * out).sum(axis=0, keepdims=True))
    return np.divide(out, norm, out=np.zeros_like(out), where=norm > 1e-12)


def compact_timeseries(
    cifti_img: nib.Cifti2Image,
    run_slices: list[slice],
    reference_grayordinates: np.ndarray,
    membership: sp.csr_matrix,
) -> tuple[np.ndarray, np.ndarray]:
    reference_runs = []
    parcel_runs = []
    for run_index, run_slice in enumerate(run_slices, start=1):
        data = np.asanyarray(cifti_img.dataobj[run_slice, :]).astype(np.float32, copy=False)
        reference_runs.append(np.asarray(data[:, reference_grayordinates], dtype=np.float32))
        parcel_ts = np.asarray((membership.T @ data.T).T)
        parcel_runs.append(parcel_ts.astype(np.float32))
        print(f"[pfm-cv] extracted run {run_index}: timepoints={data.shape[0]}")
    return np.vstack(reference_runs), np.vstack(parcel_runs)


def fc_fingerprints(
    reference_ts: np.ndarray,
    parcel_ts: np.ndarray,
    indices: np.ndarray,
) -> np.ndarray:
    ref = l2_normalize_time_columns(reference_ts[indices])
    parcel = l2_normalize_time_columns(parcel_ts[indices])
    return np.clip(parcel.T @ ref, -0.999999, 0.999999).astype(np.float32)


def random_third_folds(
    run_slices: list[slice],
    repeats: int,
    block_length: int,
    seed: int,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(seed)
    folds = []
    for repeat in range(1, repeats + 1):
        thirds: list[list[np.ndarray]] = [[], [], []]
        for run_slice in run_slices:
            run_indices = np.arange(run_slice.start, run_slice.stop, dtype=np.int64)
            blocks = [run_indices[i : i + block_length] for i in range(0, run_indices.size, block_length)]
            order = rng.permutation(len(blocks))
            offset = int(rng.integers(0, 3))
            for position, block_index in enumerate(order):
                thirds[(position + offset) % 3].append(blocks[int(block_index)])
        third_indices = [np.sort(np.concatenate(group)) for group in thirds]
        for heldout in range(3):
            test = third_indices[heldout]
            train = np.sort(np.concatenate([third_indices[x] for x in range(3) if x != heldout]))
            folds.append(
                {
                    "repeat": repeat,
                    "heldout_third": heldout + 1,
                    "train_indices": train,
                    "test_indices": test,
                }
            )
    return folds


def normalized_pattern(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=np.float64)
    x = x - x.mean()
    norm = float(np.sqrt(np.dot(x, x)))
    return x / norm if norm > 1e-12 else np.zeros_like(x)


def convex_two_fit(y: np.ndarray, t0: np.ndarray, t1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    direction = t1 - t0
    denom = float(np.dot(direction, direction))
    w1 = 0.0 if denom <= 1e-12 else float(np.clip(np.dot(y - t0, direction) / denom, 0.0, 1.0))
    weights = np.array([1.0 - w1, w1], dtype=np.float64)
    return weights, (weights[0] * t0) + (weights[1] * t1)


def convex_three_fit(
    y: np.ndarray,
    t0: np.ndarray,
    t1: np.ndarray,
    t2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    templates = np.vstack((t0, t1, t2))
    candidates = []

    design = np.column_stack((t1 - t0, t2 - t0))
    try:
        w12 = np.linalg.solve(design.T @ design, design.T @ (y - t0))
        w = np.array([1.0 - w12.sum(), w12[0], w12[1]], dtype=np.float64)
        if np.all(w >= 0.0):
            candidates.append(w)
    except np.linalg.LinAlgError:
        pass

    for edge in ((0, 1), (0, 2), (1, 2)):
        weights2, _ = convex_two_fit(y, templates[edge[0]], templates[edge[1]])
        w = np.zeros(3, dtype=np.float64)
        w[edge[0]], w[edge[1]] = weights2
        candidates.append(w)

    best_weights = min(candidates, key=lambda w: float(np.sum((y - (w @ templates)) ** 2)))
    return best_weights, best_weights @ templates


def cosine_score(y: np.ndarray, prediction: np.ndarray) -> float:
    pred = normalized_pattern(prediction)
    return float(np.dot(y, pred)) if np.any(pred) else 0.0


def select_and_score_models(
    train: np.ndarray,
    test: np.ndarray,
    templates: np.ndarray,
    assigned_network: int,
) -> dict[str, object]:
    y_train = normalized_pattern(train)
    y_test = normalized_pattern(test)
    template_patterns = np.vstack([normalized_pattern(x) for x in templates])
    t0 = template_patterns[assigned_network]
    score1_train = cosine_score(y_train, t0)
    score1_test = cosine_score(y_test, t0)

    best2 = None
    candidate2_gains = np.full(template_patterns.shape[0], np.nan, dtype=np.float64)
    for alt in range(template_patterns.shape[0]):
        if alt == assigned_network:
            continue
        weights, prediction = convex_two_fit(y_train, t0, template_patterns[alt])
        score = cosine_score(y_train, prediction)
        score_test = cosine_score(y_test, prediction)
        candidate2_gains[alt] = score_test - score1_test
        candidate = (score, alt, weights, prediction, score_test)
        if best2 is None or candidate[0] > best2[0]:
            best2 = candidate
    assert best2 is not None
    _, alt1, weights2, prediction2, score2_test = best2

    best3 = None
    candidate3_gains = np.full(template_patterns.shape[0], np.nan, dtype=np.float64)
    for alt2 in range(template_patterns.shape[0]):
        if alt2 in (assigned_network, alt1):
            continue
        weights, prediction = convex_three_fit(
            y_train,
            t0,
            template_patterns[alt1],
            template_patterns[alt2],
        )
        score = cosine_score(y_train, prediction)
        score_test = cosine_score(y_test, prediction)
        candidate3_gains[alt2] = score_test - score2_test
        candidate = (score, alt2, weights, prediction, score_test)
        if best3 is None or candidate[0] > best3[0]:
            best3 = candidate
    assert best3 is not None
    _, alt2, weights3, prediction3, score3_test = best3

    return {
        "network_1": assigned_network,
        "network_2": int(alt1),
        "network_3": int(alt2),
        "weights_2": weights2,
        "weights_3": weights3,
        "score_1_train": score1_train,
        "score_1_test": score1_test,
        "score_2_test": score2_test,
        "score_3_test": score3_test,
        "delta_2_vs_1": score2_test - score1_test,
        "delta_3_vs_2": score3_test - score2_test,
        "candidate_2_test_gains": candidate2_gains,
        "candidate_3_test_gains": candidate3_gains,
    }


def repeat_sign_flip_max_null(
    gains: np.ndarray,
    fold_repeats: np.ndarray,
    domains: np.ndarray,
    alpha: float,
) -> tuple[dict[str, float], np.ndarray, int]:
    """Exact repeat-level sign-flip max-statistic null, stratified by domain."""
    repeat_values = np.unique(fold_repeats)
    if repeat_values.size > 15:
        raise ValueError("Exact sign-flip enumeration is limited to 15 repeats")
    n_patterns = 1 << repeat_values.size
    null_medians = np.empty((n_patterns, gains.shape[0]), dtype=np.float64)
    for pattern in range(n_patterns):
        signs = {
            int(repeat): (1.0 if pattern & (1 << bit) else -1.0)
            for bit, repeat in enumerate(repeat_values)
        }
        fold_sign = np.asarray([signs[int(repeat)] for repeat in fold_repeats], dtype=np.float64)
        null_medians[pattern] = np.median(gains * fold_sign[np.newaxis, :], axis=1)

    observed = np.median(gains, axis=1)
    thresholds = {}
    p_fwer = np.ones(gains.shape[0], dtype=np.float64)
    for domain in np.unique(domains):
        domain_mask = domains == domain
        max_distribution = np.max(null_medians[:, domain_mask], axis=1)
        ordered_max = np.sort(max_distribution)
        cutoff_index = min(
            ordered_max.size - 1,
            max(0, int(np.ceil((1.0 - alpha) * ordered_max.size)) - 1),
        )
        thresholds[str(domain)] = float(ordered_max[cutoff_index])
        p_fwer[domain_mask] = np.mean(
            max_distribution[:, np.newaxis] >= observed[domain_mask][np.newaxis, :],
            axis=0,
        )
    return thresholds, p_fwer, n_patterns


def benjamini_hochberg_by_domain(
    p_values: np.ndarray,
    domains: np.ndarray,
) -> np.ndarray:
    """Benjamini-Hochberg q-values computed separately by anatomical domain."""
    q_values = np.ones(p_values.size, dtype=np.float64)
    for domain in np.unique(domains):
        domain_indices = np.where(domains == domain)[0]
        order = domain_indices[np.argsort(p_values[domain_indices], kind="stable")]
        running = 1.0
        for rank_index in range(order.size - 1, -1, -1):
            parcel_index = int(order[rank_index])
            running = min(
                running,
                float(p_values[parcel_index]) * order.size / float(rank_index + 1),
            )
            q_values[parcel_index] = running
    return q_values


def candidate_identity_null(
    candidate_gains: np.ndarray,
    selected_networks: np.ndarray,
    excluded_networks: np.ndarray,
    fold_repeats: np.ndarray,
    iterations: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Conditional null for reproducible added-network identity.

    Within each random-third repeat, a cyclic permutation maps every selected
    added-network identity to another admissible candidate. This preserves the
    parcel's full held-out gain landscape—including generic gains from added
    model flexibility—while breaking the correspondence between the network
    selected on training data and its held-out identity. The test statistic is
    the median held-out gain across folds.
    """
    if candidate_gains.ndim != 3:
        raise ValueError("candidate_gains must have parcel x fold x network dimensions")
    n_parcels, n_folds, n_networks = candidate_gains.shape
    if selected_networks.shape != (n_parcels, n_folds):
        raise ValueError("selected_networks shape does not match candidate_gains")
    if excluded_networks.shape[:2] != (n_parcels, n_folds):
        raise ValueError("excluded_networks shape does not match candidate_gains")
    if fold_repeats.size != n_folds:
        raise ValueError("fold repeat labels do not match candidate_gains")
    if iterations < 1:
        raise ValueError("candidate-identity null iterations must be positive")

    n_excluded = excluded_networks.shape[2]
    n_candidates = n_networks - n_excluded
    network_ids = np.arange(n_networks, dtype=np.int32)
    valid = np.ones((n_parcels, n_folds, n_networks), dtype=bool)
    for exclusion_index in range(n_excluded):
        excluded = excluded_networks[:, :, exclusion_index]
        valid[
            np.arange(n_parcels)[:, np.newaxis],
            np.arange(n_folds)[np.newaxis, :],
            excluded,
        ] = False
    valid_ids = np.sort(
        np.where(valid, network_ids[np.newaxis, np.newaxis, :], n_networks),
        axis=2,
    )[:, :, :n_candidates]
    selected_positions = np.argmax(
        valid_ids == selected_networks[:, :, np.newaxis],
        axis=2,
    )
    if not np.all(
        np.take_along_axis(
            valid_ids,
            selected_positions[:, :, np.newaxis],
            axis=2,
        )[:, :, 0]
        == selected_networks
    ):
        raise ValueError("A selected network is not admissible under the supplied exclusions")

    observed_gains = np.take_along_axis(
        candidate_gains,
        selected_networks[:, :, np.newaxis],
        axis=2,
    )[:, :, 0]
    if not np.isfinite(observed_gains).all():
        raise ValueError("Observed selected-network gains contain non-finite values")
    observed_statistic = np.median(observed_gains, axis=1)

    rng = np.random.default_rng(seed)
    exceedances = np.zeros(n_parcels, dtype=np.int64)
    repeat_values = np.unique(fold_repeats)
    for _ in range(iterations):
        null_gains = np.empty((n_parcels, n_folds), dtype=np.float64)
        for repeat in repeat_values:
            fold_indices = np.where(fold_repeats == repeat)[0]
            offsets = rng.integers(0, n_candidates, size=n_parcels)
            mapped_positions = (
                selected_positions[:, fold_indices] + offsets[:, np.newaxis]
            ) % n_candidates
            mapped_networks = np.take_along_axis(
                valid_ids[:, fold_indices, :],
                mapped_positions[:, :, np.newaxis],
                axis=2,
            )[:, :, 0]
            null_gains[:, fold_indices] = np.take_along_axis(
                candidate_gains[:, fold_indices, :],
                mapped_networks[:, :, np.newaxis],
                axis=2,
            )[:, :, 0]
        null_statistic = np.median(null_gains, axis=1)
        exceedances += null_statistic >= observed_statistic

    p_values = (exceedances + 1.0) / float(iterations + 1)
    return p_values, observed_statistic


def minimum_repeat_mean_gain(gains: np.ndarray, fold_repeats: np.ndarray) -> np.ndarray:
    """Lowest mean held-out gain across random-third perturbation repeats."""
    repeat_means = np.column_stack(
        [np.mean(gains[:, fold_repeats == repeat], axis=1) for repeat in np.unique(fold_repeats)]
    )
    return np.min(repeat_means, axis=1)


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Fixed-parcel random-third PFM network-mixture evaluation")
    ap.add_argument("--in-cifti", required=True, type=Path)
    ap.add_argument("--fd-txt", required=True, type=Path)
    ap.add_argument("--raw-run-lengths", type=parse_int_csv)
    ap.add_argument(
        "--scan-info-json",
        type=Path,
        help="Concatenation ScanInfo.json; supplies run lengths when --raw-run-lengths is omitted",
    )
    ap.add_argument("--fd-threshold", type=float, default=0.3)
    ap.add_argument("--parcels-dlabel", required=True, type=Path)
    ap.add_argument("--priors-mat", required=True, type=Path)
    ap.add_argument("--distance-npy", required=True, type=Path)
    ap.add_argument("--reference-count", type=int, default=1500)
    ap.add_argument("--reference-seed", type=int, default=12345)
    ap.add_argument("--local-exclusion-mm", type=float, default=30.0)
    ap.add_argument("--split-repeats", type=int, default=8)
    ap.add_argument("--split-block-length", type=int, default=25)
    ap.add_argument("--split-seed", type=int, default=12345)
    ap.add_argument("--min-consistency-fraction", type=float, default=0.8)
    ap.add_argument("--min-positive-gain-fraction", type=float, default=0.8)
    ap.add_argument("--selection-correction", choices=["fdr", "fwer"], default="fdr")
    ap.add_argument("--fdr-alpha", type=float, default=0.05)
    ap.add_argument("--parcel-null-iterations", type=int, default=2000)
    ap.add_argument("--parcel-null-seed", type=int, default=12345)
    ap.add_argument("--min-second-network-weight", type=float, default=0.20)
    ap.add_argument("--min-third-network-weight", type=float, default=0.15)
    ap.add_argument("--null-alpha", type=float, default=0.05)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--outfile-prefix", default="ParcelRunSplitCV")
    return ap


def main() -> int:
    args = build_arg_parser().parse_args()
    if not 1 <= args.split_repeats <= 15 or args.split_block_length < 2:
        raise ValueError("split repeats must be in [1, 15] and block length must be at least 2")
    if not 0.0 < args.min_consistency_fraction <= 1.0:
        raise ValueError("--min-consistency-fraction must be in (0, 1]")
    if not 0.0 < args.min_positive_gain_fraction <= 1.0:
        raise ValueError("--min-positive-gain-fraction must be in (0, 1]")
    if not 0.0 < args.null_alpha < 1.0:
        raise ValueError("--null-alpha must be in (0, 1)")
    if not 0.0 < args.fdr_alpha < 1.0:
        raise ValueError("--fdr-alpha must be in (0, 1)")
    if args.parcel_null_iterations < 1:
        raise ValueError("--parcel-null-iterations must be positive")
    if not 0.0 <= args.min_second_network_weight <= 1.0:
        raise ValueError("--min-second-network-weight must be in [0, 1]")
    if not 0.0 <= args.min_third_network_weight <= 1.0:
        raise ValueError("--min-third-network-weight must be in [0, 1]")

    args.outdir.mkdir(parents=True, exist_ok=True)
    cifti_img = nib.load(str(args.in_cifti))
    parcel_img = nib.load(str(args.parcels_dlabel))
    if cifti_img.shape[1] != parcel_img.shape[1]:
        raise ValueError("Input and parcel CIFTIs have different grayordinate counts")
    parcel_axis = parcel_img.header.get_axis(0)
    if not isinstance(parcel_axis, LabelAxis):
        raise ValueError("--parcels-dlabel must be a dense label CIFTI")
    parcel_ids = np.asanyarray(parcel_img.dataobj).reshape(-1).astype(np.int32)
    brain_axis = cifti_img.header.get_axis(1)
    cortex_idx = cortical_grayordinates(brain_axis)
    priors_fc, network_names = load_network_priors(args.priors_mat, cortex_idx.size)

    if args.raw_run_lengths is not None:
        raw_run_lengths = args.raw_run_lengths
        run_length_source = "--raw-run-lengths"
    elif args.scan_info_json is not None:
        raw_run_lengths = load_raw_run_lengths(args.scan_info_json)
        run_length_source = str(args.scan_info_json.resolve())
    else:
        raise ValueError("Provide --scan-info-json or --raw-run-lengths")

    fd = np.loadtxt(args.fd_txt).reshape(-1)
    run_slices = censored_run_slices(fd, raw_run_lengths, args.fd_threshold)
    if run_slices[-1].stop != cifti_img.shape[0]:
        raise ValueError(
            f"Censored run counts end at {run_slices[-1].stop}, input has {cifti_img.shape[0]} timepoints"
        )

    membership, parcel_values = parcel_membership(parcel_ids)
    rng = np.random.default_rng(args.reference_seed)
    reference_positions = np.sort(
        rng.choice(cortex_idx.size, size=min(args.reference_count, cortex_idx.size), replace=False)
    )
    reference_gray = cortex_idx[reference_positions]
    reference_templates = priors_fc[reference_positions, :].T
    reference_ts, parcel_ts = compact_timeseries(
        cifti_img,
        run_slices,
        reference_gray,
        membership,
    )

    distance = np.load(args.distance_npy, mmap_mode="r")
    keep_masks = []
    for parcel_id in parcel_values:
        members = np.where(parcel_ids == int(parcel_id))[0]
        local_distance = np.asarray(distance[np.ix_(members, reference_gray)])
        keep = np.min(local_distance, axis=0) > args.local_exclusion_mm
        if int(keep.sum()) < 100:
            raise ValueError(f"Parcel {parcel_id} retains only {int(keep.sum())} FC references")
        keep_masks.append(keep)

    labels = parcel_axis.label[0]
    assigned_networks = []
    parcel_labels = []
    for parcel_id in parcel_values:
        label = labels[int(parcel_id)][0]
        assigned, _ = parcel_network_from_label(int(parcel_id), label, network_names)
        if assigned < 0:
            raise ValueError(f"Cannot map parcel label to a network: {label}")
        assigned_networks.append(assigned)
        parcel_labels.append(label)

    folds = random_third_folds(
        run_slices,
        args.split_repeats,
        args.split_block_length,
        args.split_seed,
    )
    n_folds = len(folds)
    min_consistent_folds = int(np.ceil(args.min_consistency_fraction * n_folds))
    min_positive_folds = int(np.ceil(args.min_positive_gain_fraction * n_folds))

    fold_rows = []
    parcel_fold_results: list[list[dict[str, object]]] = [[] for _ in parcel_values]
    for fold_index, fold in enumerate(folds, start=1):
        train_indices = fold["train_indices"]
        test_indices = fold["test_indices"]
        train_fc = fc_fingerprints(reference_ts, parcel_ts, train_indices)
        test_fc = fc_fingerprints(reference_ts, parcel_ts, test_indices)
        for parcel_col, parcel_id in enumerate(parcel_values):
            keep = keep_masks[parcel_col]
            result = select_and_score_models(
                train_fc[parcel_col, keep],
                test_fc[parcel_col, keep],
                reference_templates[:, keep],
                assigned_networks[parcel_col],
            )
            parcel_fold_results[parcel_col].append(result)
            fold_rows.append(
                {
                    "parcel_id": int(parcel_id),
                    "fold": fold_index,
                    "repeat": int(fold["repeat"]),
                    "heldout_third": int(fold["heldout_third"]),
                    "n_train_timepoints": int(train_indices.size),
                    "n_test_timepoints": int(test_indices.size),
                    "assigned_network_id": int(result["network_1"]) + 1,
                    "second_network_id": int(result["network_2"]) + 1,
                    "third_network_id": int(result["network_3"]) + 1,
                    "second_network": network_names[int(result["network_2"])],
                    "third_network": network_names[int(result["network_3"])],
                    "score_1_test": float(result["score_1_test"]),
                    "score_2_test": float(result["score_2_test"]),
                    "score_3_test": float(result["score_3_test"]),
                    "delta_2_vs_1": float(result["delta_2_vs_1"]),
                    "delta_3_vs_2": float(result["delta_3_vs_2"]),
                    "weights_2": ";".join(f"{x:.8g}" for x in result["weights_2"]),
                    "weights_3": ";".join(f"{x:.8g}" for x in result["weights_3"]),
                }
            )
        print(
            f"[pfm-cv] fold {fold_index}/{n_folds}: repeat={fold['repeat']} "
            f"heldout_third={fold['heldout_third']} train={train_indices.size} test={test_indices.size}"
        )

    cortex_mask = np.zeros(parcel_ids.size, dtype=bool)
    cortex_mask[cortex_idx] = True
    parcel_domains = []
    for parcel_id in parcel_values:
        parcel_members = parcel_ids == int(parcel_id)
        has_cortex = bool(np.any(cortex_mask & parcel_members))
        has_subcortex = bool(np.any(~cortex_mask & parcel_members))
        if has_cortex and has_subcortex:
            raise ValueError(f"Parcel {parcel_id} spans cortex and subcortex")
        parcel_domains.append("cortex" if has_cortex else "subcortex")
    parcel_domains = np.asarray(parcel_domains)

    fold_repeats = np.asarray([int(fold["repeat"]) for fold in folds], dtype=np.int32)
    delta2_matrix = np.asarray(
        [[float(result["delta_2_vs_1"]) for result in results] for results in parcel_fold_results]
    )
    delta3_matrix = np.asarray(
        [[float(result["delta_3_vs_2"]) for result in results] for results in parcel_fold_results]
    )
    null_threshold2, p_fwer2, null_patterns = repeat_sign_flip_max_null(
        delta2_matrix,
        fold_repeats,
        parcel_domains,
        args.null_alpha,
    )
    null_threshold3, p_fwer3, null_patterns3 = repeat_sign_flip_max_null(
        delta3_matrix,
        fold_repeats,
        parcel_domains,
        args.null_alpha,
    )
    if null_patterns3 != null_patterns:
        raise RuntimeError("Null enumerations unexpectedly differ")

    candidate2_matrix = np.asarray(
        [
            [np.asarray(result["candidate_2_test_gains"]) for result in results]
            for results in parcel_fold_results
        ]
    )
    candidate3_matrix = np.asarray(
        [
            [np.asarray(result["candidate_3_test_gains"]) for result in results]
            for results in parcel_fold_results
        ]
    )
    selected2_matrix = np.asarray(
        [[int(result["network_2"]) for result in results] for results in parcel_fold_results],
        dtype=np.int32,
    )
    selected3_matrix = np.asarray(
        [[int(result["network_3"]) for result in results] for results in parcel_fold_results],
        dtype=np.int32,
    )
    primary_matrix = np.broadcast_to(
        np.asarray(assigned_networks, dtype=np.int32)[:, np.newaxis],
        selected2_matrix.shape,
    )
    p_parcel2, parcel_null_statistic2 = candidate_identity_null(
        candidate2_matrix,
        selected2_matrix,
        primary_matrix[:, :, np.newaxis],
        fold_repeats,
        args.parcel_null_iterations,
        args.parcel_null_seed,
    )
    p_parcel3, parcel_null_statistic3 = candidate_identity_null(
        candidate3_matrix,
        selected3_matrix,
        np.stack((primary_matrix, selected2_matrix), axis=2),
        fold_repeats,
        args.parcel_null_iterations,
        args.parcel_null_seed + 1,
    )
    if not np.allclose(parcel_null_statistic2, np.median(delta2_matrix, axis=1)):
        raise RuntimeError("Two-network candidate-null statistic does not match observed gains")
    if not np.allclose(parcel_null_statistic3, np.median(delta3_matrix, axis=1)):
        raise RuntimeError("Three-network candidate-null statistic does not match observed gains")
    q_fdr2 = benjamini_hochberg_by_domain(p_parcel2, parcel_domains)
    q_fdr3 = benjamini_hochberg_by_domain(p_parcel3, parcel_domains)
    repeat_lower2 = minimum_repeat_mean_gain(delta2_matrix, fold_repeats)
    repeat_lower3 = minimum_repeat_mean_gain(delta3_matrix, fold_repeats)
    median_weight2 = np.asarray(
        [np.median([float(result["weights_2"][1]) for result in results]) for results in parcel_fold_results]
    )
    median_weight3 = np.asarray(
        [np.median([float(result["weights_3"][2]) for result in results]) for results in parcel_fold_results]
    )

    summary_rows = []
    scalar_maps = np.zeros((18, parcel_ids.size), dtype=np.float32)
    for parcel_col, parcel_id in enumerate(parcel_values):
        results = parcel_fold_results[parcel_col]
        domain = str(parcel_domains[parcel_col])
        second_ids = [int(x["network_2"]) for x in results]
        third_sets = [tuple(sorted((int(x["network_2"]), int(x["network_3"])))) for x in results]
        second_mode, second_consistency = Counter(second_ids).most_common(1)[0]
        third_mode, third_consistency = Counter(third_sets).most_common(1)[0]
        delta2 = np.asarray([float(x["delta_2_vs_1"]) for x in results])
        delta3 = np.asarray([float(x["delta_3_vs_2"]) for x in results])
        positive2 = int(np.sum(delta2 > 0.0))
        positive3 = int(np.sum(delta3 > 0.0))
        passes_null2 = bool(
            float(np.median(delta2)) > null_threshold2[domain]
            and float(p_fwer2[parcel_col]) <= args.null_alpha
        )
        passes_null3 = bool(
            float(np.median(delta3)) > null_threshold3[domain]
            and float(p_fwer3[parcel_col]) <= args.null_alpha
        )
        passes_fdr2 = bool(float(q_fdr2[parcel_col]) <= args.fdr_alpha)
        passes_fdr3 = bool(float(q_fdr3[parcel_col]) <= args.fdr_alpha)
        passes_repeat_gain2 = bool(float(repeat_lower2[parcel_col]) > 0.0)
        passes_repeat_gain3 = bool(float(repeat_lower3[parcel_col]) > 0.0)
        passes_weight2 = bool(
            float(median_weight2[parcel_col]) >= args.min_second_network_weight
        )
        passes_weight3 = bool(
            float(median_weight3[parcel_col]) >= args.min_third_network_weight
        )
        correction_pass2 = passes_fdr2 if args.selection_correction == "fdr" else passes_null2
        correction_pass3 = passes_fdr3 if args.selection_correction == "fdr" else passes_null3
        select2 = bool(
            second_consistency >= min_consistent_folds
            and positive2 >= min_positive_folds
            and correction_pass2
            and passes_repeat_gain2
            and passes_weight2
        )
        select3 = bool(
            select2
            and third_consistency >= min_consistent_folds
            and positive3 >= min_positive_folds
            and correction_pass3
            and passes_repeat_gain3
            and passes_weight3
            and second_mode in third_mode
        )
        fwer_select2 = bool(
            second_consistency >= min_consistent_folds
            and positive2 >= min_positive_folds
            and passes_null2
        )
        fwer_select3 = bool(
            fwer_select2
            and third_consistency >= min_consistent_folds
            and positive3 >= min_positive_folds
            and passes_null3
            and second_mode in third_mode
        )
        selected = [assigned_networks[parcel_col]]
        if select2:
            selected.append(second_mode)
        if select3:
            selected.append(next(x for x in third_mode if x != second_mode))
        row = {
            "parcel_id": int(parcel_id),
            "parcel_label": parcel_labels[parcel_col],
            "domain": domain,
            "n_grayordinates": int(np.sum(parcel_ids == int(parcel_id))),
            "assigned_network_id": assigned_networks[parcel_col] + 1,
            "assigned_network": network_names[assigned_networks[parcel_col]],
            "selected_network_count": len(selected),
            "selected_network_ids": ";".join(str(x + 1) for x in selected),
            "selected_networks": ";".join(network_names[x] for x in selected),
            "selection_correction": args.selection_correction,
            "fwer_high_confidence_network_count": 3 if fwer_select3 else 2 if fwer_select2 else 1,
            "second_network_mode_id": second_mode + 1,
            "second_network_mode": network_names[second_mode],
            "second_network_consistent_folds": second_consistency,
            "second_network_positive_gain_folds": positive2,
            "median_delta_2_vs_1": float(np.median(delta2)),
            "null_threshold_delta_2_vs_1": null_threshold2[domain],
            "fwer_p_delta_2_vs_1": float(p_fwer2[parcel_col]),
            "passes_null_2_vs_1": int(passes_null2),
            "candidate_identity_p_2_vs_1": float(p_parcel2[parcel_col]),
            "candidate_identity_fdr_q_2_vs_1": float(q_fdr2[parcel_col]),
            "passes_fdr_2_vs_1": int(passes_fdr2),
            "minimum_repeat_mean_delta_2_vs_1": float(repeat_lower2[parcel_col]),
            "passes_repeat_gain_2_vs_1": int(passes_repeat_gain2),
            "median_second_network_weight": float(median_weight2[parcel_col]),
            "passes_second_network_weight": int(passes_weight2),
            "third_network_set_mode_ids": ";".join(str(x + 1) for x in third_mode),
            "third_network_set_mode": ";".join(network_names[x] for x in third_mode),
            "third_network_set_consistent_folds": third_consistency,
            "third_network_positive_gain_folds": positive3,
            "median_delta_3_vs_2": float(np.median(delta3)),
            "null_threshold_delta_3_vs_2": null_threshold3[domain],
            "fwer_p_delta_3_vs_2": float(p_fwer3[parcel_col]),
            "passes_null_3_vs_2": int(passes_null3),
            "candidate_identity_p_3_vs_2": float(p_parcel3[parcel_col]),
            "candidate_identity_fdr_q_3_vs_2": float(q_fdr3[parcel_col]),
            "passes_fdr_3_vs_2": int(passes_fdr3),
            "minimum_repeat_mean_delta_3_vs_2": float(repeat_lower3[parcel_col]),
            "passes_repeat_gain_3_vs_2": int(passes_repeat_gain3),
            "median_third_network_weight": float(median_weight3[parcel_col]),
            "passes_third_network_weight": int(passes_weight3),
        }
        summary_rows.append(row)
        mask = parcel_ids == int(parcel_id)
        scalar_maps[0, mask] = len(selected)
        scalar_maps[1, mask] = float(np.median(delta2))
        scalar_maps[2, mask] = positive2 / float(n_folds)
        scalar_maps[3, mask] = second_consistency / float(n_folds)
        scalar_maps[4, mask] = float(np.median(delta3))
        scalar_maps[5, mask] = positive3 / float(n_folds)
        scalar_maps[6, mask] = third_consistency / float(n_folds)
        scalar_maps[7, mask] = float(np.median([x["score_1_test"] for x in results]))
        scalar_maps[8, mask] = float(p_fwer2[parcel_col])
        scalar_maps[9, mask] = float(p_fwer3[parcel_col])
        scalar_maps[10, mask] = float(p_parcel2[parcel_col])
        scalar_maps[11, mask] = float(q_fdr2[parcel_col])
        scalar_maps[12, mask] = float(repeat_lower2[parcel_col])
        scalar_maps[13, mask] = float(median_weight2[parcel_col])
        scalar_maps[14, mask] = float(p_parcel3[parcel_col])
        scalar_maps[15, mask] = float(q_fdr3[parcel_col])
        scalar_maps[16, mask] = float(repeat_lower3[parcel_col])
        scalar_maps[17, mask] = float(median_weight3[parcel_col])

    prefix = args.outfile_prefix
    summary_csv = args.outdir / f"{prefix}_ParcelSummary.csv"
    fold_csv = args.outdir / f"{prefix}_FoldResults.csv"
    dscalar = args.outdir / f"{prefix}_Diagnostics.dscalar.nii"
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    with fold_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fold_rows[0]))
        writer.writeheader()
        writer.writerows(fold_rows)
    save_dscalar(
        parcel_img,
        scalar_maps,
        [
            "CV selected network count",
            "Median held-out gain: 2 vs 1",
            "Positive-gain fold fraction: 2 vs 1",
            "Second-network identity consistency",
            "Median held-out gain: 3 vs 2",
            "Positive-gain fold fraction: 3 vs 2",
            "Three-network set consistency",
            "Median held-out one-network score",
            "Max-statistic FWER p: 2 vs 1",
            "Max-statistic FWER p: 3 vs 2",
            "Candidate-identity null p: 2 vs 1",
            "Candidate-identity FDR q: 2 vs 1",
            "Minimum repeat-mean gain: 2 vs 1",
            "Median second-network weight",
            "Candidate-identity null p: 3 vs 2",
            "Candidate-identity FDR q: 3 vs 2",
            "Minimum repeat-mean gain: 3 vs 2",
            "Median third-network weight",
        ],
        dscalar,
    )
    counts = Counter(int(row["selected_network_count"]) for row in summary_rows)
    fwer_counts = Counter(
        int(row["fwer_high_confidence_network_count"]) for row in summary_rows
    )
    summary = {
        "method": "fixed parcels; repeated blocked random thirds; train on two thirds and test on one",
        "input": str(args.in_cifti.resolve()),
        "raw_run_lengths": raw_run_lengths,
        "run_length_source": run_length_source,
        "run_timepoints_after_censoring": [slc.stop - slc.start for slc in run_slices],
        "reference_count": int(reference_positions.size),
        "local_exclusion_mm": args.local_exclusion_mm,
        "split_repeats": args.split_repeats,
        "split_block_length_timepoints": args.split_block_length,
        "split_seed": args.split_seed,
        "n_folds": n_folds,
        "selection_rule": {
            "active_correction": args.selection_correction,
            "minimum_identity_consistency_fraction": args.min_consistency_fraction,
            "minimum_positive_gain_fraction": args.min_positive_gain_fraction,
            "minimum_identity_consistent_folds": min_consistent_folds,
            "minimum_positive_gain_folds": min_positive_folds,
            "require_positive_mean_gain_in_every_randomization_repeat": True,
            "minimum_second_network_median_weight": args.min_second_network_weight,
            "minimum_third_network_median_weight": args.min_third_network_weight,
            "parcelwise_null": {
                "conditional_null_calibrated": True,
                "method": (
                    "repeat-preserving cyclic permutation of added-network identity within "
                    "each parcel's complete held-out candidate-gain landscape"
                ),
                "interpretation": (
                    "tests whether the network identity selected on training data generalizes "
                    "better than an arbitrary admissible added network while preserving generic "
                    "gains from added model flexibility"
                ),
                "iterations": args.parcel_null_iterations,
                "seed_2_vs_1": args.parcel_null_seed,
                "seed_3_vs_2": args.parcel_null_seed + 1,
                "fdr_method": "Benjamini-Hochberg separately within cortex and subcortex",
                "fdr_alpha": args.fdr_alpha,
            },
            "high_confidence_fwer_diagnostic": {
                "minimum_median_heldout_gain_by_domain": {
                    "2_vs_1": null_threshold2,
                    "3_vs_2": null_threshold3,
                },
                "null_method": (
                    "exact repeat-level sign flips of held-out gains; one-sided max statistic "
                    "separately within cortex and subcortex"
                ),
                "null_alpha": args.null_alpha,
                "null_sign_patterns": null_patterns,
            },
        },
        "parcel_counts_by_selected_network_count": {str(k): v for k, v in sorted(counts.items())},
        "fwer_high_confidence_counts_by_network_count": {
            str(k): v for k, v in sorted(fwer_counts.items())
        },
        "outputs": {
            "parcel_summary_csv": str(summary_csv.resolve()),
            "fold_results_csv": str(fold_csv.resolve()),
            "diagnostics_dscalar": str(dscalar.resolve()),
        },
    }
    with (args.outdir / f"{prefix}_Summary.json").open("w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
    print(
        f"[pfm-cv] selected network counts ({args.selection_correction}): "
        f"{dict(sorted(counts.items()))}"
    )
    print(f"[pfm-cv] FWER high-confidence counts: {dict(sorted(fwer_counts.items()))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
