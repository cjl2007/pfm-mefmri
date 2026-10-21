#!/usr/bin/env python3
"""Pure Python RidgePFM mapping for PFM."""

from __future__ import annotations

import argparse
import json
import hashlib
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import scipy.io as sio
from nibabel.cifti2.cifti2_axes import SeriesAxis


def zscore_rows(x: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    mu = x.mean(axis=1, keepdims=True)
    sd = x.std(axis=1, keepdims=True)
    sd = np.where(sd > eps, sd, 1.0)
    return (x - mu) / sd


def bool_int(value: str) -> int:
    ivalue = int(value)
    if ivalue not in (0, 1):
        raise argparse.ArgumentTypeError("expected 0 or 1")
    return ivalue


def clamp_prob_smooth(p: np.ndarray, floor: float = 1e-6, alpha: float = 1e-5) -> np.ndarray:
    """Match the legacy MATLAB prior safety transform."""
    out = np.asarray(p, dtype=np.float32)
    if alpha > 0:
        out = (1.0 - float(alpha)) * out + float(alpha) / out.shape[1]
    return np.clip(out, float(floor), 1.0 - float(floor))


def spatial_prior_term(p: np.ndarray, transform: str) -> np.ndarray:
    if transform == "log":
        return np.log(np.maximum(p, 1e-6))
    if transform == "logit":
        q = clamp_prob_smooth(p)
        return np.log(q / (1.0 - q))
    raise ValueError(f"Unknown spatial prior transform: {transform}")


def save_scalar_like(ref_img: nib.Cifti2Image, values: np.ndarray, out_path: Path) -> None:
    out = values.reshape(-1, 1).astype(np.float32, copy=False)
    axes = [ref_img.header.get_axis(i) for i in range(ref_img.ndim)]
    series = SeriesAxis(start=0.0, step=1.0, size=1)
    hdr = nib.Cifti2Header.from_axes((series, axes[1]))
    nib.save(nib.Cifti2Image(out.T, hdr, nifti_header=ref_img.nifti_header), str(out_path))


def save_prob_dtseries(ref_img: nib.Cifti2Image, prob: np.ndarray, out_path: Path) -> None:
    # prob: n_gray x n_net -> dense timeseries laid out as maps x grayordinates
    axes = [ref_img.header.get_axis(i) for i in range(ref_img.ndim)]
    series = SeriesAxis(start=0.0, step=1.0, size=prob.shape[1])
    hdr = nib.Cifti2Header.from_axes((series, axes[1]))
    nib.save(nib.Cifti2Image(prob.T.astype(np.float32), hdr, nifti_header=ref_img.nifti_header), str(out_path))


def write_label_list(out_txt: Path, labels: list[str], colors: np.ndarray) -> None:
    with out_txt.open("w") as f:
        for i, name in enumerate(labels, start=1):
            rgb = np.clip(np.rint(colors[i - 1] * 255.0), 0, 255).astype(int)
            f.write(f"{name}\n")
            f.write(f"{i} {rgb[0]} {rgb[1]} {rgb[2]} 255\n")


def parse_structures(axis, structures_csv: str) -> np.ndarray:
    requested = {s.strip().upper() for s in structures_csv.split(",") if s.strip()}
    keep = np.zeros((axis.size,), dtype=bool)
    for name, slc, _ in axis.iter_structures():
        stop = slc.stop if slc.stop is not None else axis.size
        tag = name.replace("CIFTI_STRUCTURE_", "").upper()
        if not requested or tag in requested:
            keep[slc.start:stop] = True
    return keep


def inject_subcortical_spatial_priors(
    axis,
    spatial_full: np.ndarray,
    subcort_priors_nii: str,
    n_net: int,
) -> np.ndarray:
    informative = np.zeros((axis.size,), dtype=bool)
    if not subcort_priors_nii:
        return informative
    vol = nib.load(subcort_priors_nii)
    if axis.affine is None or not np.allclose(vol.affine, axis.affine, atol=1e-4):
        raise ValueError("Subcortical priors affine must match the CIFTI volume grid")
    arr = vol.get_fdata(dtype=np.float32)
    if arr.ndim == 3:
        arr = arr[..., np.newaxis]
    if arr.ndim != 4:
        raise ValueError(f"Expected 3D/4D NIfTI for subcortical priors: {subcort_priors_nii}, got {arr.shape}")
    if arr.shape[3] < n_net:
        raise ValueError(
            f"Subcortical priors maps ({arr.shape[3]}) must be >= network count ({n_net}): {subcort_priors_nii}"
        )

    used = 0
    for name, slc, bm in axis.iter_structures():
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            continue
        vox = np.asarray(bm.voxel, dtype=np.int64)
        if vox.size == 0:
            continue
        if (
            (vox[:, 0] < 0).any()
            or (vox[:, 1] < 0).any()
            or (vox[:, 2] < 0).any()
            or (vox[:, 0] >= arr.shape[0]).any()
            or (vox[:, 1] >= arr.shape[1]).any()
            or (vox[:, 2] >= arr.shape[2]).any()
        ):
            raise ValueError(f"Subcortical voxel indices are out of bounds for priors volume: {subcort_priors_nii}")
        block = arr[vox[:, 0], vox[:, 1], vox[:, 2], :n_net]
        spatial_full[slc, :] = block
        informative[slc] = np.ptp(block, axis=1) > 1e-12
        used += int(vox.shape[0])

    print(f"[ridge] loaded ACPC subcortical priors: {subcort_priors_nii} (rows={used})")
    return informative


def target_specific_ridge(
    design: np.ndarray,
    fingerprints: np.ndarray,
    retained: list[np.ndarray],
    lam: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit each target with the design left after its local edge exclusion.

    ``fingerprints`` is feature x target and must already contain zeros for
    excluded features. ``retained`` gives the matching feature mask per target.
    The returned fit R2 is calculated only across retained FC features.
    """
    n_targets = fingerprints.shape[1]
    n_networks = design.shape[1]
    beta = np.zeros((n_targets, n_networks), dtype=np.float64)
    r2 = np.zeros((n_targets,), dtype=np.float32)
    eye = np.eye(n_networks, dtype=np.float64)
    for target, keep in enumerate(retained):
        if not np.any(keep):
            continue
        ak = design[keep, :]
        y = fingerprints[keep, target]
        gram = ak.T @ ak
        projection = ak.T @ y
        coef = np.linalg.solve(gram + lam * eye, projection)
        beta[target, :] = coef
        residual = y - ak @ coef
        ss_res = float(residual @ residual)
        centered = y - y.mean()
        ss_tot = float(centered @ centered)
        r2[target] = np.float32(1.0 - ss_res / max(ss_tot, 1e-8))
    return beta, r2


def solve_target_grams(
    projections: np.ndarray, grams: np.ndarray, lam: float
) -> np.ndarray:
    """Solve target-specific ridge systems from sufficient statistics."""
    eye = np.eye(projections.shape[1], dtype=np.float64)
    return np.linalg.solve(grams + lam * eye, projections[..., np.newaxis])[..., 0]


def main() -> int:
    ap = argparse.ArgumentParser(description="Python RidgePFM")
    ap.add_argument("--in-cifti", required=True)
    ap.add_argument("--distance-npy", required=True, help="npy uint8 distance matrix")
    ap.add_argument("--priors-mat", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--outfile", default="RidgePFM_VTX")
    ap.add_argument("--fc-weight", type=float, default=1.0)
    ap.add_argument("--fc-demean", type=bool_int, default=0, help="demean each target FC fingerprint after local edge exclusion")
    ap.add_argument("--positive-fc-only", type=bool_int, default=0, help="fit ridge model using only positive target FC edges")
    ap.add_argument("--spatial-weight", type=float, default=0.05)
    ap.add_argument(
        "--spatial-prior-transform",
        choices=("logit", "log"),
        default="logit",
        help="transform used for spatial prior fusion; logit matches legacy V4/V5 MATLAB RidgePFM",
    )
    ap.add_argument("--lambda", dest="lam", type=float, default=1.0)
    ap.add_argument("--local-exclusion-mm", type=float, default=10.0)
    ap.add_argument(
        "--fit-r2-mode",
        choices=("exact", "skip"),
        default="exact",
        help="exact computes full FC-profile fit; skip retains exact labels and writes NaN fit values",
    )
    ap.add_argument("--brain-structures-csv", default="")
    ap.add_argument("--subcort-priors-nii", default="", help="ACPC-space 4D subcortical priors NIfTI")
    ap.add_argument("--chunk", type=int, default=256)
    ap.add_argument("--frame-indices", help="NPY row indices for an explicitly aligned session/split")
    ap.add_argument("--fc-moments", help="Validated retained-feature FC energy/row-sum cache (NPZ), bound to input path/mtime/size and exclusion radius")
    ap.add_argument("--local-kernel", choices=('auto','numpy'), default='auto', help='Optional compiled exact local corrections; auto falls back to NumPy without a compiler')
    ap.add_argument("--left-surf", required=True)
    ap.add_argument("--right-surf", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    img = nib.load(args.in_cifti)
    data_tg = np.asanyarray(img.dataobj).astype(np.float32, copy=False)  # time x gray
    if args.frame_indices:
        frames = np.load(args.frame_indices)
        if frames.ndim != 1 or frames.dtype.kind not in 'iu' or len(frames) < 3 or np.any(np.diff(frames) <= 0) or frames[0] < 0 or frames[-1] >= len(data_tg):
            raise ValueError("Invalid frame indices")
        data_tg = data_tg[frames]
    data = data_tg.T  # gray x time
    axis = img.header.get_axis(1)

    # Cortex indices
    cortex_idx = []
    for name, slc, _ in axis.iter_structures():
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            stop = slc.stop if slc.stop is not None else axis.size
            cortex_idx.extend(range(slc.start, stop))
    cortex_idx = np.asarray(cortex_idx, dtype=np.int32)
    n_gray, n_time = data.shape
    n_cort = cortex_idx.size
    cortex_pos = np.full((n_gray,), -1, dtype=np.int32)
    cortex_pos[cortex_idx] = np.arange(n_cort, dtype=np.int32)

    pri = sio.loadmat(args.priors_mat, squeeze_me=True, struct_as_record=False)["Priors"]
    pri_fc = np.asarray(pri.FC, dtype=np.float32)
    if pri_fc.shape[0] != n_cort:
        raise ValueError(f"Priors.FC rows ({pri_fc.shape[0]}) must match cortical count ({n_cort})")
    n_net = pri_fc.shape[1]
    pri_spatial = np.asarray(pri.Spatial, dtype=np.float32)
    if pri_spatial.shape[0] != n_cort:
        pri_spatial = pri_spatial[:n_cort, :]
    if pri_spatial.shape[0] != n_cort or pri_spatial.shape[1] < n_net:
        raise ValueError(
            f"Priors.Spatial shape {pri_spatial.shape} cannot supply "
            f"{n_cort} cortical rows x {n_net} networks"
        )
    pri_spatial = pri_spatial[:, :n_net]
    labels_raw = np.asarray(pri.NetworkLabels).ravel()
    labels = [str(x).strip() for x in labels_raw[:n_net]]
    colors = np.asarray(pri.NetworkColors, dtype=np.float32)[:n_net, :]

    # Selected grayordinates by structure
    selected = parse_structures(axis, args.brain_structures_csv)
    if not selected.any():
        selected[:] = True
    valid_signal = np.isfinite(data).all(axis=1) & (np.std(data, axis=1) > 1e-8)
    if not valid_signal[cortex_idx].all():
        raise ValueError("Cortical FC features contain nonfinite or constant time series")
    good_idx = np.where(selected & valid_signal)[0]

    # Prepare design
    x_cort = zscore_rows(data[cortex_idx, :].astype(np.float64))
    a = pri_fc.astype(np.float64)
    a = a - a.mean(axis=0, keepdims=True)
    a /= np.maximum(np.sqrt((a * a).sum(axis=0, keepdims=True)), 1e-8)

    # Distance matrix (for local exclusion on cortical targets only). Keep it
    # memory mapped and read cortical target rows sequentially; advanced slicing
    # the complete cortex would copy several gigabytes for standard CIFTIs.
    dist = np.load(args.distance_npy, mmap_mode="r")
    if dist.ndim != 2 or dist.shape[0] < n_gray or dist.shape[1] < n_gray:
        raise ValueError(f"Distance matrix {dist.shape} does not cover {n_gray} grayordinates")

    # Spatial priors full
    spatial_full = np.full((n_gray, n_net), 0.5, dtype=np.float32)
    spatial_full[cortex_idx, :] = pri_spatial
    spatial_informative = np.zeros((n_gray,), dtype=bool)
    spatial_informative[cortex_idx] = np.ptp(pri_spatial, axis=1) > 1e-12
    spatial_informative |= inject_subcortical_spatial_priors(
        axis, spatial_full, args.subcort_priors_nii, n_net
    )
    spatial_full = clamp_prob_smooth(spatial_full)

    label_idx = np.zeros((n_gray,), dtype=np.int32)
    prob_best = np.zeros((n_gray,), dtype=np.float32)
    r2 = np.full((n_gray,), np.nan, dtype=np.float32)
    prob_all = np.zeros((n_gray, n_net), dtype=np.float32)
    temporal_templates = x_cort.T @ a
    global_gram = a.T @ a
    moments = None
    if args.fc_moments:
        if args.frame_indices or args.fc_demean or args.positive_fc_only:
            raise ValueError('FC moments require unmodified full-input signed FC')
        moments = np.load(args.fc_moments)
        source = Path(args.in_cifti).resolve(); stat = source.stat()
        expected = dict(source=str(source), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                        exclusion=args.local_exclusion_mm,
                        distance=str(Path(args.distance_npy).resolve()), distance_mtime=Path(args.distance_npy).stat().st_mtime_ns)
        if json.loads(str(moments['binding'])) != expected:
            raise ValueError('FC moment cache does not match input/exclusion binding')
        if moments['energy'].shape != (n_gray,) or moments['row_sum'].shape != (n_gray,):
            raise ValueError('FC moment cache shape mismatch')

    if args.fit_r2_mode == "skip" and (
        int(args.fc_demean) == 1 or int(args.positive_fc_only) == 1
    ):
        raise ValueError("fit-r2-mode=skip is incompatible with FC demeaning or positive-FC-only")

    native_terms = None
    local_indices = None
    if args.local_kernel == 'auto' and (args.fit_r2_mode == 'skip' or moments is not None):
        from pfm_ridge_local import local_terms
        local_indices = [np.flatnonzero(np.asarray(dist[gi, cortex_idx]) <= args.local_exclusion_mm) for gi in cortex_idx]
        native_terms = local_terms(x_cort, a, local_indices, outdir/'native_cache')
        print('[ridge] local correction kernel:', 'compiled' if native_terms is not None else 'numpy', flush=True)

    for i in range(0, good_idx.size, args.chunk):
        j = min(i + args.chunk, good_idx.size)
        g = good_idx[i:j]
        y = zscore_rows(data[g, :].astype(np.float64))
        grams = np.broadcast_to(global_gram, (g.size, n_net, n_net)).copy()
        retained = []
        excluded = []
        for c, gi in enumerate(g):
            keep = np.ones((n_cort,), dtype=bool)
            pos = int(cortex_pos[gi])
            if pos >= 0:
                # Distance rows follow the target order and are read sequentially.
                if native_terms is not None:
                    local = local_indices[pos]
                    keep[local] = False
                    grams[c] -= native_terms[1][pos]
                else:
                    keep &= np.asarray(dist[gi, cortex_idx]) > args.local_exclusion_mm
                    local = np.flatnonzero(~keep)
                    grams[c] -= a[local].T @ a[local]
            else:
                local = np.empty((0,), dtype=np.int64)
            excluded.append(local)
            retained.append(keep)

        if args.fit_r2_mode == "exact" and moments is None:
            # Rows use population-SD normalization, so division by n_time gives
            # Pearson FC and matches the checked validation engine.
            m = (x_cort @ y.T) / max(n_time, 1)  # n_cort x chunk
            # Demean FC fingerprints after local exclusion; useful without GSR/MGTR.
            for c, keep in enumerate(retained):
                if int(args.fc_demean) == 1:
                    if keep.any():
                        m[:, c] -= float(m[keep, c].mean())
                    else:
                        m[:, c] = 0.0
                m[~keep, c] = 0.0

            if int(args.positive_fc_only) == 1:
                beta = np.zeros((g.size, n_net), dtype=np.float64)
                pos_retained = []
                for c in range(g.size):
                    keep_pos = retained[c] & (m[:, c] > 0)
                    pos_retained.append(keep_pos)
                    if not keep_pos.any():
                        continue
                    ak = a[keep_pos, :]
                    beta[c, :] = np.linalg.solve(
                        ak.T @ ak + args.lam * np.eye(n_net), ak.T @ m[keep_pos, c]
                    )
                _, r2_chunk = target_specific_ridge(a, m, pos_retained, args.lam)
            else:
                projection = m.T @ a
                beta = solve_target_grams(projection, grams, args.lam)
                ss = np.sum(m * m, axis=0)
                counts = np.asarray([np.sum(k) for k in retained])
                total = ss - np.sum(m, axis=0)**2 / np.maximum(counts, 1)
                residual = ss - 2*np.sum(beta*projection, axis=1) + np.einsum('bi,bij,bj->b', beta, grams, beta)
                r2_chunk = (1 - residual / np.maximum(total, 1e-8)).astype(np.float32)
                r2_chunk[counts == 0] = 0
            r2[g] = r2_chunk
        else:
            # Exact fast path for labels: project directly onto the small network
            # dictionary and subtract only locally excluded cortical edges.
            projection = (y @ temporal_templates) / max(n_time, 1)
            for c, local in enumerate(excluded):
                if native_terms is not None and cortex_pos[g[c]] >= 0:
                    projection[c] -= native_terms[0][cortex_pos[g[c]]]
                elif local.size:
                    local_fc = (x_cort[local] @ y[c]) / max(n_time, 1)
                    projection[c] -= local_fc @ a[local]
            beta = solve_target_grams(projection, grams, args.lam)
            if moments is not None:
                energy = moments['energy'][g]
                count = np.asarray([np.sum(k) for k in retained])
                total = energy - moments['row_sum'][g]**2 / np.maximum(count, 1)
                residual = energy - 2*np.sum(beta*projection, axis=1) + np.einsum('bi,bij,bj->b', beta, grams, beta)
                r2[g] = 1 - residual / np.maximum(total, 1e-8)
        beta_mu = beta.mean(axis=1, keepdims=True)
        beta_sd = np.maximum(beta.std(axis=1, keepdims=True), 1e-8)
        beta_z = (beta - beta_mu) / beta_sd

        score = args.fc_weight * beta_z + args.spatial_weight * spatial_prior_term(
            spatial_full[g, :], args.spatial_prior_transform
        )
        score -= score.max(axis=1, keepdims=True)
        p = np.exp(score)
        p /= np.maximum(p.sum(axis=1, keepdims=True), 1e-8)

        best = np.argmax(p, axis=1)
        label_idx[g] = best + 1
        prob_best[g] = p[np.arange(p.shape[0]), best]
        prob_all[g, :] = p

        if (i // args.chunk + 1) % 10 == 0:
            print(f"[ridge] processed {j}/{good_idx.size}")

    # Write outputs
    label_tmp = outdir / "Tmp_labels.dtseries.nii"
    save_scalar_like(img, label_idx.astype(np.float32), label_tmp)
    labfile = outdir / "LabelListFile.txt"
    write_label_list(labfile, labels, colors)
    dlabel = outdir / f"{args.outfile}.dlabel.nii"
    subprocess.run(
        ["wb_command", "-cifti-label-import", str(label_tmp), str(labfile), str(dlabel), "-discard-others"],
        check=True,
    )
    subprocess.run(
        ["wb_command", "-cifti-label-to-border", str(dlabel), "-border", args.left_surf, str(outdir / f"{args.outfile}.L.border")],
        check=True,
    )
    subprocess.run(
        ["wb_command", "-cifti-label-to-border", str(dlabel), "-border", args.right_surf, str(outdir / f"{args.outfile}.R.border")],
        check=True,
    )

    save_scalar_like(img, r2, outdir / f"{args.outfile}_TemplateFit_R2.dtseries.nii")
    save_scalar_like(
        img,
        prob_best,
        outdir / f"{args.outfile}_TopFusedScoreWeight.dtseries.nii",
    )
    save_scalar_like(
        img,
        spatial_informative.astype(np.float32),
        outdir / f"{args.outfile}_SpatialPriorAvailable.dtseries.nii",
    )
    # This softmax is a normalized fused-score distribution. It is useful for
    # ranking competing labels, but it is not calibrated P(true network).
    save_prob_dtseries(
        img,
        prob_all,
        outdir / f"{args.outfile}_FusedScoreSoftmax.dtseries.nii",
    )
    provenance = dict(parameters=vars(args), network_names=labels,
                      cortical_prior_sha256=hashlib.sha256(Path(args.priors_mat).read_bytes()).hexdigest(),
                      input_bytes=Path(args.in_cifti).stat().st_size,
                      input_mtime_ns=Path(args.in_cifti).stat().st_mtime_ns,
                      mapped_grayordinates=int(len(good_idx)),
                      local_kernel='compiled' if native_terms is not None else 'numpy',
                      interpretation='WTA and normalized fused scores; softmax weights are not calibrated assignment probabilities.')
    (outdir / f'{args.outfile}_Mapping.json').write_text(json.dumps(provenance, indent=2)+'\n')
    try:
        if label_tmp.exists():
            label_tmp.unlink()
        if labfile.exists():
            labfile.unlink()
    except Exception:
        pass
    print(f"[ridge] wrote outputs to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
