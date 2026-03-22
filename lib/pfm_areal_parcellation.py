#!/usr/bin/env python3
"""Python areal parcellation pass with MATLAB-like ridge-fusion post-processing.

Implements the same high-level flow used in pfm_areal_parcellation.m:
1) Build cortical adjacency from neighbor table.
2) Build far-only FC fingerprints from cortical time series + distance matrix.
3) Compute local FC-fingerprint gradient proxy.
4) Clean WTA via majority smoothing + small-component pruning.
5) Sub-parcellate each cleaned WTA patch (kmeans mode) with diffusion refinement.
6) Remove tiny parcels and conservatively merge weak seams.
7) Write dlabel + borders.
"""

from __future__ import annotations

import argparse
import subprocess
from collections import Counter, deque
from pathlib import Path
from typing import Iterable

import nibabel as nib
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
from nibabel.cifti2.cifti2_axes import SeriesAxis
from scipy.cluster.vq import kmeans2
from scipy.sparse.linalg import factorized


def save_scalar_like(ref_img: nib.Cifti2Image, values: np.ndarray, out_path: Path) -> None:
    out = values.reshape(-1, 1).astype(np.float32, copy=False)
    axes = [ref_img.header.get_axis(i) for i in range(ref_img.ndim)]
    series = SeriesAxis(start=0.0, step=1.0, size=1)
    hdr = nib.Cifti2Header.from_axes((series, axes[1]))
    nib.save(nib.Cifti2Image(out.T, hdr, nifti_header=ref_img.nifti_header), str(out_path))


def load_neighbors_cortex(neighbors_mat: str) -> np.ndarray:
    nmat = sio.loadmat(neighbors_mat, squeeze_me=True, struct_as_record=False)
    nkey = [k for k in nmat.keys() if not k.startswith("__")][0]
    arr = np.asarray(nmat[nkey], dtype=np.int64)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D neighbors matrix, got {arr.shape}")
    # Legacy mats often store 1-based self index in col 0.
    if arr.shape[1] >= 2 and np.all(arr[:, 0] == np.arange(1, arr.shape[0] + 1)):
        arr = arr[:, 1:]
    return arr


def cortical_grayordinates(axis) -> np.ndarray:
    idx = []
    for name, slc, _ in axis.iter_structures():
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            stop = slc.stop if slc.stop is not None else axis.size
            idx.extend(range(slc.start, stop))
    return np.asarray(idx, dtype=np.int64)


def build_full_graph(axis, neighbors_cortex: np.ndarray) -> sp.csr_matrix:
    n_gray = axis.size
    cortex_idx = cortical_grayordinates(axis)
    nc = cortex_idx.size
    if neighbors_cortex.shape[0] != nc:
        raise ValueError(
            f"Cortical neighbor rows ({neighbors_cortex.shape[0]}) must match cortical grayordinates ({nc})"
        )

    rows = []
    cols = []

    # Cortical edges from provided surface neighbor table.
    for r in range(nc):
        g_r = int(cortex_idx[r])
        for nb in neighbors_cortex[r]:
            nb = int(nb)
            if nb <= 0:
                continue
            j = nb - 1
            if 0 <= j < nc:
                rows.append(g_r)
                cols.append(int(cortex_idx[j]))

    # Subcortical edges via 6-neighborhood in voxel index space.
    for name, slc, bm in axis.iter_structures():
        if name in ("CIFTI_STRUCTURE_CORTEX_LEFT", "CIFTI_STRUCTURE_CORTEX_RIGHT"):
            continue
        stop = slc.stop if slc.stop is not None else n_gray
        gidx = np.arange(slc.start, stop, dtype=np.int64)
        vox = np.asarray(bm.voxel, dtype=np.int64)
        if vox.shape[0] != gidx.size:
            continue
        lut = {tuple(v.tolist()): int(g) for v, g in zip(vox, gidx)}
        for v, g in zip(vox, gidx):
            x, y, z = int(v[0]), int(v[1]), int(v[2])
            for dv in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                key = (x + dv[0], y + dv[1], z + dv[2])
                nb = lut.get(key)
                if nb is not None:
                    rows.append(int(g))
                    cols.append(int(nb))

    if not rows:
        raise ValueError("No edges built from neighbor table/brain models")
    data = np.ones((len(rows),), dtype=np.uint8)
    w = sp.csr_matrix((data, (rows, cols)), shape=(n_gray, n_gray), dtype=np.uint8)
    return ((w + w.T) > 0).astype(np.uint8).tocsr()


def row_l2_demean(x: np.ndarray) -> np.ndarray:
    x = x.astype(np.float64, copy=False)
    x = x - x.mean(axis=1, keepdims=True)
    d = np.sqrt((x * x).sum(axis=1, keepdims=True)) + 1e-8
    return x / d


def select_reference_vertices(dist_mm: np.ndarray, min_mm: float, ref_max: int, rngseed: int) -> np.ndarray:
    dist_mm = np.asarray(dist_mm, dtype=np.float64)
    if dist_mm.ndim != 2 or dist_mm.shape[0] != dist_mm.shape[1]:
        raise ValueError(f"Distance matrix must be square, got {dist_mm.shape}")
    n = dist_mm.shape[0]
    rng = np.random.default_rng(int(rngseed))

    ref = np.zeros((max(1, int(ref_max)),), dtype=np.int64)
    ref[0] = int(rng.integers(0, n))
    used = 1
    mindist = dist_mm[:, ref[0]].copy()
    mindist[~np.isfinite(mindist)] = -np.inf

    while used < ref.size:
        j = int(np.argmax(mindist))
        mx = float(mindist[j])
        if mx < min_mm:
            break
        ref[used] = j
        used += 1
        dj = dist_mm[:, j]
        dj = np.where(np.isfinite(dj), dj, -np.inf)
        mindist = np.minimum(mindist, dj)

    ref = np.unique(ref[:used])

    # Thinning pass to enforce spacing.
    keep = np.ones((ref.size,), dtype=bool)
    for i in range(ref.size):
        if not keep[i]:
            continue
        too_close = np.where(dist_mm[ref[i], ref] < min_mm)[0]
        too_close = too_close[too_close > i]
        keep[too_close] = False
    ref = ref[keep]
    return ref.astype(np.int64)


def compute_local_gradient_fcfp(w: sp.csr_matrix, f: np.ndarray) -> np.ndarray:
    deg = np.asarray(w.sum(axis=1)).ravel().astype(np.float64)
    deg[deg == 0] = 1.0
    nb = w @ f
    mu = (f * nb).sum(axis=1) / deg
    gmin = float(mu.min())
    gmax = float(mu.max())
    return ((mu - gmin) / (gmax - gmin + 1e-8)).astype(np.float32)


def mode_nonzero(vals: Iterable[int], fallback: int = 0) -> int:
    arr = [int(v) for v in vals if int(v) != 0]
    if not arr:
        return int(fallback)
    return int(Counter(arr).most_common(1)[0][0])


def adjacency_list(w: sp.csr_matrix) -> list[np.ndarray]:
    w = w.tocsr()
    out = []
    for i in range(w.shape[0]):
        st, en = w.indptr[i], w.indptr[i + 1]
        out.append(w.indices[st:en])
    return out


def clean_wta_majority(wta: np.ndarray, adj: list[np.ndarray], iters: int) -> np.ndarray:
    out = wta.astype(np.int32, copy=True)
    for _ in range(max(0, int(iters))):
        new = out.copy()
        for v in range(out.size):
            m = mode_nonzero(out[adj[v]], fallback=int(out[v]))
            if m != 0:
                new[v] = m
        out = new
    return out


def cc_labels_mask(mask: np.ndarray, adj: list[np.ndarray]) -> np.ndarray:
    n = mask.size
    cc = np.zeros((n,), dtype=np.int32)
    seen = np.zeros((n,), dtype=bool)
    cid = 0
    for v in np.where(mask)[0]:
        if seen[v]:
            continue
        cid += 1
        q = deque([int(v)])
        seen[v] = True
        cc[v] = cid
        while q:
            u = q.popleft()
            for nb in adj[u]:
                nb = int(nb)
                if (not seen[nb]) and mask[nb]:
                    seen[nb] = True
                    cc[nb] = cid
                    q.append(nb)
    return cc


def prune_small_wta_components(
    wta: np.ndarray,
    adj: list[np.ndarray],
    min_comp: int,
    max_pass: int,
    verbose: bool,
) -> np.ndarray:
    out = wta.astype(np.int32, copy=True)
    for p in range(max(1, int(max_pass))):
        changed = False
        for lab in np.unique(out):
            if lab <= 0:
                continue
            mask = out == lab
            if not mask.any():
                continue
            cc = cc_labels_mask(mask, adj)
            for cid in np.unique(cc[cc > 0]):
                verts = np.where(cc == cid)[0]
                if verts.size < int(min_comp):
                    for v in verts:
                        m = mode_nonzero([out[int(nb)] for nb in adj[int(v)] if out[int(nb)] not in (0, lab)], fallback=0)
                        if m > 0:
                            out[v] = m
                            changed = True
        if verbose:
            print(f"[wta] prune pass {p + 1}: changed={changed}")
        if not changed:
            break
    return out


def wta_connected_components(wta_clean: np.ndarray, adj: list[np.ndarray], active_mask: np.ndarray) -> tuple[np.ndarray, int]:
    n = wta_clean.size
    patch = np.zeros((n,), dtype=np.int32)
    visited = np.zeros((n,), dtype=bool)
    pid = 0
    for v in np.where((wta_clean > 0) & active_mask)[0]:
        if visited[v]:
            continue
        pid += 1
        lab = int(wta_clean[v])
        q = deque([int(v)])
        visited[v] = True
        patch[v] = pid
        while q:
            u = q.popleft()
            for nb in adj[u]:
                nb = int(nb)
                if (not visited[nb]) and active_mask[nb] and int(wta_clean[nb]) == lab:
                    visited[nb] = True
                    patch[nb] = pid
                    q.append(nb)
    return patch, pid


def suggest_k_bounds_per_patch(vidx: np.ndarray, va: np.ndarray | None, opts) -> tuple[int, int, int]:
    min_sz_verts = max(8, int(round(opts.min_parcel_size)))
    use_area = va is not None and np.isfinite(va[vidx]).all()
    if use_area:
        patch_area = float(np.sum(va[vidx]))
        k_area = max(1, int(np.floor(patch_area / max(1.0, float(opts.min_parcel_area_mm2)))))
        k_by_size = max(1, int(np.floor(vidx.size / min_sz_verts)))
        k_ceiling = min(int(opts.kmax_cap), max(1, min(k_area, k_by_size)))
    else:
        k_ceiling = min(int(opts.kmax_cap), max(1, int(np.floor(vidx.size / min_sz_verts))))
    k0 = max(1, k_ceiling)
    halfwin = int(max(0, opts.k_search_halfwin))
    kmin_p = max(int(opts.kmin_per_patch), k0 - halfwin)
    kmax_p = max(kmin_p, min(k_ceiling, k0 + halfwin))
    return int(k0), int(kmin_p), int(kmax_p)


def pca_reduce_rows(x: np.ndarray, n_comp: int) -> np.ndarray:
    x = x.astype(np.float64, copy=False)
    x = x - x.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    d = min(n_comp, u.shape[1])
    return (u[:, :d] * s[:d]).astype(np.float64)


def enforce_min_cluster_size(lbl: np.ndarray, xfeat: np.ndarray, min_sz: int) -> np.ndarray:
    lbl = lbl.astype(np.int32, copy=True).ravel()
    if lbl.size == 0:
        return lbl
    # make contiguous 1..K
    _, lbl = np.unique(lbl, return_inverse=True)
    lbl = lbl + 1
    k = int(lbl.max())
    if k <= 1:
        return lbl

    c = np.zeros((k, xfeat.shape[1]), dtype=np.float64)
    for i in range(1, k + 1):
        m = lbl == i
        if m.any():
            c[i - 1] = xfeat[m].mean(axis=0)

    for i in range(1, k + 1):
        idx = np.where(lbl == i)[0]
        if 0 < idx.size < int(min_sz):
            d = ((xfeat[idx, None, :] - c[None, :, :]) ** 2).sum(axis=2)
            d[:, i - 1] = np.inf
            to = np.argmin(d, axis=1) + 1
            lbl[idx] = to.astype(np.int32)

    _, lbl = np.unique(lbl, return_inverse=True)
    return (lbl + 1).astype(np.int32)


def split_kmeans_fcfp(vidx: np.ndarray, w: sp.csr_matrix, f: np.ndarray, k: int, opts) -> np.ndarray:
    n_tot = w.shape[0]
    ploc = np.zeros((n_tot,), dtype=np.int32)
    if vidx.size == 0:
        return ploc

    fi = f[vidx, :].astype(np.float64, copy=False)
    if fi.shape[1] > int(opts.kmeans_dim_fc_all):
        xk = pca_reduce_rows(fi, int(opts.kmeans_dim_fc_all))
    else:
        xk = fi

    min_sz = max(8, int(round(opts.min_parcel_size)))
    k = max(1, min(int(k), vidx.size))
    k = min(k, max(1, int(np.floor(vidx.size / min_sz))))

    cent, lbl0 = kmeans2(
        data=xk,
        k=k,
        iter=100,
        minit="++",
        missing="warn",
    )
    _ = cent
    lbl = enforce_min_cluster_size(lbl0.astype(np.int32) + 1, xk, min_sz)

    wloc = w[vidx[:, None], vidx].astype(np.float64).tocsr()
    deg = np.asarray(wloc.sum(axis=1)).ravel()
    l = sp.diags(deg, offsets=0) - wloc
    a = sp.eye(l.shape[0], dtype=np.float64, format="csr") + float(opts.lambda_spatial) * l
    a = a + 1e-6 * sp.eye(a.shape[0], dtype=np.float64, format="csr")
    solve = factorized(a.tocsc())

    prev = np.zeros_like(lbl)
    for _ in range(max(1, int(opts.diffuse_iters))):
        u = np.zeros((lbl.size, int(lbl.max())), dtype=np.float64)
        u[np.arange(lbl.size), lbl - 1] = 1.0
        u = solve(u)
        lbl = (np.argmax(u, axis=1) + 1).astype(np.int32)
        if np.array_equal(lbl, prev):
            break
        prev = lbl.copy()

    ploc[vidx] = lbl
    return ploc


def cluster_contrast_fcfp(lbl: np.ndarray, fi: np.ndarray) -> float:
    if lbl.size == 0:
        return -np.inf
    _, lblu = np.unique(lbl.astype(np.int32), return_inverse=True)
    lblu = lblu + 1
    k = int(lblu.max())
    if k < 2:
        return -np.inf
    c = np.zeros((k, fi.shape[1]), dtype=np.float64)
    for i in range(1, k + 1):
        c[i - 1] = fi[lblu == i].mean(axis=0)
    fi_n = fi / (np.sqrt((fi * fi).sum(axis=1, keepdims=True)) + 1e-8)
    c_n = c / (np.sqrt((c * c).sum(axis=1, keepdims=True)) + 1e-8)
    s = fi_n @ c_n.T
    own = s[np.arange(s.shape[0]), lblu - 1].copy()
    s[np.arange(s.shape[0]), lblu - 1] = -np.inf
    other = np.max(s, axis=1)
    return float(np.mean(own - other))


def auto_k_in_range_fcfp(vidx: np.ndarray, w: sp.csr_matrix, f: np.ndarray, k0: int, kmin_p: int, kmax_p: int, opts) -> int:
    best_k = int(kmin_p)
    best_s = -np.inf
    for k in range(int(kmin_p), int(kmax_p) + 1):
        lbl = split_kmeans_fcfp(vidx, w, f, int(k), opts)
        s = cluster_contrast_fcfp(lbl[vidx], f[vidx, :])
        if s > best_s:
            best_s = s
            best_k = int(k)
    if not np.isfinite(best_s):
        best_k = max(1, min(int(kmax_p), int(k0)))
    min_sz = max(8, int(round(opts.min_parcel_size)))
    k_max_by_size = max(1, int(np.floor(vidx.size / min_sz)))
    return max(1, min(best_k, k_max_by_size))


def parcel_neighbors(verts: np.ndarray, ploc: np.ndarray, patch_mask: np.ndarray, adj: list[np.ndarray]) -> np.ndarray:
    in_a = np.zeros((ploc.size,), dtype=bool)
    in_a[verts] = True
    nb = []
    for u in verts:
        for w in adj[int(u)]:
            w = int(w)
            if patch_mask[w] and (not in_a[w]) and ploc[w] > 0:
                nb.append(int(ploc[w]))
    if not nb:
        return np.empty((0,), dtype=np.int32)
    return np.unique(np.asarray(nb, dtype=np.int32))


def boundary_strength(averts: np.ndarray, bverts: np.ndarray, adj: list[np.ndarray], grad: np.ndarray) -> float:
    if averts.size == 0 or bverts.size == 0:
        return np.inf
    aset = set(int(v) for v in averts)
    bset = set(int(v) for v in bverts)
    vals = []
    for u in aset:
        for w in adj[u]:
            w = int(w)
            if w in bset:
                vals.append(0.5 * (float(grad[u]) + float(grad[w])))
    if not vals:
        return np.inf
    return float(np.mean(vals))


def remove_tiny_and_soft_merge(
    ploc: np.ndarray,
    patch_mask: np.ndarray,
    adj: list[np.ndarray],
    grad: np.ndarray,
    va: np.ndarray | None,
    opts,
) -> np.ndarray:
    out = ploc.astype(np.int32, copy=True)
    use_area = va is not None and va.size == out.size and np.isfinite(va).all()
    min_area = float(opts.min_parcel_area_mm2)
    min_sz = max(8, int(round(opts.min_parcel_size)))

    changed = True
    while changed:
        changed = False
        u = np.unique(out[patch_mask])
        u = u[u > 0]
        for uu in u:
            verts = np.where((out == uu) & patch_mask)[0]
            if verts.size == 0:
                continue
            tiny = float(np.sum(va[verts])) < min_area if use_area else verts.size < min_sz
            if tiny:
                nb = parcel_neighbors(verts, out, patch_mask, adj)
                if nb.size == 0:
                    continue
                best = int(nb[0])
                best_s = boundary_strength(verts, np.where((out == best) & patch_mask)[0], adj, grad)
                for cand in nb[1:]:
                    cand = int(cand)
                    s = boundary_strength(verts, np.where((out == cand) & patch_mask)[0], adj, grad)
                    if s < best_s:
                        best_s = s
                        best = cand
                out[verts] = best
                changed = True

    merge_again = True
    while merge_again:
        merge_again = False
        u = np.unique(out[patch_mask])
        u = u[u > 0]
        for lab in u:
            av = np.where((out == lab) & patch_mask)[0]
            nb = parcel_neighbors(av, out, patch_mask, adj)
            for cand in nb:
                cand = int(cand)
                s = boundary_strength(av, np.where((out == cand) & patch_mask)[0], adj, grad)
                if s < float(opts.merge_thresh):
                    out[(out == cand) & patch_mask] = int(lab)
                    merge_again = True
                    break
            if merge_again:
                break
    return out


def load_priors(priors_mat: str):
    if not priors_mat:
        return None
    try:
        pri = sio.loadmat(priors_mat, squeeze_me=True, struct_as_record=False)["Priors"]
    except Exception as exc:
        print(f"[areal] WARNING: failed to load priors mat {priors_mat}: {exc}")
        return None
    return pri


def write_parcel_label_list(labfile: Path, parcel_net: np.ndarray, pri) -> None:
    have_names = pri is not None and hasattr(pri, "NetworkLabels")
    have_colors = pri is not None and hasattr(pri, "NetworkColors")
    labels = np.asarray(pri.NetworkLabels).ravel() if have_names else None
    colors = np.asarray(pri.NetworkColors, dtype=np.float64) if have_colors else None

    with labfile.open("w") as f:
        for i, net_id in enumerate(parcel_net.tolist(), start=1):
            net_id = int(net_id)
            if have_names and 1 <= net_id <= labels.size:
                base = str(labels[net_id - 1]).strip()
            else:
                base = f"Net_{net_id}"
            if have_colors and 1 <= net_id <= colors.shape[0]:
                rgb = np.clip(np.rint(255.0 * colors[net_id - 1, :3]), 0, 255).astype(int)
            else:
                rgb = np.array([200, 200, 200], dtype=int)
            f.write(f"{base}_{i}\n")
            f.write(f"{i} {rgb[0]} {rgb[1]} {rgb[2]} 255\n")


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="MATLAB-style areal parcellation from ridge WTA labels")
    ap.add_argument("--in-cifti", required=True)
    ap.add_argument("--wta-dlabel", required=True)
    ap.add_argument("--neighbors-mat", required=True)
    ap.add_argument("--distance-npy", required=True)
    ap.add_argument("--priors-mat", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--outfile", default="RidgeFusion_VTX+ArealParcellation")
    ap.add_argument("--left-surf", required=True)
    ap.add_argument("--right-surf", required=True)

    ap.add_argument("--method", default="kmeans", choices=["kmeans"])
    ap.add_argument("--kmeans-dim-fc-all", type=int, default=64)
    ap.add_argument("--kmax-cap", type=int, default=8)
    ap.add_argument("--kmin-per-patch", type=int, default=1)
    ap.add_argument("--k-search-halfwin", type=int, default=1)
    ap.add_argument("--k-fixed", type=int, default=0)
    ap.add_argument("--kmeans-replicates", type=int, default=3)

    ap.add_argument("--min-parcel-area-mm2", type=float, default=30.0)
    ap.add_argument("--min-parcel-size", type=int, default=30)
    ap.add_argument("--diffuse-iters", type=int, default=6)
    ap.add_argument("--lambda-spatial", type=float, default=8.0)
    ap.add_argument("--wta-smooth-iters", type=int, default=1)
    ap.add_argument("--min-wta-comp-size", type=int, default=20)
    ap.add_argument("--max-wta-prune-pass", type=int, default=2)
    ap.add_argument("--merge-thresh", type=float, default=0.08)

    ap.add_argument("--fc-min-distance-mm", type=float, default=30.0)
    ap.add_argument("--fc-ref-max", type=int, default=10000)
    ap.add_argument("--fc-ref-pca-dim", type=int, default=0)
    ap.add_argument("--ref-seed", type=int, default=1)
    ap.add_argument("--verbose", type=int, default=1)
    return ap


def main() -> int:
    args = build_arg_parser().parse_args()
    verbose = int(args.verbose) == 1

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cimg = nib.load(args.in_cifti)
    axis = cimg.header.get_axis(1)
    n_gray = axis.size

    limg = nib.load(args.wta_dlabel)
    wta = np.asanyarray(limg.dataobj).squeeze().astype(np.int32, copy=False)
    if wta.size != n_gray:
        raise ValueError(f"WTA label size ({wta.size}) != CIFTI grayordinates ({n_gray})")

    cortex_idx = cortical_grayordinates(axis)
    active_mask = wta > 0

    data_tg = np.asanyarray(cimg.dataobj).astype(np.float32, copy=False)
    x = row_l2_demean(data_tg.T)

    neighbors_cortex = load_neighbors_cortex(args.neighbors_mat)
    w = build_full_graph(axis, neighbors_cortex)
    adj = adjacency_list(w)

    dist = np.load(args.distance_npy, mmap_mode="r")
    dist_c = np.asarray(dist[np.ix_(cortex_idx, cortex_idx)], dtype=np.float64)

    ref_idx = select_reference_vertices(
        dist_c,
        min_mm=float(args.fc_min_distance_mm),
        ref_max=int(args.fc_ref_max),
        rngseed=int(args.ref_seed),
    )
    if verbose:
        print(f"[areal] selected {ref_idx.size} far-only FC references")

    y = row_l2_demean(x[cortex_idx[ref_idx], :])
    if int(args.fc_ref_pca_dim) > 0 and int(args.fc_ref_pca_dim) < y.shape[1]:
        # Keep optional path for parity; defaults keep this disabled.
        y = pca_reduce_rows(y.T, int(args.fc_ref_pca_dim)).T
        xf = pca_reduce_rows(x.T, int(args.fc_ref_pca_dim)).T
    else:
        y = y.T
        xf = x

    f = xf @ y
    f = f / (np.sqrt((f * f).sum(axis=1, keepdims=True)) + 1e-8)

    grad = compute_local_gradient_fcfp(w, f)

    if verbose:
        print(
            f"[wta] smoothing {args.wta_smooth_iters} iters + pruning CCs < {args.min_wta_comp_size} "
            f"(max {args.max_wta_prune_pass} passes)..."
        )
    wta_clean = clean_wta_majority(wta, adj, int(args.wta_smooth_iters))
    wta_clean = prune_small_wta_components(
        wta_clean,
        adj,
        min_comp=int(args.min_wta_comp_size),
        max_pass=int(args.max_wta_prune_pass),
        verbose=verbose,
    )

    patch_id, n_patches = wta_connected_components(wta_clean, adj, active_mask)
    parcel_labels = np.zeros((n_gray,), dtype=np.int32)
    parcel_net = []
    parcel_ctr = 0
    kchosen = np.zeros((n_patches,), dtype=np.float32)

    for p in range(1, n_patches + 1):
        patch_mask = patch_id == p
        if not patch_mask.any():
            continue
        vidx = np.where(patch_mask)[0]
        vals = wta_clean[vidx]
        vals = vals[vals > 0]
        if vals.size == 0:
            continue
        netp = mode_nonzero(vals, fallback=0)

        k0, kmin_p, kmax_p = suggest_k_bounds_per_patch(vidx, None, args)
        if int(args.k_fixed) > 0:
            k_use = max(kmin_p, min(kmax_p, int(args.k_fixed)))
        else:
            k_use = auto_k_in_range_fcfp(
                vidx=vidx,
                w=w,
                f=f,
                k0=k0,
                kmin_p=kmin_p,
                kmax_p=kmax_p,
                opts=args,
            )
        kchosen[p - 1] = float(k_use)

        ploc = split_kmeans_fcfp(vidx, w, f, int(k_use), args)

        ploc = remove_tiny_and_soft_merge(ploc, patch_mask, adj, grad, None, args)

        u = np.unique(ploc[patch_mask])
        u = u[u > 0]
        for uu in u:
            parcel_ctr += 1
            sel = patch_mask & (ploc == uu)
            parcel_labels[sel] = parcel_ctr
            parcel_net.append(int(netp))

        if verbose:
            print(f"[{p:03d}/{n_patches:03d}] verts={vidx.size} net={netp} k={int(k_use)}")

    parcel_net_arr = np.asarray(parcel_net, dtype=np.int32)
    if verbose:
        print(f"[summary] patches={n_patches} parcels={int(parcel_labels.max())}")

    pri = load_priors(args.priors_mat)
    if args.priors_mat and pri is None:
        print("[areal] WARNING: proceeding without network labels/colors; dlabel colors will use fallback gray")
    elif pri is not None and (not hasattr(pri, "NetworkLabels") or not hasattr(pri, "NetworkColors")):
        print("[areal] WARNING: priors missing NetworkLabels/NetworkColors; dlabel colors will use fallback gray")

    tmp = outdir / "Tmp_Areal.dtseries.nii"
    save_scalar_like(cimg, parcel_labels.astype(np.float32), tmp)
    out_dlabel = outdir / f"{args.outfile}.dlabel.nii"
    labfile = outdir / "LabelListFile_Areal.txt"
    write_parcel_label_list(labfile, parcel_net_arr, pri)

    subprocess.run(
        ["wb_command", "-cifti-label-import", str(tmp), str(labfile), str(out_dlabel), "-discard-others"],
        check=True,
    )
    subprocess.run(
        ["wb_command", "-cifti-label-to-border", str(out_dlabel), "-border", args.left_surf, str(outdir / f"{args.outfile}.L.border")],
        check=True,
    )
    subprocess.run(
        ["wb_command", "-cifti-label-to-border", str(out_dlabel), "-border", args.right_surf, str(outdir / f"{args.outfile}.R.border")],
        check=True,
    )

    if tmp.exists():
        tmp.unlink()
    if labfile.exists():
        labfile.unlink()

    print(f"[areal] wrote {out_dlabel}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
