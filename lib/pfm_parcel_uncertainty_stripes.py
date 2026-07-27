#!/usr/bin/env python3
"""Render parcel-level Ridge-Fusion ambiguity as whole-brain network stripes.

The input areal parcellation supplies each parcel's assigned network.  The
Ridge-Fusion probability maps are averaged within each parcel and used to find
the strongest alternative network. Qualifying cortical parcels are rendered
with alternating latitude bands on the spherical surface. Qualifying
subcortical parcels are rendered both as alternating CIFTI voxels and as an
upsampled label volume with stripes inside each native voxel. These follow the
visual conventions of the legacy ``make_striped_cifti.m`` and
``subcort_network_stripes_nifti.m`` helpers.

This is deliberately a visualization/diagnostic pass: it does not replace the
probability maps and it does not claim that model uncertainty proves biological
network overlap.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
import scipy.io as sio
from nibabel.cifti2.cifti2_axes import LabelAxis, ScalarAxis


@dataclass
class ParcelResult:
    parcel_id: int
    parcel_label: str
    n_grayordinates: int
    n_cortical_vertices: int
    n_subcortical_voxels: int
    structures: str
    assigned_network_id: int
    assigned_network: str
    assigned_probability: float
    alternative_network_id: int
    alternative_network: str
    alternative_probability: float
    alternative_to_assigned_ratio: float
    second_alternative_network_id: int
    second_alternative_network: str
    second_alternative_probability: float
    second_alternative_to_assigned_ratio: float
    normalized_entropy: float
    stripe_network_count: int
    stripe_network_ids: str
    stripe_networks: str
    stripe_network_probabilities: str
    mixed_cortical_by_rule: bool
    mixed_subcortical_by_rule: bool
    mixed_by_rule: bool
    rendered_as_stripes: bool
    assignment_source: str


def load_priors(priors_mat: Path, n_networks: int) -> tuple[list[str], np.ndarray]:
    pri = sio.loadmat(str(priors_mat), squeeze_me=True, struct_as_record=False)["Priors"]
    names = [str(x).strip() for x in np.asarray(pri.NetworkLabels).ravel()]
    colors = np.asarray(pri.NetworkColors, dtype=np.float64)
    if len(names) < n_networks or colors.shape[0] < n_networks:
        raise ValueError(
            f"Priors contain {len(names)} names/{colors.shape[0]} colors but "
            f"the probability CIFTI contains {n_networks} maps"
        )
    names = names[:n_networks]
    colors = colors[:n_networks, :3]
    if np.nanmax(colors) > 1.0:
        colors = colors / 255.0
    return names, np.clip(colors, 0.0, 1.0)


def surface_coordinates(surface_path: Path) -> np.ndarray:
    img = nib.load(str(surface_path))
    pointset_code = nib.nifti1.intent_codes.code["NIFTI_INTENT_POINTSET"]
    for darray in img.darrays:
        if int(darray.intent) == int(pointset_code):
            coords = np.asarray(darray.data, dtype=np.float64)
            if coords.ndim == 2 and coords.shape[1] == 3:
                return coords
    raise ValueError(f"No POINTSET coordinates found in surface: {surface_path}")


def grayordinate_geometry(
    brain_axis,
    left_sphere: Path,
    right_sphere: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return cortical mask, latitude, structure names, and volume IJK."""
    n_gray = int(brain_axis.size)
    is_cortex = np.zeros(n_gray, dtype=bool)
    latitude = np.full(n_gray, np.nan, dtype=np.float64)
    structure = np.full(n_gray, "", dtype=object)
    voxel_ijk = np.full((n_gray, 3), -1, dtype=np.int32)
    spheres = {
        "CIFTI_STRUCTURE_CORTEX_LEFT": surface_coordinates(left_sphere),
        "CIFTI_STRUCTURE_CORTEX_RIGHT": surface_coordinates(right_sphere),
    }

    for name, slc, model in brain_axis.iter_structures():
        stop = slc.stop if slc.stop is not None else n_gray
        gidx = np.arange(slc.start, stop, dtype=np.int64)
        structure[gidx] = name.replace("CIFTI_STRUCTURE_", "")
        if name not in spheres:
            voxels = np.asarray(model.voxel, dtype=np.int64)
            if voxels.shape != (gidx.size, 3):
                raise ValueError(f"Voxel count mismatch for {name}: {voxels.shape} != ({gidx.size}, 3)")
            voxel_ijk[gidx] = voxels.astype(np.int32)
            continue
        vertices = np.asarray(model.vertex, dtype=np.int64)
        if vertices.size != gidx.size:
            raise ValueError(f"Vertex count mismatch for {name}: {vertices.size} != {gidx.size}")
        coords = spheres[name][vertices]
        # MATLAB cart2sph elevation (theta), as used by make_striped_cifti.m.
        latitude[gidx] = np.arctan2(coords[:, 2], np.hypot(coords[:, 0], coords[:, 1]))
        is_cortex[gidx] = True
    return is_cortex, latitude, structure, voxel_ijk


def parcel_network_from_label(parcel_id: int, parcel_label: str, network_names: list[str]) -> tuple[int, str]:
    suffix = f"_{parcel_id}"
    base = parcel_label[: -len(suffix)] if parcel_label.endswith(suffix) else parcel_label
    try:
        return network_names.index(base), "parcel_label_table"
    except ValueError:
        return -1, "probability_argmax_fallback"


def normalized_entropy(probability: np.ndarray) -> float:
    p = np.asarray(probability, dtype=np.float64)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    return float(-(p * np.log(p)).sum() / np.log(probability.size))


def load_parcel_support_csv(path: Path) -> dict[int, list[int]]:
    """Read one-based network sets from CV or parcel-stripe summaries."""
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return {}
    if "selected_network_ids" in rows[0]:
        id_field = "selected_network_ids"
        count_field = "selected_network_count"
    elif "stripe_network_ids" in rows[0]:
        id_field = "stripe_network_ids"
        count_field = "stripe_network_count"
    else:
        raise ValueError("Parcel support CSV needs selected_network_ids or stripe_network_ids")
    support = {}
    for row in rows:
        count = int(row[count_field])
        ids = [int(x) for x in row[id_field].split(";") if x]
        if count >= 2 and len(ids) >= count:
            support[int(row["parcel_id"])] = ids[:count]
    return support


def stripe_assignments(
    latitude: np.ndarray,
    width_radians: float,
    primary_id: int,
    alternative_id: int,
) -> np.ndarray:
    """Backward-compatible two-network surface stripe helper."""
    return multi_stripe_assignments(latitude, width_radians, [primary_id, alternative_id])


def multi_stripe_assignments(
    latitude: np.ndarray,
    width_radians: float,
    network_ids: list[int],
) -> np.ndarray:
    """Assign stable global latitude stripes for two or more networks."""
    ids = np.asarray(network_ids, dtype=np.int32)
    if ids.size < 1:
        raise ValueError("At least one network ID is required")
    if latitude.size < ids.size:
        return np.full(latitude.shape, ids[0], dtype=np.int32)
    bands = np.floor((latitude + (np.pi / 2.0)) / width_radians).astype(np.int64)
    phase = bands % ids.size

    # Very narrow parcels can land entirely within one global band.  Preserve
    # all supported colors with ordered local bands in that case.
    if np.unique(phase).size < ids.size:
        order = np.argsort(latitude, kind="stable")
        phase = np.zeros(latitude.size, dtype=np.int64)
        local_phase = np.minimum((np.arange(latitude.size) * ids.size) // latitude.size, ids.size - 1)
        phase[order] = local_phase

    return ids[phase]


def stripe_plane_axes(view: str) -> tuple[int, int]:
    if view == "axial":
        return 0, 1
    if view == "coronal":
        return 0, 2
    if view == "sagittal":
        return 1, 2
    raise ValueError(f"Unknown subcortical stripe view: {view}")


def voxel_stripe_assignments(
    voxel_ijk: np.ndarray,
    width_voxels: int,
    view: str,
    primary_id: int,
    alternative_id: int,
) -> np.ndarray:
    """Backward-compatible two-network volume stripe helper."""
    return multi_voxel_stripe_assignments(
        voxel_ijk,
        width_voxels,
        view,
        [primary_id, alternative_id],
    )


def multi_voxel_stripe_assignments(
    voxel_ijk: np.ndarray,
    width_voxels: int,
    view: str,
    network_ids: list[int],
) -> np.ndarray:
    """Assign diagonal stripes to existing CIFTI voxels for multiple networks."""
    ids = np.asarray(network_ids, dtype=np.int32)
    if ids.size < 1:
        raise ValueError("At least one network ID is required")
    if voxel_ijk.shape[0] < ids.size:
        return np.full(voxel_ijk.shape[0], ids[0], dtype=np.int32)
    axis_a, axis_b = stripe_plane_axes(view)
    coordinate = voxel_ijk[:, axis_a] + voxel_ijk[:, axis_b]
    bands = np.floor_divide(coordinate, int(width_voxels))
    phase = bands % ids.size
    if np.unique(phase).size < ids.size:
        order = np.lexsort((voxel_ijk[:, 2], voxel_ijk[:, 1], voxel_ijk[:, 0], coordinate))
        phase = np.zeros(voxel_ijk.shape[0], dtype=np.int64)
        local_phase = np.minimum(
            (np.arange(voxel_ijk.shape[0]) * ids.size) // voxel_ijk.shape[0],
            ids.size - 1,
        )
        phase[order] = local_phase
    return ids[phase]


def write_workbench_label_list(
    out_path: Path,
    network_names: list[str],
    network_colors: np.ndarray,
) -> None:
    with out_path.open("w") as f:
        for network_id, (name, rgb) in enumerate(zip(network_names, network_colors), start=1):
            color = np.clip(np.rint(255.0 * rgb), 0, 255).astype(int)
            f.write(f"{name}\n{network_id} {color[0]} {color[1]} {color[2]} 255\n")


def highres_affine_and_factors(
    affine: np.ndarray,
    output_resolution_mm: float,
) -> tuple[np.ndarray, np.ndarray]:
    voxel_sizes = np.sqrt((np.asarray(affine[:3, :3], dtype=np.float64) ** 2).sum(axis=0))
    factors = np.rint(voxel_sizes / output_resolution_mm).astype(np.int32)
    if np.any(factors < 1):
        raise ValueError("Subcortical output resolution cannot exceed a native voxel dimension")
    achieved = voxel_sizes / factors
    if not np.allclose(achieved, output_resolution_mm, atol=1e-4):
        raise ValueError(
            f"Requested {output_resolution_mm:g} mm does not divide native voxel sizes {voxel_sizes.tolist()}"
        )
    highres_affine = np.asarray(affine, dtype=np.float64).copy()
    highres_affine[:3, :3] = affine[:3, :3] @ np.diag(1.0 / factors)
    highres_affine[:3, 3] = affine[:3, 3] - highres_affine[:3, :3] @ ((factors - 1.0) / 2.0)
    return highres_affine, factors


def build_subcortical_highres_volume(
    volume_shape: tuple[int, int, int],
    affine: np.ndarray,
    voxel_ijk: np.ndarray,
    stripe_network_maps: np.ndarray,
    stripe_network_count_map: np.ndarray,
    mixed_subcortical_mask: np.ndarray,
    output_resolution_mm: float,
    stripe_width_mm: float,
    view: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Upsample labels and place both colors inside each selected native voxel."""
    highres_affine, factors = highres_affine_and_factors(affine, output_resolution_mm)
    primary_map = stripe_network_maps[0]
    if int(primary_map.max()) > np.iinfo(np.uint8).max:
        raise ValueError("Subcortical label volume currently supports at most 255 network IDs")

    native = np.zeros(volume_shape, dtype=np.uint8)
    has_label = (voxel_ijk[:, 0] >= 0) & (primary_map > 0)
    voxels = voxel_ijk[has_label]
    native[voxels[:, 0], voxels[:, 1], voxels[:, 2]] = primary_map[has_label].astype(np.uint8)
    highres = native
    for axis, factor in enumerate(factors):
        highres = np.repeat(highres, int(factor), axis=axis)

    axis_a, axis_b = stripe_plane_axes(view)
    width_subvoxels = max(1, int(round(stripe_width_mm / output_resolution_mm)))
    local_grid = np.indices(tuple(int(x) for x in factors), dtype=np.int32)
    for gray_index in np.where(mixed_subcortical_mask)[0]:
        network_count = int(stripe_network_count_map[gray_index])
        network_ids = stripe_network_maps[:network_count, gray_index].astype(np.uint8)
        voxel = voxel_ijk[gray_index]
        base = voxel * factors
        global_a = local_grid[axis_a] + int(base[axis_a])
        global_b = local_grid[axis_b] + int(base[axis_b])
        phase = ((global_a + global_b) // width_subvoxels) % network_count
        if np.unique(phase).size < network_count:
            flat_order = np.argsort((global_a + global_b).ravel(), kind="stable")
            flat_phase = np.zeros(phase.size, dtype=np.int64)
            local_phase = np.minimum(
                (np.arange(phase.size) * network_count) // phase.size,
                network_count - 1,
            )
            flat_phase[flat_order] = local_phase
            phase = flat_phase.reshape(phase.shape)
        block = network_ids[phase]
        slices = tuple(slice(int(base[d]), int(base[d] + factors[d])) for d in range(3))
        highres[slices] = block
    return highres, highres_affine


def save_subcortical_label_volume(
    values: np.ndarray,
    affine: np.ndarray,
    network_names: list[str],
    network_colors: np.ndarray,
    out_path: Path,
    label_list_path: Path,
    wb_command: str,
) -> None:
    tmp_path = out_path.with_name(f"Tmp_{out_path.name}")
    image = nib.Nifti1Image(values, affine)
    image.header.set_xyzt_units("mm")
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, str(tmp_path))
    write_workbench_label_list(label_list_path, network_names, network_colors)
    try:
        subprocess.run(
            [
                wb_command,
                "-volume-label-import",
                str(tmp_path),
                str(label_list_path),
                str(out_path),
                "-discard-others",
                "-unlabeled-value",
                "0",
            ],
            check=True,
        )
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def save_dlabel(
    reference_img: nib.Cifti2Image,
    values: np.ndarray,
    network_names: list[str],
    network_colors: np.ndarray,
    out_path: Path,
    map_name: str = "Parcel multi-network stripes",
) -> None:
    table: dict[int, tuple[str, tuple[float, float, float, float]]] = {
        0: ("Unlabeled", (0.0, 0.0, 0.0, 0.0))
    }
    for network_id, (name, rgb) in enumerate(zip(network_names, network_colors), start=1):
        table[network_id] = (name, (float(rgb[0]), float(rgb[1]), float(rgb[2]), 1.0))
    brain_axis = reference_img.header.get_axis(1)
    header = nib.Cifti2Header.from_axes((LabelAxis([map_name], [table]), brain_axis))
    image = nib.Cifti2Image(
        values.reshape(1, -1).astype(np.float32),
        header,
        nifti_header=reference_img.nifti_header,
    )
    nib.save(image, str(out_path))


def save_dscalar(
    reference_img: nib.Cifti2Image,
    maps: np.ndarray,
    map_names: list[str],
    out_path: Path,
) -> None:
    brain_axis = reference_img.header.get_axis(1)
    header = nib.Cifti2Header.from_axes((ScalarAxis(map_names), brain_axis))
    image = nib.Cifti2Image(maps.astype(np.float32), header, nifti_header=reference_img.nifti_header)
    nib.save(image, str(out_path))


def write_csv(rows: list[ParcelResult], out_path: Path) -> None:
    fieldnames = list(asdict(rows[0]).keys()) if rows else [f.name for f in ParcelResult.__dataclass_fields__.values()]
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Create parcel-level multi-network stripe diagnostics")
    ap.add_argument("--parcels-dlabel", required=True, type=Path)
    ap.add_argument("--prob-maps-cifti", required=True, type=Path)
    ap.add_argument(
        "--parcel-support-csv",
        type=Path,
        help="Optional CV-selected parcel/network sets; overrides the probability threshold rule",
    )
    ap.add_argument("--priors-mat", required=True, type=Path)
    ap.add_argument("--left-sphere", required=True, type=Path)
    ap.add_argument("--right-sphere", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--outfile-prefix", default="RidgeFusion_ArealMultiNetwork")
    ap.add_argument("--secondary-prob-min", type=float, default=0.20)
    ap.add_argument("--secondary-to-assigned-ratio-min", type=float, default=2.0 / 3.0)
    ap.add_argument("--max-networks-per-parcel", type=int, choices=[2, 3], default=3)
    ap.add_argument("--min-parcel-vertices", type=int, default=20)
    ap.add_argument("--stripe-width-radians", type=float, default=1.0 / 40.0)
    ap.add_argument("--min-subcortical-parcel-voxels", type=int, default=5)
    ap.add_argument("--subcortical-native-stripe-width-voxels", type=int, default=1)
    ap.add_argument("--subcortical-output-resolution-mm", type=float, default=0.5)
    ap.add_argument("--subcortical-stripe-width-mm", type=float, default=0.5)
    ap.add_argument(
        "--subcortical-stripe-view",
        choices=["axial", "coronal", "sagittal"],
        default="axial",
    )
    ap.add_argument("--wb-command", default="wb_command")
    return ap


def main() -> int:
    args = build_arg_parser().parse_args()
    if not 0.0 <= args.secondary_prob_min <= 1.0:
        raise ValueError("--secondary-prob-min must be between 0 and 1")
    if args.secondary_to_assigned_ratio_min < 0.0:
        raise ValueError("--secondary-to-assigned-ratio-min must be nonnegative")
    if args.min_parcel_vertices < 2:
        raise ValueError("--min-parcel-vertices must be at least 2")
    if args.stripe_width_radians <= 0.0:
        raise ValueError("--stripe-width-radians must be positive")
    if args.min_subcortical_parcel_voxels < 2:
        raise ValueError("--min-subcortical-parcel-voxels must be at least 2")
    if args.subcortical_native_stripe_width_voxels < 1:
        raise ValueError("--subcortical-native-stripe-width-voxels must be at least 1")
    if args.subcortical_output_resolution_mm <= 0.0:
        raise ValueError("--subcortical-output-resolution-mm must be positive")
    if args.subcortical_stripe_width_mm <= 0.0:
        raise ValueError("--subcortical-stripe-width-mm must be positive")

    args.outdir.mkdir(parents=True, exist_ok=True)
    parcel_img = nib.load(str(args.parcels_dlabel))
    parcel_axis = parcel_img.header.get_axis(0)
    if not isinstance(parcel_axis, LabelAxis) or parcel_img.shape[0] != 1:
        raise ValueError("--parcels-dlabel must be a single-map CIFTI dlabel")
    parcel_ids = np.asanyarray(parcel_img.dataobj).reshape(-1).astype(np.int32)
    brain_axis = parcel_img.header.get_axis(1)

    prob_img = nib.load(str(args.prob_maps_cifti))
    probabilities = np.asanyarray(prob_img.dataobj).astype(np.float64, copy=False).T
    if probabilities.shape[0] != parcel_ids.size:
        raise ValueError(
            f"Probability grayordinates ({probabilities.shape[0]}) do not match parcels ({parcel_ids.size})"
        )
    if np.nanmin(probabilities) < -1e-6:
        raise ValueError("Probability CIFTI contains negative values")
    row_sum = probabilities.sum(axis=1, keepdims=True)
    valid_probability = np.isfinite(probabilities).all(axis=1) & (row_sum[:, 0] > 0)
    probabilities = np.divide(
        probabilities,
        row_sum,
        out=np.zeros_like(probabilities),
        where=row_sum > 0,
    )

    network_names, network_colors = load_priors(args.priors_mat, probabilities.shape[1])
    is_cortex, latitude, structures, voxel_ijk = grayordinate_geometry(
        brain_axis,
        args.left_sphere,
        args.right_sphere,
    )
    is_subcortex = voxel_ijk[:, 0] >= 0
    label_table = parcel_axis.label[0]
    external_support = (
        load_parcel_support_csv(args.parcel_support_csv)
        if args.parcel_support_csv is not None
        else None
    )

    stripe_network_maps = np.zeros(
        (args.max_networks_per_parcel, parcel_ids.size),
        dtype=np.int32,
    )
    stripe_network_count_map = np.ones(parcel_ids.size, dtype=np.int32)
    stripe_map = np.zeros(parcel_ids.size, dtype=np.int32)
    mixed_subcortical_mask = np.zeros(parcel_ids.size, dtype=bool)
    scalar_maps = np.zeros((9, parcel_ids.size), dtype=np.float32)
    rows: list[ParcelResult] = []

    for parcel_id in np.unique(parcel_ids[parcel_ids > 0]):
        parcel_id = int(parcel_id)
        in_parcel = parcel_ids == parcel_id
        cortical = in_parcel & is_cortex
        subcortical = in_parcel & is_subcortex
        usable = in_parcel & valid_probability
        parcel_label = label_table.get(parcel_id, (f"Parcel_{parcel_id}", (0, 0, 0, 0)))[0]
        assigned_idx, assignment_source = parcel_network_from_label(parcel_id, parcel_label, network_names)

        if usable.any():
            parcel_probability = probabilities[usable].mean(axis=0)
        else:
            parcel_probability = np.zeros(probabilities.shape[1], dtype=np.float64)
        if assigned_idx < 0:
            assigned_idx = int(np.argmax(parcel_probability))

        probability_alternative_order = [
            int(i) for i in np.argsort(parcel_probability)[::-1] if int(i) != assigned_idx
        ]
        if external_support is not None:
            selected_ids = external_support.get(parcel_id, [])
            selected_indices = [network_id - 1 for network_id in selected_ids]
            if any(index < 0 or index >= probabilities.shape[1] for index in selected_indices):
                raise ValueError(f"Parcel {parcel_id} support contains an invalid network ID")
            if selected_indices and assigned_idx not in selected_indices:
                raise ValueError(
                    f"Parcel {parcel_id} external support omits assigned network {assigned_idx + 1}"
                )
            selected_indices = (
                [assigned_idx] + [index for index in selected_indices if index != assigned_idx]
            )[: args.max_networks_per_parcel]
            stripe_network_indices = selected_indices if len(selected_indices) >= 2 else [assigned_idx]
            support_rule = len(stripe_network_indices) >= 2
            assignment_source = f"{assignment_source}+external_parcel_support"
        else:
            qualifying_alternatives = [
                network_idx
                for network_idx in probability_alternative_order
                if parcel_probability[network_idx] >= args.secondary_prob_min
                and parcel_probability[network_idx]
                / max(float(parcel_probability[assigned_idx]), np.finfo(float).eps)
                >= args.secondary_to_assigned_ratio_min
            ][: args.max_networks_per_parcel - 1]
            stripe_network_indices = [assigned_idx] + qualifying_alternatives
            support_rule = len(qualifying_alternatives) > 0

        selected_alternatives = stripe_network_indices[1:]
        alternative_order = selected_alternatives + [
            index for index in probability_alternative_order if index not in selected_alternatives
        ]
        alternative_idx = alternative_order[0]
        second_alternative_idx = alternative_order[1]
        assigned_probability = float(parcel_probability[assigned_idx])
        alternative_probability = float(parcel_probability[alternative_idx])
        ratio = alternative_probability / max(assigned_probability, np.finfo(float).eps)
        second_alternative_probability = float(parcel_probability[second_alternative_idx])
        second_ratio = second_alternative_probability / max(
            assigned_probability,
            np.finfo(float).eps,
        )
        entropy = normalized_entropy(parcel_probability)
        stripe_network_ids = [network_idx + 1 for network_idx in stripe_network_indices]
        mixed_cortical = bool(
            cortical.sum() >= args.min_parcel_vertices
            and support_rule
        )
        mixed_subcortical = bool(
            subcortical.sum() >= args.min_subcortical_parcel_voxels
            and support_rule
        )
        mixed = mixed_cortical or mixed_subcortical

        for rank, network_id in enumerate(stripe_network_ids):
            stripe_network_maps[rank, in_parcel] = network_id
        stripe_network_count_map[in_parcel] = len(stripe_network_ids)
        stripe_map[in_parcel] = assigned_idx + 1
        rendered = False
        if mixed_cortical:
            cidx = np.where(cortical)[0]
            stripe_map[cidx] = multi_stripe_assignments(
                latitude[cidx],
                args.stripe_width_radians,
                stripe_network_ids,
            )
            rendered = bool(np.unique(stripe_map[cidx]).size == len(stripe_network_ids))
        if mixed_subcortical:
            sidx = np.where(subcortical)[0]
            stripe_map[sidx] = multi_voxel_stripe_assignments(
                voxel_ijk[sidx],
                args.subcortical_native_stripe_width_voxels,
                args.subcortical_stripe_view,
                stripe_network_ids,
            )
            mixed_subcortical_mask[sidx] = True
            rendered = rendered or bool(
                np.unique(stripe_map[sidx]).size == len(stripe_network_ids)
            )

        scalar_maps[0, in_parcel] = assigned_probability
        scalar_maps[1, in_parcel] = alternative_probability
        scalar_maps[2, in_parcel] = ratio
        scalar_maps[3, in_parcel] = second_alternative_probability
        scalar_maps[4, in_parcel] = second_ratio
        scalar_maps[5, in_parcel] = entropy
        scalar_maps[6, in_parcel] = float(len(stripe_network_ids) if mixed else 1)
        scalar_maps[7, in_parcel] = float(mixed)
        scalar_maps[8, in_parcel] = float(parcel_id)

        structure_names = sorted({str(x) for x in structures[in_parcel] if str(x)})
        rows.append(
            ParcelResult(
                parcel_id=parcel_id,
                parcel_label=parcel_label,
                n_grayordinates=int(in_parcel.sum()),
                n_cortical_vertices=int(cortical.sum()),
                n_subcortical_voxels=int(subcortical.sum()),
                structures=";".join(structure_names),
                assigned_network_id=assigned_idx + 1,
                assigned_network=network_names[assigned_idx],
                assigned_probability=assigned_probability,
                alternative_network_id=alternative_idx + 1,
                alternative_network=network_names[alternative_idx],
                alternative_probability=alternative_probability,
                alternative_to_assigned_ratio=ratio,
                second_alternative_network_id=second_alternative_idx + 1,
                second_alternative_network=network_names[second_alternative_idx],
                second_alternative_probability=second_alternative_probability,
                second_alternative_to_assigned_ratio=second_ratio,
                normalized_entropy=entropy,
                stripe_network_count=len(stripe_network_ids) if mixed else 1,
                stripe_network_ids=";".join(str(x) for x in stripe_network_ids if mixed)
                if mixed
                else str(assigned_idx + 1),
                stripe_networks=";".join(network_names[x] for x in stripe_network_indices)
                if mixed
                else network_names[assigned_idx],
                stripe_network_probabilities=";".join(
                    f"{parcel_probability[x]:.9g}" for x in stripe_network_indices
                )
                if mixed
                else f"{assigned_probability:.9g}",
                mixed_cortical_by_rule=mixed_cortical,
                mixed_subcortical_by_rule=mixed_subcortical,
                mixed_by_rule=mixed,
                rendered_as_stripes=rendered,
                assignment_source=assignment_source,
            )
        )

    prefix = args.outfile_prefix
    dlabel_path = args.outdir / f"{prefix}_NetworkStripes.dlabel.nii"
    dscalar_path = args.outdir / f"{prefix}_ParcelAmbiguity.dscalar.nii"
    csv_path = args.outdir / f"{prefix}_ParcelSummary.csv"
    json_path = args.outdir / f"{prefix}_Summary.json"
    resolution_tag = f"{args.subcortical_output_resolution_mm:g}".replace(".", "p")
    subcortical_volume_path = args.outdir / f"{prefix}_SubcorticalNetworkStripes_{resolution_tag}mm.nii.gz"
    subcortical_label_list_path = args.outdir / f"{prefix}_SubcorticalNetworkStripes_LabelTable.txt"
    save_dlabel(parcel_img, stripe_map, network_names, network_colors, dlabel_path)
    save_dscalar(
        parcel_img,
        scalar_maps,
        [
            "Assigned network probability",
            "Strongest alternative probability",
            "Alternative / assigned probability",
            "Second alternative probability",
            "Second alternative / assigned probability",
            "Normalized network entropy",
            "Stripe network count",
            "Mixed parcel flag",
            "Original parcel ID",
        ],
        dscalar_path,
    )
    write_csv(rows, csv_path)
    subcortical_volume, subcortical_affine = build_subcortical_highres_volume(
        tuple(int(x) for x in brain_axis.volume_shape),
        np.asarray(brain_axis.affine, dtype=np.float64),
        voxel_ijk,
        stripe_network_maps,
        stripe_network_count_map,
        mixed_subcortical_mask,
        args.subcortical_output_resolution_mm,
        args.subcortical_stripe_width_mm,
        args.subcortical_stripe_view,
    )
    save_subcortical_label_volume(
        subcortical_volume,
        subcortical_affine,
        network_names,
        network_colors,
        subcortical_volume_path,
        subcortical_label_list_path,
        args.wb_command,
    )

    mixed_rows = [row for row in rows if row.mixed_by_rule]
    mixed_cortical_rows = [row for row in rows if row.mixed_cortical_by_rule]
    mixed_subcortical_rows = [row for row in rows if row.mixed_subcortical_by_rule]
    pair_counts = Counter(
        " + ".join(sorted((row.assigned_network, row.alternative_network))) for row in mixed_rows
    )
    combination_counts = Counter(
        " + ".join(sorted(row.stripe_networks.split(";"))) for row in mixed_rows
    )
    summary = {
        "inputs": {
            "parcels_dlabel": str(args.parcels_dlabel.resolve()),
            "prob_maps_cifti": str(args.prob_maps_cifti.resolve()),
            "priors_mat": str(args.priors_mat.resolve()),
            "left_sphere": str(args.left_sphere.resolve()),
            "right_sphere": str(args.right_sphere.resolve()),
            "parcel_support_csv": (
                str(args.parcel_support_csv.resolve()) if args.parcel_support_csv is not None else None
            ),
        },
        "criterion": {
            "selection_source": (
                "external_parcel_support_csv"
                if args.parcel_support_csv is not None
                else "parcel_mean_probability_threshold"
            ),
            "secondary_probability_min": args.secondary_prob_min,
            "secondary_to_assigned_ratio_min": args.secondary_to_assigned_ratio_min,
            "max_networks_per_parcel": args.max_networks_per_parcel,
            "min_cortical_parcel_vertices": args.min_parcel_vertices,
            "stripe_width_radians": args.stripe_width_radians,
            "min_subcortical_parcel_voxels": args.min_subcortical_parcel_voxels,
            "subcortical_native_stripe_width_voxels": args.subcortical_native_stripe_width_voxels,
            "subcortical_output_resolution_mm": args.subcortical_output_resolution_mm,
            "subcortical_stripe_width_mm": args.subcortical_stripe_width_mm,
            "subcortical_stripe_view": args.subcortical_stripe_view,
        },
        "counts": {
            "parcels_total": len(rows),
            "parcels_mixed": len(mixed_rows),
            "parcels_mixed_cortical": len(mixed_cortical_rows),
            "parcels_mixed_subcortical": len(mixed_subcortical_rows),
            "parcels_striped_with_two_networks": sum(
                row.stripe_network_count == 2 for row in mixed_rows
            ),
            "parcels_striped_with_three_networks": sum(
                row.stripe_network_count == 3 for row in mixed_rows
            ),
            "parcels_rendered_as_stripes": sum(row.rendered_as_stripes for row in rows),
            "cortical_vertices_in_mixed_parcels": sum(
                row.n_cortical_vertices for row in mixed_cortical_rows
            ),
            "subcortical_voxels_in_mixed_parcels": sum(
                row.n_subcortical_voxels for row in mixed_subcortical_rows
            ),
            "grayordinates_with_zero_probability": int((~valid_probability).sum()),
            "parcel_assignment_fallbacks": sum(
                row.assignment_source == "probability_argmax_fallback" for row in rows
            ),
        },
        "mixed_network_pair_counts": dict(sorted(pair_counts.items(), key=lambda item: (-item[1], item[0]))),
        "mixed_network_combination_counts": dict(
            sorted(combination_counts.items(), key=lambda item: (-item[1], item[0]))
        ),
        "outputs": {
            "network_stripes_dlabel": str(dlabel_path.resolve()),
            "parcel_ambiguity_dscalar": str(dscalar_path.resolve()),
            "parcel_summary_csv": str(csv_path.resolve()),
            "subcortical_network_stripes_volume": str(subcortical_volume_path.resolve()),
            "subcortical_label_table": str(subcortical_label_list_path.resolve()),
        },
        "interpretation": (
            "Stripes denote parcel-level multi-network support under the external held-out/null rule; "
            "they do not by themselves establish biological network overlap."
            if args.parcel_support_csv is not None
            else
            "Stripes denote parcel-level ambiguity under the stated Ridge-Fusion probability rule; "
            "they do not by themselves establish biological network overlap."
        ),
    }
    with json_path.open("w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    print(
        f"[pfm-stripes] wrote {dlabel_path.name}; "
        f"striped {len(mixed_rows)}/{len(rows)} parcels "
        f"({summary['counts']['cortical_vertices_in_mixed_parcels']} cortical vertices, "
        f"{summary['counts']['subcortical_voxels_in_mixed_parcels']} subcortical voxels)"
    )
    print(f"[pfm-stripes] subcortical volume: {subcortical_volume_path.name}")
    print(f"[pfm-stripes] diagnostics: {dscalar_path.name}, {csv_path.name}, {json_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
