#!/usr/bin/env python3
"""Translate PFM multi-network support into clustered vertex/voxel stripes."""

from __future__ import annotations

import argparse
import csv
import json
from collections import deque
from pathlib import Path

import nibabel as nib
import numpy as np

from pfm_areal_parcellation import adjacency_list, build_full_graph, load_neighbors_cortex
from pfm_parcel_uncertainty_stripes import (
    build_subcortical_highres_volume,
    grayordinate_geometry,
    load_priors,
    multi_stripe_assignments,
    multi_voxel_stripe_assignments,
    save_dlabel,
    save_dscalar,
    save_subcortical_label_volume,
)


def load_parcel_support(path: Path) -> dict[int, list[int]]:
    support = {}
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return support
    if "selected_network_ids" in rows[0]:
        id_field = "selected_network_ids"
        count_field = "selected_network_count"
    elif "stripe_network_ids" in rows[0]:
        id_field = "stripe_network_ids"
        count_field = "stripe_network_count"
    else:
        raise ValueError("Parcel support CSV needs selected_network_ids or stripe_network_ids")
    for row in rows:
        ids = [int(x) for x in row[id_field].split(";") if x]
        count = int(row[count_field])
        if count >= 2 and len(ids) >= 2:
            support[int(row["parcel_id"])] = ids[:count]
    return support


def network_code(network_ids: list[int]) -> int:
    code = 0
    for network_id in network_ids:
        code |= 1 << (int(network_id) - 1)
    return code


def code_networks(code: int, n_networks: int) -> list[int]:
    return [network_id for network_id in range(1, n_networks + 1) if code & (1 << (network_id - 1))]


def connected_components_for_code(
    indices: np.ndarray,
    code: int,
    support_code: np.ndarray,
    adjacency: list[np.ndarray],
) -> list[np.ndarray]:
    allowed = np.zeros(support_code.size, dtype=bool)
    allowed[indices] = True
    seen = np.zeros(support_code.size, dtype=bool)
    components = []
    for seed in indices:
        if seen[seed]:
            continue
        queue = deque([int(seed)])
        seen[seed] = True
        members = []
        while queue:
            vertex = queue.popleft()
            members.append(vertex)
            for neighbor in adjacency[vertex]:
                neighbor = int(neighbor)
                if not seen[neighbor] and allowed[neighbor] and int(support_code[neighbor]) == code:
                    seen[neighbor] = True
                    queue.append(neighbor)
        components.append(np.asarray(members, dtype=np.int64))
    return components


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Cluster direct or parcel-gated vertexwise PFM ambiguity")
    ap.add_argument("--mode", choices=["direct", "parcel_gated"], required=True)
    ap.add_argument("--wta-dlabel", required=True, type=Path)
    ap.add_argument("--prob-maps-cifti", required=True, type=Path)
    ap.add_argument("--parcels-dlabel", type=Path)
    ap.add_argument("--parcel-support-csv", type=Path)
    ap.add_argument("--neighbors-mat", required=True, type=Path)
    ap.add_argument("--priors-mat", required=True, type=Path)
    ap.add_argument("--left-sphere", required=True, type=Path)
    ap.add_argument("--right-sphere", required=True, type=Path)
    ap.add_argument("--secondary-prob-min", type=float, default=0.20)
    ap.add_argument("--secondary-to-top-ratio-min", type=float, default=2.0 / 3.0)
    ap.add_argument("--max-networks-per-zone", type=int, choices=[2, 3], default=3)
    ap.add_argument("--min-cortical-cluster-vertices", type=int, default=20)
    ap.add_argument("--min-subcortical-cluster-voxels", type=int, default=5)
    ap.add_argument("--stripe-width-radians", type=float, default=1.0 / 40.0)
    ap.add_argument("--subcortical-native-stripe-width-voxels", type=int, default=1)
    ap.add_argument("--subcortical-output-resolution-mm", type=float, default=0.5)
    ap.add_argument("--subcortical-stripe-width-mm", type=float, default=0.5)
    ap.add_argument("--subcortical-stripe-view", choices=["axial", "coronal", "sagittal"], default="axial")
    ap.add_argument("--wb-command", default="wb_command")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--outfile-prefix", required=True)
    return ap


def main() -> int:
    args = build_arg_parser().parse_args()
    if args.mode == "parcel_gated" and (args.parcels_dlabel is None or args.parcel_support_csv is None):
        raise ValueError("parcel_gated mode requires --parcels-dlabel and --parcel-support-csv")
    if not 0.0 <= args.secondary_prob_min <= 1.0:
        raise ValueError("--secondary-prob-min must be in [0, 1]")

    args.outdir.mkdir(parents=True, exist_ok=True)
    wta_img = nib.load(str(args.wta_dlabel))
    wta = np.asanyarray(wta_img.dataobj).reshape(-1).astype(np.int32)
    prob_img = nib.load(str(args.prob_maps_cifti))
    probability = np.asanyarray(prob_img.dataobj).astype(np.float64).T
    probability = np.divide(
        probability,
        probability.sum(axis=1, keepdims=True),
        out=np.zeros_like(probability),
        where=probability.sum(axis=1, keepdims=True) > 0,
    )
    if probability.shape[0] != wta.size:
        raise ValueError("WTA and probability CIFTIs have different grayordinate counts")
    network_names, network_colors = load_priors(args.priors_mat, probability.shape[1])
    brain_axis = wta_img.header.get_axis(1)
    is_cortex, latitude, structures, voxel_ijk = grayordinate_geometry(
        brain_axis,
        args.left_sphere,
        args.right_sphere,
    )
    is_subcortex = voxel_ijk[:, 0] >= 0

    parcel_ids = np.zeros(wta.size, dtype=np.int32)
    parcel_support: dict[int, list[int]] = {}
    if args.mode == "parcel_gated":
        parcel_img = nib.load(str(args.parcels_dlabel))
        parcel_ids = np.asanyarray(parcel_img.dataobj).reshape(-1).astype(np.int32)
        if parcel_ids.size != wta.size:
            raise ValueError("Parcel and WTA CIFTIs have different grayordinate counts")
        parcel_support = load_parcel_support(args.parcel_support_csv)

    support_code = np.zeros(wta.size, dtype=np.int64)
    local_count = np.ones(wta.size, dtype=np.int32)
    top_probability = np.max(probability, axis=1)
    sorted_indices = np.argsort(probability, axis=1)[:, ::-1]
    second_probability = probability[np.arange(wta.size), sorted_indices[:, 1]]
    third_probability = probability[np.arange(wta.size), sorted_indices[:, 2]]

    for gray_index in np.where(wta > 0)[0]:
        if args.mode == "direct":
            candidates = [int(x) + 1 for x in sorted_indices[gray_index, : args.max_networks_per_zone]]
        else:
            candidates = parcel_support.get(int(parcel_ids[gray_index]), [])
            candidates = sorted(
                candidates,
                key=lambda network_id: probability[gray_index, network_id - 1],
                reverse=True,
            )[: args.max_networks_per_zone]
        if len(candidates) < 2:
            continue
        baseline = float(probability[gray_index, candidates[0] - 1])
        supported = [
            network_id
            for network_id in candidates
            if probability[gray_index, network_id - 1] >= args.secondary_prob_min
            and probability[gray_index, network_id - 1] / max(baseline, np.finfo(float).eps)
            >= args.secondary_to_top_ratio_min
        ]
        if len(supported) >= 2:
            supported = sorted(supported[: args.max_networks_per_zone])
            support_code[gray_index] = network_code(supported)
            local_count[gray_index] = len(supported)

    neighbors = load_neighbors_cortex(str(args.neighbors_mat))
    graph = build_full_graph(brain_axis, neighbors)
    adjacency = adjacency_list(graph)
    striped_map = wta.copy()
    retained_mask = np.zeros(wta.size, dtype=bool)
    retained_count = np.ones(wta.size, dtype=np.int32)
    stripe_network_maps = np.zeros((args.max_networks_per_zone, wta.size), dtype=np.int32)
    stripe_network_maps[0] = wta
    cluster_rows = []
    cluster_id = 0

    for code in np.unique(support_code[support_code > 0]):
        code = int(code)
        indices = np.where(support_code == code)[0]
        network_ids = code_networks(code, probability.shape[1])
        if args.mode == "parcel_gated":
            index_groups = [indices[parcel_ids[indices] == parcel_id] for parcel_id in np.unique(parcel_ids[indices])]
        else:
            index_groups = [indices]
        components = []
        for index_group in index_groups:
            components.extend(
                connected_components_for_code(index_group, code, support_code, adjacency)
            )
        for component in components:
            cluster_id += 1
            cortical = bool(np.all(is_cortex[component]))
            subcortical = bool(np.all(is_subcortex[component]))
            minimum = (
                args.min_cortical_cluster_vertices
                if cortical
                else args.min_subcortical_cluster_voxels
            )
            retained = bool((cortical or subcortical) and component.size >= minimum)
            structure_names = sorted({str(x) for x in structures[component] if str(x)})
            parcel_values = sorted({int(x) for x in parcel_ids[component] if int(x) > 0})
            cluster_rows.append(
                {
                    "cluster_id": cluster_id,
                    "mode": args.mode,
                    "network_count": len(network_ids),
                    "network_ids": ";".join(str(x) for x in network_ids),
                    "networks": ";".join(network_names[x - 1] for x in network_ids),
                    "n_grayordinates": int(component.size),
                    "domain": "cortex" if cortical else "subcortex" if subcortical else "mixed_or_invalid",
                    "structures": ";".join(structure_names),
                    "parcel_ids": ";".join(str(x) for x in parcel_values),
                    "retained": retained,
                }
            )
            if not retained:
                continue
            retained_mask[component] = True
            retained_count[component] = len(network_ids)
            for rank, network_id in enumerate(network_ids):
                stripe_network_maps[rank, component] = network_id
            if cortical:
                striped_map[component] = multi_stripe_assignments(
                    latitude[component],
                    args.stripe_width_radians,
                    network_ids,
                )
            else:
                striped_map[component] = multi_voxel_stripe_assignments(
                    voxel_ijk[component],
                    args.subcortical_native_stripe_width_voxels,
                    args.subcortical_stripe_view,
                    network_ids,
                )

    prefix = args.outfile_prefix
    dlabel_path = args.outdir / f"{prefix}_NetworkStripes.dlabel.nii"
    dscalar_path = args.outdir / f"{prefix}_Diagnostics.dscalar.nii"
    csv_path = args.outdir / f"{prefix}_ClusterSummary.csv"
    json_path = args.outdir / f"{prefix}_Summary.json"
    save_dlabel(
        wta_img,
        striped_map,
        network_names,
        network_colors,
        dlabel_path,
        map_name=f"{args.mode.replace('_', ' ')} clustered multi-network stripes",
    )
    save_dscalar(
        wta_img,
        np.vstack(
            (
                top_probability,
                second_probability,
                second_probability / np.maximum(top_probability, np.finfo(float).eps),
                third_probability,
                third_probability / np.maximum(top_probability, np.finfo(float).eps),
                local_count,
                retained_mask.astype(np.float32),
                retained_count,
            )
        ),
        [
            "Top local network probability",
            "Second local network probability",
            "Second / top local probability",
            "Third local network probability",
            "Third / top local probability",
            "Raw locally supported network count",
            "Retained multi-network cluster mask",
            "Retained stripe network count",
        ],
        dscalar_path,
    )
    with csv_path.open("w", newline="") as f:
        fieldnames = list(cluster_rows[0]) if cluster_rows else [
            "cluster_id", "mode", "network_count", "network_ids", "networks",
            "n_grayordinates", "domain", "structures", "parcel_ids", "retained",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(cluster_rows)

    resolution_tag = f"{args.subcortical_output_resolution_mm:g}".replace(".", "p")
    volume_path = args.outdir / f"{prefix}_SubcorticalNetworkStripes_{resolution_tag}mm.nii.gz"
    label_list_path = args.outdir / f"{prefix}_SubcorticalNetworkStripes_LabelTable.txt"
    volume, affine = build_subcortical_highres_volume(
        tuple(int(x) for x in brain_axis.volume_shape),
        np.asarray(brain_axis.affine),
        voxel_ijk,
        stripe_network_maps,
        retained_count,
        retained_mask & is_subcortex,
        args.subcortical_output_resolution_mm,
        args.subcortical_stripe_width_mm,
        args.subcortical_stripe_view,
    )
    save_subcortical_label_volume(
        volume,
        affine,
        network_names,
        network_colors,
        volume_path,
        label_list_path,
        args.wb_command,
    )

    retained_rows = [row for row in cluster_rows if row["retained"]]
    summary = {
        "mode": args.mode,
        "criterion": {
            "secondary_probability_min": args.secondary_prob_min,
            "secondary_to_top_ratio_min": args.secondary_to_top_ratio_min,
            "max_networks_per_zone": args.max_networks_per_zone,
            "min_cortical_cluster_vertices": args.min_cortical_cluster_vertices,
            "min_subcortical_cluster_voxels": args.min_subcortical_cluster_voxels,
        },
        "counts": {
            "raw_multi_network_grayordinates": int(np.sum(support_code > 0)),
            "raw_clusters": len(cluster_rows),
            "retained_clusters": len(retained_rows),
            "retained_cortical_clusters": sum(row["domain"] == "cortex" for row in retained_rows),
            "retained_subcortical_clusters": sum(row["domain"] == "subcortex" for row in retained_rows),
            "retained_grayordinates": int(retained_mask.sum()),
            "retained_two_network_grayordinates": int(np.sum(retained_mask & (retained_count == 2))),
            "retained_three_network_grayordinates": int(np.sum(retained_mask & (retained_count == 3))),
            "eligible_parcels": len(parcel_support) if args.mode == "parcel_gated" else None,
        },
        "outputs": {
            "network_stripes_dlabel": str(dlabel_path.resolve()),
            "diagnostics_dscalar": str(dscalar_path.resolve()),
            "cluster_summary_csv": str(csv_path.resolve()),
            "subcortical_network_stripes_volume": str(volume_path.resolve()),
        },
    }
    with json_path.open("w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
    print(f"[pfm-vertex-stripes] {args.mode}: {summary['counts']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
