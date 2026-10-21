#!/usr/bin/env python3
"""Create RidgePFM areal probability maps and review scenes."""

from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import nibabel as nib
import numpy as np

from pfm_infomap_labeler import add_review_scene_border_files, relpath_or_name, retarget_review_scene_primary_overlays


def bool_int(value: str) -> int:
    ivalue = int(value)
    if ivalue not in (0, 1):
        raise argparse.ArgumentTypeError("expected 0 or 1")
    return ivalue


def save_dtseries(ref_img: nib.Cifti2Image, maps: np.ndarray, out_path: Path) -> None:
    from nibabel.cifti2.cifti2_axes import SeriesAxis

    arr = np.asarray(maps, dtype=np.float32)
    if arr.ndim == 1:
        arr = arr[None, :]
    axis = ref_img.header.get_axis(1)
    series = SeriesAxis(start=0.0, step=1.0, size=arr.shape[0], unit="SECOND")
    hdr = nib.Cifti2Header.from_axes((series, axis))
    nib.save(nib.Cifti2Image(arr, hdr, nifti_header=ref_img.nifti_header), str(out_path))


def make_areal_probability_maps(parcels_dlabel: Path, prob_maps: Path, out_path: Path) -> None:
    parcel_img = nib.load(str(parcels_dlabel))
    prob_img = nib.load(str(prob_maps))
    parcel_data = np.asanyarray(parcel_img.dataobj)
    prob_data = np.asanyarray(prob_img.dataobj).astype(np.float32, copy=False)

    if parcel_data.ndim != 2 or parcel_data.shape[0] < 1:
        raise ValueError(f"Expected dlabel with at least one map: {parcels_dlabel}")
    parcels = np.rint(parcel_data[0]).astype(np.int32)
    if prob_data.ndim != 2:
        raise ValueError(f"Expected 2D probability dtseries: {prob_maps}")
    if prob_data.shape[1] != parcels.size:
        raise ValueError(
            f"Probability maps grayordinates ({prob_data.shape[1]}) do not match parcels ({parcels.size})"
        )

    out = np.zeros_like(prob_data, dtype=np.float32)
    for parcel_id in sorted(int(x) for x in np.unique(parcels) if int(x) > 0):
        mask = parcels == parcel_id
        if not np.any(mask):
            continue
        out[:, mask] = np.nanmean(prob_data[:, mask], axis=1, keepdims=True)
    save_dtseries(prob_img, out, out_path)


def _selection_map_xml(name: str, values: list[bool], indent: str) -> str:
    lines = [f'{indent}<ObjectMap Type="enumeratedType" Name="{name}">']
    for idx, value in enumerate(values):
        lines.append(f'{indent}    <Value Key="{idx}">{"SELECTED" if value else "UNSELECTED"}</Value>')
    lines.append(f"{indent}</ObjectMap>")
    return "\n".join(lines) + "\n"


def _selection_array_xml(name: str, values: list[bool], indent: str) -> str:
    lines = [f'{indent}<ObjectArray Type="enumeratedType" Name="{name}" Length="{len(values)}">']
    for idx, value in enumerate(values):
        lines.append(f'{indent}    <Element Index="{idx}">{"SELECTED" if value else "UNSELECTED"}</Element>')
    lines.append(f"{indent}</ObjectArray>")
    return "\n".join(lines) + "\n"


def _bool_map_xml(name: str, values: list[bool], indent: str) -> str:
    lines = [f'{indent}<ObjectMap Type="boolean" Name="{name}">']
    for idx, value in enumerate(values):
        lines.append(f'{indent}    <Value Key="{idx}">{"true" if value else "false"}</Value>')
    lines.append(f"{indent}</ObjectMap>")
    return "\n".join(lines) + "\n"


def _bool_array_xml(name: str, values: list[bool], indent: str) -> str:
    lines = [f'{indent}<ObjectArray Type="boolean" Name="{name}" Length="{len(values)}">']
    for idx, value in enumerate(values):
        lines.append(f'{indent}    <Element Index="{idx}">{"true" if value else "false"}</Element>')
    lines.append(f"{indent}</ObjectArray>")
    return "\n".join(lines) + "\n"


def _border_display_helper_xml(selected_tabs: list[bool], selected_groups: list[bool], indent: str) -> str:
    collapsed_tabs = [False] * len(selected_tabs)
    collapsed_groups = [False] * len(selected_groups)
    return (
        f'{indent}<Object Type="class" Class="DisplayGroupAndTabItemHelper" Name="m_displayGroupAndTabItemHelper" Version="1">\n'
        + _selection_map_xml("m_selectedInTab", selected_tabs, indent + "    ")
        + _bool_map_xml("m_expandedStatusInTab", collapsed_tabs, indent + "    ")
        + _selection_array_xml("m_selectedInDisplayGroup", selected_groups, indent + "    ")
        + _bool_array_xml("m_expandedStatusInDisplayGroup", collapsed_groups, indent + "    ")
        + f"{indent}</Object>\n"
    )


def _replace_enum_map_values(block: str, name: str, values: list[str]) -> str:
    pattern = re.compile(
        rf'(<ObjectMap Type="enumeratedType" Name="{re.escape(name)}"[^>]*>)(.*?)(</ObjectMap>)',
        re.DOTALL,
    )

    def repl(match: re.Match) -> str:
        indent = re.search(r"\n([ \t]*)<Value\b", match.group(2))
        value_indent = indent.group(1) if indent else "    "
        body = "".join(
            f'\n{value_indent}<Value Key="{idx}">{value}</Value>'
            for idx, value in enumerate(values)
        )
        return match.group(1) + body + "\n" + match.group(3)

    out, count = pattern.subn(repl, block, count=1)
    if count > 0:
        return out
    lines = ['    <ObjectMap Type="enumeratedType" Name="m_displayGroup">']
    for idx, value in enumerate(values):
        lines.append(f'        <Value Key="{idx}">{value}</Value>')
    lines.append("    </ObjectMap>")
    insert = "\n".join(lines) + "\n"
    return re.sub(
        r'(<Object Type="class" Class="DisplayPropertiesBorders" Name="displayPropertiesBorders" Version="1">\n)',
        r"\1" + insert,
        block,
        count=1,
    )


def _replace_bool_map_values(block: str, name: str, values: list[bool]) -> str:
    pattern = re.compile(
        rf'(<ObjectMap Type="boolean" Name="{re.escape(name)}"[^>]*>)(.*?)(</ObjectMap>)',
        re.DOTALL,
    )

    def repl(match: re.Match) -> str:
        indent = re.search(r"\n([ \t]*)<Value\b", match.group(2))
        value_indent = indent.group(1) if indent else "    "
        body = "".join(
            f'\n{value_indent}<Value Key="{idx}">{"true" if value else "false"}</Value>'
            for idx, value in enumerate(values)
        )
        return match.group(1) + body + "\n" + match.group(3)

    return pattern.sub(repl, block, count=1)


def _replace_bool_array_values(block: str, name: str, values: list[bool]) -> str:
    pattern = re.compile(
        rf'(<ObjectArray Type="boolean" Name="{re.escape(name)}"[^>]*>)(.*?)(</ObjectArray>)',
        re.DOTALL,
    )

    def repl(match: re.Match) -> str:
        indent = re.search(r"\n([ \t]*)<Element\b", match.group(2))
        value_indent = indent.group(1) if indent else "    "
        body = "".join(
            f'\n{value_indent}<Element Index="{idx}">{"true" if value else "false"}</Element>'
            for idx, value in enumerate(values)
        )
        opening = re.sub(r'Length="\d+"', f'Length="{len(values)}"', match.group(1))
        return opening + body + "\n" + match.group(3)

    return pattern.sub(repl, block, count=1)


def set_ridge_border_display_groups(
    text: str,
    *,
    network_border_names: list[str],
    areal_border_names: list[str],
) -> str:
    """Use tab-specific border selection so vertexwise and areal borders do not mix."""
    border_re = re.compile(
        r'(<Object Type="class" Class="DisplayPropertiesBorders" Name="displayPropertiesBorders" Version="1">.*?</Object>)',
        re.DOTALL,
    )

    def set_display_group_block(match: re.Match) -> str:
        block = match.group(1)
        block = _replace_enum_map_values(
            block,
            "m_displayGroup",
            ["DISPLAY_GROUP_TAB", "DISPLAY_GROUP_TAB", "DISPLAY_GROUP_TAB", "DISPLAY_GROUP_TAB"],
        )
        block = _replace_bool_map_values(block, "m_displayStatusInTab", [True, True, True, True])
        block = _replace_bool_array_values(block, "m_displayStatusInDisplayGroup", [False, False, False, False, False])
        return block

    text = border_re.sub(set_display_group_block, text, count=1)

    def set_file_helper(filename: str, selected_tabs: list[bool], selected_groups: list[bool], text: str) -> str:
        pattern = re.compile(
            rf'(<Object Type="class" Class="CaretDataFile" Name="{re.escape(filename)}" Version="1">)(.*?)(</Object>)',
            re.DOTALL,
        )

        def repl(match: re.Match) -> str:
            body = re.sub(
                r'\s*<Object Type="class" Class="DisplayGroupAndTabItemHelper" Name="m_displayGroupAndTabItemHelper" Version="1">.*?</Object>',
                "",
                match.group(2),
                flags=re.DOTALL,
            )
            indent_match = re.search(r"\n([ \t]*)<Object Type=\"pathName\"", body)
            helper_indent = indent_match.group(1) if indent_match else "    "
            helper = _border_display_helper_xml(selected_tabs, selected_groups, helper_indent)
            return match.group(1) + "\n" + helper + body + match.group(3)

        return pattern.sub(repl, text, count=1)

    for filename in network_border_names:
        text = set_file_helper(filename, [True, True, False, False], [False, False, False, False, False], text)
    for filename in areal_border_names:
        text = set_file_helper(filename, [False, False, True, True], [False, False, False, False, False], text)
    return text


def _set_object_text(root: ET.Element, name: str, value: str) -> None:
    for obj in root.iter("Object"):
        if obj.get("Name") == name:
            obj.text = value


def _trim_integer_array(root: ET.Element, name: str, values: list[int]) -> None:
    for array in root.iter("ObjectArray"):
        if array.get("Name") != name or array.get("Type") != "integer":
            continue
        for child in list(array):
            array.remove(child)
        array.set("Length", str(len(values)))
        for idx, value in enumerate(values):
            elem = ET.SubElement(array, "Element", {"Index": str(idx)})
            elem.text = str(value)


def _trim_browser_tabs(root: ET.Element, keep_count: int = 2) -> None:
    for array in root.iter("ObjectArray"):
        if array.get("Name") != "m_browserTabs":
            continue
        children = list(array)
        if len(children) <= keep_count:
            array.set("Length", str(len(children)))
            continue
        for child in children[keep_count:]:
            array.remove(child)
        array.set("Length", str(keep_count))


def _tile_configuration_1x2_v2() -> str:
    return """\
<TileTabsConfigurationTwo Version="2">
    <Name></Name>
    <UniqueIdentifier>{bb636b14-73ee-4dc7-b7ae-098cc7b5f753}</UniqueIdentifier>
    <CenteringCorrection>false</CenteringCorrection>
    <Columns>
        <Element ContentType="TAB" StretchType="WEIGHT" PercentStretch="50.00" WeightStretch="1.00"/>
        <Element ContentType="TAB" StretchType="WEIGHT" PercentStretch="50.00" WeightStretch="1.00"/>
    </Columns>
    <Rows>
        <Element ContentType="TAB" StretchType="WEIGHT" PercentStretch="20.00" WeightStretch="1.00"/>
    </Rows>
</TileTabsConfigurationTwo>"""


def _tile_configuration_1x2_v1() -> str:
    return """\
<TileTabsConfiguration>
    <Version Number="1"/>
    <Name></Name>
    <UniqueIdentifier>{bb636b14-73ee-4dc7-b7ae-098cc7b5f753}</UniqueIdentifier>
    <RowStretchFactors SelectedColumnCount="1" TotalColumnCount="1">1</RowStretchFactors>
    <ColumnStretchFactors SelectedRowCount="2" TotalRowCount="2">1 1</ColumnStretchFactors>
</TileTabsConfiguration>"""


def trim_review_scene_to_two_tabs(text: str) -> str:
    """Drop inherited tabs 3/4 so split review scenes cannot mix border families."""
    try:
        root = ET.fromstring(text)
    except ET.ParseError:
        return text

    _trim_browser_tabs(root, keep_count=2)
    _trim_integer_array(root, "m_sceneTabIndices", [0, 1])
    _trim_integer_array(root, "tabIndices", [0, 1])
    _set_object_text(root, "m_sceneSelectedTabIndex", "0")
    _set_object_text(root, "selectedTabIndex", "0")
    _set_object_text(root, "m_customTileTabsConfigurationLatest", "\n" + _tile_configuration_1x2_v2())
    _set_object_text(root, "m_customTileTabsConfiguration", "\n" + _tile_configuration_1x2_v1())
    _set_object_text(root, "m_sceneTileTabsConfiguration", "\n" + _tile_configuration_1x2_v1())
    return ET.tostring(root, encoding="unicode")


def rewrite_ridge_scene_template(
    text: str,
    *,
    subject: str,
    scene_name: str,
    ridge_dir: Path,
    left_surf: str,
    right_surf: str,
    wta_file: str,
    prob_file: str,
    border_files: list[tuple[str, str]],
) -> str:
    replacements = {
        "__SCENE_NAME__": scene_name,
        "__SUBJECT__": subject,
        "HC05_manual_labeling": scene_name,
        "Bipartite_PhysicalCommunities+AlgorithmicLabeling.dlabel.nii": wta_file,
        "Bipartite_PhysicalCommunities+AlgorithmicLabeling_InfoMapCommunities.dlabel.nii": wta_file,
        "Bipartite_PhysicalCommunities+AlgorithmicLabeling_FC_WholeBrain.dtseries.nii": prob_file,
        "Bipartite_PhysicalCommunities+AlgorithmicLabeling_FC_btwn_InfoMapCommunities.dtseries.nii": prob_file,
        "Bipartite_PhysicalCommunities+AlgorithmicLabeling_FC_WholeBrain.dynconn.nii": prob_file,
        "Bipartite_PhysicalCommunities+AlgorithmicLabeling_FC_btwn_InfoMapCommunities.dynconn.nii": prob_file,
        "../../../../../../../outputs/data/Priors/Priors_FC.dtseries.nii": prob_file,
        "../../../../../../../outputs/data/Priors/Priors_FC.dynconn.nii": prob_file,
    }
    for old, new in replacements.items():
        text = text.replace(old, new)

    text = retarget_review_scene_primary_overlays(
        text,
        [wta_file, prob_file],
    )

    left_path = Path(left_surf) if left_surf else None
    right_path = Path(right_surf) if right_surf else None
    surf_dir = left_path.parent if left_path else right_path.parent if right_path else None
    if subject and surf_dir is not None:
        for hemi, source_path in (("L", left_path), ("R", right_path)):
            for surf in ("inflated", "very_inflated"):
                basename = f"{subject}.{hemi}.{surf}.32k_fs_LR.surf.gii"
                rel = relpath_or_name(surf_dir / basename, ridge_dir)
                text = text.replace(f"HC05.{hemi}.{surf}.32k_fs_LR.surf.gii", basename)
                text = text.replace(f"../../../../anat/T1w/fsaverage_LR32k/{basename}", rel)
                text = text.replace(f"../../../../anat/T1w/fsaverage_LR32k/HC05.{hemi}.{surf}.32k_fs_LR.surf.gii", rel)
    if subject:
        text = text.replace("HC05", subject)
    text = add_review_scene_border_files(text, border_files)
    text = trim_review_scene_to_two_tabs(text)
    return text


def write_one_review_scene(
    args: argparse.Namespace,
    *,
    scene_outfile: str,
    scene_name_suffix: str,
    wta_path: Path,
    prob_path: Path,
    border_paths: list[Path],
) -> Path:
    ridge_dir = Path(args.outdir)
    template = Path(args.scene_template)
    wta_file = wta_path.name
    prob_file = prob_path.name
    border_files = [
        (border_paths[0].name, "CORTEX_LEFT"),
        (border_paths[1].name, "CORTEX_RIGHT"),
    ]
    required = [
        wta_path,
        prob_path,
        border_paths[0],
        border_paths[1],
        template,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing RidgePFM scene inputs: " + ", ".join(missing))

    scene_name = args.scene_name or f"{args.subject}_ridgepfm_{scene_name_suffix}"
    text = template.read_text(encoding="utf-8", errors="surrogateescape")
    text = rewrite_ridge_scene_template(
        text,
        subject=args.subject,
        scene_name=scene_name,
        ridge_dir=ridge_dir,
        left_surf=args.left_surf,
        right_surf=args.right_surf,
        wta_file=wta_file,
        prob_file=prob_file,
        border_files=border_files,
    )
    out_path = ridge_dir / scene_outfile
    out_path.write_text(text, encoding="utf-8", errors="surrogateescape")
    return out_path


def write_review_scenes(args: argparse.Namespace, areal_score_maps: Path) -> list[Path]:
    return [
        write_one_review_scene(
            args,
            scene_outfile=args.vertex_scene_outfile,
            scene_name_suffix="vertexwise_review",
            wta_path=Path(args.network_dlabel),
            prob_path=Path(args.score_maps),
            border_paths=[Path(args.network_left_border), Path(args.network_right_border)],
        ),
        write_one_review_scene(
            args,
            scene_outfile=args.areal_scene_outfile,
            scene_name_suffix="areal_parcellation_review",
            wta_path=Path(args.areal_dlabel),
            prob_path=areal_score_maps,
            border_paths=[Path(args.areal_left_border), Path(args.areal_right_border)],
        ),
    ]


def main() -> int:
    ap = argparse.ArgumentParser(description="Create RidgePFM review outputs")
    ap.add_argument("--network-dlabel", required=True)
    ap.add_argument(
        "--score-maps",
        "--prob-maps",
        dest="score_maps",
        required=True,
        help="normalized fused-score maps; --prob-maps is a deprecated alias",
    )
    ap.add_argument("--areal-dlabel", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--outfile-prefix", default="RidgePFM_VTX+ArealParcellation")
    ap.add_argument("--write-scene", type=bool_int, default=1)
    ap.add_argument("--scene-template", required=True)
    ap.add_argument("--scene-outfile", default="")
    ap.add_argument("--vertex-scene-outfile", default="RidgePFM_Vertexwise_Review.scene")
    ap.add_argument("--areal-scene-outfile", default="RidgePFM_ArealParcellation_Review.scene")
    ap.add_argument("--scene-name", default="")
    ap.add_argument("--subject", default="")
    ap.add_argument("--left-surf", required=True)
    ap.add_argument("--right-surf", required=True)
    ap.add_argument("--network-left-border", required=True)
    ap.add_argument("--network-right-border", required=True)
    ap.add_argument("--areal-left-border", required=True)
    ap.add_argument("--areal-right-border", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    areal_score_maps = outdir / f"{args.outfile_prefix}_FusedScoreSoftmax.dtseries.nii"
    make_areal_probability_maps(Path(args.areal_dlabel), Path(args.score_maps), areal_score_maps)
    print(f"[ridge-review] wrote {areal_score_maps}")
    if int(args.write_scene) == 1:
        if args.scene_outfile:
            args.vertex_scene_outfile = args.scene_outfile
        for scene_path in write_review_scenes(args, areal_score_maps):
            print(f"[ridge-review] wrote {scene_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
