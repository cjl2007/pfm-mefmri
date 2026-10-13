#!/usr/bin/env python3
"""Focused tests for RidgePFM review outputs."""

from __future__ import annotations

import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.cifti2.cifti2_axes import BrainModelAxis, LabelAxis, SeriesAxis


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from pfm_ridge_fusion_review import make_areal_probability_maps, rewrite_ridge_scene_template, trim_review_scene_to_two_tabs  # noqa: E402


def dense_axis(n_gray: int) -> BrainModelAxis:
    return BrainModelAxis.from_mask(np.ones(n_gray, dtype=bool), name="CORTEX_LEFT")


class RidgePFMReviewTests(unittest.TestCase):
    def test_areal_probability_maps_are_parcel_averages(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            axis = dense_axis(4)
            labels = {0: ("???", (0, 0, 0, 0)), 1: ("A", (1, 0, 0, 1)), 2: ("B", (0, 1, 0, 1))}
            parcel_hdr = nib.Cifti2Header.from_axes((LabelAxis(["parcels"], labels), axis))
            parcel_img = nib.Cifti2Image(np.asarray([[1, 1, 2, 2]], dtype=np.float32), parcel_hdr)
            parcel_path = tmp / "parcels.dlabel.nii"
            nib.save(parcel_img, str(parcel_path))

            prob = np.asarray(
                [
                    [0.2, 0.4, 0.8, 0.6],
                    [0.8, 0.6, 0.2, 0.4],
                ],
                dtype=np.float32,
            )
            prob_hdr = nib.Cifti2Header.from_axes((SeriesAxis(0, 1, 2), axis))
            prob_img = nib.Cifti2Image(prob, prob_hdr)
            prob_path = tmp / "prob.dtseries.nii"
            out_path = tmp / "areal_prob.dtseries.nii"
            nib.save(prob_img, str(prob_path))

            make_areal_probability_maps(parcel_path, prob_path, out_path)
            out = np.asanyarray(nib.load(str(out_path)).dataobj)
            expected = np.asarray(
                [
                    [0.3, 0.3, 0.7, 0.7],
                    [0.7, 0.7, 0.3, 0.3],
                ],
                dtype=np.float32,
            )
            self.assertTrue(np.allclose(out, expected))

    def test_split_scenes_use_only_requested_panels_and_borders(self) -> None:
        template = (
            '<ObjectArray Type="class" Name="allCaretDataFiles_V2" Length="1">\n'
            '                            </ObjectArray>\n'
            '                            <Object Type="class" Class="SpecFile" Name="specFile" Version="1">\n'
            '                                <ObjectArray Type="class" Name="dataFilesArray" Length="1">\n'
            '                                </ObjectArray>\n'
            '                            </Object>\n'
            '                            <Object Type="class" Class="CaretDataFile" Name="m_sceneAnnotationFile" Version="1"/>\n'
            '                            <Object Type="class" Class="Overlay" Name="m_overlays" Version="1">'
            '<Object Type="boolean" Name="m_enabled">true</Object>'
            '<Object Type="pathName" Name="selectedMapFileNameWithPath">old1</Object>'
            '<Object Type="string" Name="selectedMapFile">old1</Object>'
            '<Object Type="string" Name="selectedMapName">old</Object>'
            '<Object Type="integer" Name="selectedMapIndex">9</Object></Object>'
            '                            <Object Type="class" Class="Overlay" Name="m_overlays" Version="1">'
            '<Object Type="boolean" Name="m_enabled">true</Object>'
            '<Object Type="pathName" Name="selectedMapFileNameWithPath">old2</Object>'
            '<Object Type="string" Name="selectedMapFile">old2</Object>'
            '<Object Type="string" Name="selectedMapName">old</Object>'
            '<Object Type="integer" Name="selectedMapIndex">9</Object></Object>'
            '                            <Object Type="class" Class="Overlay" Name="m_overlays" Version="1">'
            '<Object Type="boolean" Name="m_enabled">true</Object>'
            '<Object Type="pathName" Name="selectedMapFileNameWithPath">old3</Object>'
            '<Object Type="string" Name="selectedMapFile">old3</Object>'
            '<Object Type="string" Name="selectedMapName">old</Object>'
            '<Object Type="integer" Name="selectedMapIndex">9</Object></Object>'
            '                            <Object Type="class" Class="Overlay" Name="m_overlays" Version="1">'
            '<Object Type="boolean" Name="m_enabled">true</Object>'
            '<Object Type="pathName" Name="selectedMapFileNameWithPath">old4</Object>'
            '<Object Type="string" Name="selectedMapFile">old4</Object>'
            '<Object Type="string" Name="selectedMapName">old</Object>'
            '<Object Type="integer" Name="selectedMapIndex">9</Object></Object>'
            '                            <Object Type="class" Class="DisplayPropertiesBorders" Name="displayPropertiesBorders" Version="1">\n'
            '                                <ObjectMap Type="boolean" Name="m_displayStatusInTab"><Value Key="0">false</Value></ObjectMap>\n'
            '                                <ObjectArray Type="boolean" Name="m_displayStatusInDisplayGroup" Length="1"><Element Index="0">false</Element></ObjectArray>\n'
            '                                <ObjectMap Type="enumeratedType" Name="m_coloringTypeInTab"><Value Key="0">CUSTOM</Value></ObjectMap>\n'
            '                                <ObjectArray Type="enumeratedType" Name="m_coloringTypeInDisplayGroup" Length="1"><Element Index="0">CUSTOM</Element></ObjectArray>\n'
            '                                <ObjectMap Type="enumeratedType" Name="m_drawingTypeInTab"><Value Key="0">DRAW_AS_POINTS_SPHERES</Value></ObjectMap>\n'
            '                                <ObjectArray Type="enumeratedType" Name="m_drawingTypeInDisplayGroup" Length="1"><Element Index="0">DRAW_AS_POINTS_SPHERES</Element></ObjectArray>\n'
            '                                <ObjectMap Type="enumeratedType" Name="m_standardColorTypeInTab"><Value Key="0">RED</Value></ObjectMap>\n'
            '                                <ObjectArray Type="enumeratedType" Name="m_standardColorTypeInDisplayGroup" Length="1"><Element Index="0">RED</Element></ObjectArray>\n'
            '                            </Object>\n'
        )
        vertex = rewrite_ridge_scene_template(
            template,
            subject="S1",
            scene_name="S1_vertex",
            ridge_dir=Path("/tmp"),
            left_surf="/tmp/S1.L.inflated.32k_fs_LR.surf.gii",
            right_surf="/tmp/S1.R.inflated.32k_fs_LR.surf.gii",
            wta_file="RidgePFM_VTX.dlabel.nii",
            prob_file="RidgePFM_VTX_ProbMaps.dtseries.nii",
            border_files=[
                ("RidgePFM_VTX.L.border", "CORTEX_LEFT"),
                ("RidgePFM_VTX.R.border", "CORTEX_RIGHT"),
            ],
        )
        areal = rewrite_ridge_scene_template(
            template,
            subject="S1",
            scene_name="S1_areal",
            ridge_dir=Path("/tmp"),
            left_surf="/tmp/S1.L.inflated.32k_fs_LR.surf.gii",
            right_surf="/tmp/S1.R.inflated.32k_fs_LR.surf.gii",
            wta_file="RidgePFM_VTX+ArealParcellation.dlabel.nii",
            prob_file="RidgePFM_VTX+ArealParcellation_ProbMaps.dtseries.nii",
            border_files=[
                ("RidgePFM_VTX+ArealParcellation.L.border", "CORTEX_LEFT"),
                ("RidgePFM_VTX+ArealParcellation.R.border", "CORTEX_RIGHT"),
            ],
        )
        for name in [
            "RidgePFM_VTX.dlabel.nii",
            "RidgePFM_VTX_ProbMaps.dtseries.nii",
            "RidgePFM_VTX.L.border",
            "RidgePFM_VTX.R.border",
        ]:
            self.assertIn(name, vertex)
        for name in [
            "RidgePFM_VTX+ArealParcellation.dlabel.nii",
            "RidgePFM_VTX+ArealParcellation_ProbMaps.dtseries.nii",
            "RidgePFM_VTX+ArealParcellation.L.border",
            "RidgePFM_VTX+ArealParcellation.R.border",
        ]:
            self.assertIn(name, areal)
        self.assertNotIn("RidgePFM_VTX+ArealParcellation.L.border", vertex)
        self.assertNotIn("RidgePFM_VTX.L.border", areal)
        self.assertGreaterEqual(vertex.count("RidgePFM_VTX.dlabel.nii"), 2)
        self.assertGreaterEqual(vertex.count("RidgePFM_VTX_ProbMaps.dtseries.nii"), 2)
        self.assertGreaterEqual(areal.count("RidgePFM_VTX+ArealParcellation.dlabel.nii"), 2)
        self.assertGreaterEqual(areal.count("RidgePFM_VTX+ArealParcellation_ProbMaps.dtseries.nii"), 2)
        for out in (vertex, areal):
            self.assertIn("DRAW_AS_LINES", out)
            self.assertIn("FEATURE_COLORING_TYPE_STANDARD_COLOR", out)
            self.assertIn(">BLACK<", out)

    def test_scene_trim_removes_third_and_fourth_tabs(self) -> None:
        scene = """<SceneFile Version="3">
<ObjectArray Type="class" Name="m_browserTabs" Length="4">
<Element Index="0"><Object Type="class" Class="BrowserTabContent" Name="m_browserTabs" Version="7"/></Element>
<Element Index="1"><Object Type="class" Class="BrowserTabContent" Name="m_browserTabs" Version="7"/></Element>
<Element Index="2"><Object Type="class" Class="BrowserTabContent" Name="m_browserTabs" Version="7"/></Element>
<Element Index="3"><Object Type="class" Class="BrowserTabContent" Name="m_browserTabs" Version="7"/></Element>
</ObjectArray>
<Object Type="integer" Name="m_sceneSelectedTabIndex">3</Object>
<Object Type="string" Name="m_customTileTabsConfigurationLatest">
&lt;TileTabsConfigurationTwo Version=&quot;2&quot;&gt;
    &lt;Columns&gt;
        &lt;Element ContentType=&quot;TAB&quot;/&gt;
        &lt;Element ContentType=&quot;TAB&quot;/&gt;
        &lt;Element ContentType=&quot;TAB&quot;/&gt;
        &lt;Element ContentType=&quot;TAB&quot;/&gt;
    &lt;/Columns&gt;
    &lt;Rows&gt;&lt;Element ContentType=&quot;TAB&quot;/&gt;&lt;/Rows&gt;
&lt;/TileTabsConfigurationTwo&gt;</Object>
<ObjectArray Type="integer" Name="m_sceneTabIndices" Length="4">
<Element Index="0">0</Element><Element Index="1">1</Element><Element Index="2">2</Element><Element Index="3">3</Element>
</ObjectArray>
<ObjectArray Type="integer" Name="tabIndices" Length="4">
<Element Index="0">0</Element><Element Index="1">1</Element><Element Index="2">2</Element><Element Index="3">3</Element>
</ObjectArray>
</SceneFile>"""
        out = trim_review_scene_to_two_tabs(scene)
        root = ET.fromstring(out)
        arrays = {array.get("Name"): array for array in root.iter("ObjectArray")}
        self.assertEqual(arrays["m_browserTabs"].get("Length"), "2")
        self.assertEqual(arrays["m_sceneTabIndices"].get("Length"), "2")
        self.assertEqual(arrays["tabIndices"].get("Length"), "2")
        self.assertNotIn('Index="2"', out)
        self.assertNotIn('Index="3"', out)
        self.assertIn("PercentStretch=\"50.00\"", out)
        selected = [obj for obj in root.iter("Object") if obj.get("Name") == "m_sceneSelectedTabIndex"]
        self.assertEqual(selected[0].text, "0")


if __name__ == "__main__":
    unittest.main()
