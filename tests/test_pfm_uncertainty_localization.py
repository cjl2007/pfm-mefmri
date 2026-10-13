#!/usr/bin/env python3
"""Focused tests for temporal CV and vertexwise uncertainty localization."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from pfm_parcel_run_split_cv import (  # noqa: E402
    benjamini_hochberg_by_domain,
    candidate_identity_null,
    convex_three_fit,
    convex_two_fit,
    load_raw_run_lengths,
    random_third_folds,
    repeat_sign_flip_max_null,
)
from pfm_vertex_uncertainty_stripes import code_networks, network_code  # noqa: E402
from pfm_infomap_labeler import add_review_scene_border_files, retarget_review_scene_primary_overlays  # noqa: E402


class TemporalCVTests(unittest.TestCase):
    def test_random_thirds_are_disjoint_and_complete(self) -> None:
        run_slices = [slice(0, 31), slice(31, 70)]
        folds = random_third_folds(run_slices, repeats=2, block_length=5, seed=9)
        self.assertEqual(len(folds), 6)
        expected = set(range(70))
        for repeat in (1, 2):
            repeat_folds = [x for x in folds if x["repeat"] == repeat]
            test_sets = [set(x["test_indices"].tolist()) for x in repeat_folds]
            self.assertEqual(set.union(*test_sets), expected)
            self.assertFalse(test_sets[0] & test_sets[1])
            self.assertFalse(test_sets[0] & test_sets[2])
            self.assertFalse(test_sets[1] & test_sets[2])
            for fold in repeat_folds:
                self.assertFalse(set(fold["train_indices"]) & set(fold["test_indices"]))

    def test_two_network_convex_fit_recovers_known_mixture(self) -> None:
        t0 = np.array([1.0, 0.0, 0.0])
        t1 = np.array([0.0, 1.0, 0.0])
        y = 0.25 * t0 + 0.75 * t1
        weights, prediction = convex_two_fit(y, t0, t1)
        self.assertTrue(np.allclose(weights, [0.25, 0.75]))
        self.assertTrue(np.allclose(prediction, y))

    def test_three_network_convex_fit_recovers_known_mixture(self) -> None:
        templates = np.eye(3)
        y = np.array([0.2, 0.3, 0.5])
        weights, prediction = convex_three_fit(y, *templates)
        self.assertTrue(np.allclose(weights, y))
        self.assertTrue(np.allclose(prediction, y))

    def test_repeat_sign_flip_null_is_exact_and_domain_stratified(self) -> None:
        gains = np.asarray(
            [
                [0.4, 0.5, 0.3, 0.4, 0.5, 0.3],
                [0.1, -0.1, 0.0, 0.1, -0.1, 0.0],
                [0.2, 0.3, 0.2, 0.2, 0.3, 0.2],
            ]
        )
        repeats = np.asarray([1, 1, 1, 2, 2, 2])
        domains = np.asarray(["cortex", "cortex", "subcortex"])
        thresholds, p_fwer, n_patterns = repeat_sign_flip_max_null(
            gains,
            repeats,
            domains,
            alpha=0.25,
        )
        self.assertEqual(n_patterns, 4)
        self.assertEqual(set(thresholds), {"cortex", "subcortex"})
        self.assertTrue(np.all((p_fwer >= 0.0) & (p_fwer <= 1.0)))
        self.assertLessEqual(p_fwer[0], p_fwer[1])

    def test_scan_info_supplies_run_lengths(self) -> None:
        import json
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "ScanInfo.json"
            path.write_text(json.dumps([{"n_timepoints": 10}, {"n_timepoints": 12}]))
            self.assertEqual(load_raw_run_lengths(path), [10, 12])

    def test_candidate_identity_null_rewards_specific_generalization(self) -> None:
        repeats = np.repeat(np.arange(1, 9), 3)
        gains = np.zeros((1, repeats.size, 5), dtype=float)
        gains[:, :, 0] = np.nan
        gains[:, :, 1] = 1.0
        selected = np.ones((1, repeats.size), dtype=np.int32)
        excluded = np.zeros((1, repeats.size, 1), dtype=np.int32)
        p_specific, observed = candidate_identity_null(
            gains,
            selected,
            excluded,
            repeats,
            iterations=2000,
            seed=4,
        )
        self.assertAlmostEqual(observed[0], 1.0)
        self.assertLess(p_specific[0], 0.1)

        gains[:, :, 2:] = 1.0
        p_generic, _ = candidate_identity_null(
            gains,
            selected,
            excluded,
            repeats,
            iterations=500,
            seed=4,
        )
        self.assertEqual(p_generic[0], 1.0)

    def test_bh_is_separate_by_domain(self) -> None:
        p = np.asarray([0.01, 0.04, 0.03, 0.20])
        domains = np.asarray(["cortex", "cortex", "subcortex", "subcortex"])
        q = benjamini_hochberg_by_domain(p, domains)
        self.assertTrue(np.allclose(q, [0.02, 0.04, 0.06, 0.20]))


class VertexLocalizationTests(unittest.TestCase):
    def test_network_signature_roundtrip(self) -> None:
        ids = [2, 7, 19]
        self.assertEqual(code_networks(network_code(ids), 20), ids)


class InfomapReviewSceneTests(unittest.TestCase):
    def test_primary_overlay_panels_are_retargeted_in_order(self) -> None:
        overlay = (
            '<Object Type="class" Class="Overlay" Name="m_overlays" Version="1">'
            '<Object Type="boolean" Name="m_enabled">true</Object>'
            '<Object Type="pathName" Name="selectedMapFileNameWithPath">old.dtseries.nii</Object>'
            '<Object Type="string" Name="selectedMapFile">old.dtseries.nii</Object>'
            '<Object Type="string" Name="selectedMapName">20 seconds</Object>'
            '<Object Type="integer" Name="selectedMapIndex">19</Object>'
            "</Object>"
        )
        text = overlay * 5
        out = retarget_review_scene_primary_overlays(
            text,
            ["a.dlabel.nii", "b.dlabel.nii", "c_wholebrain.dtseries.nii", "d_between.dtseries.nii"],
        )
        self.assertEqual(out.count("old.dtseries.nii"), 0)
        self.assertEqual(out.count("<Object Type=\"integer\" Name=\"selectedMapIndex\">0</Object>"), 5)
        first_positions = [
            out.index(name)
            for name in ["a.dlabel.nii", "b.dlabel.nii", "c_wholebrain.dtseries.nii", "d_between.dtseries.nii"]
        ]
        self.assertEqual(first_positions, sorted(first_positions))

    def test_border_files_are_loaded_and_enabled(self) -> None:
        text = (
            '<ObjectArray Type="class" Name="allCaretDataFiles_V2" Length="1">\n'
            '                            </ObjectArray>\n'
            '                            <Object Type="class" Class="SpecFile" Name="specFile" Version="1">\n'
            '                                <ObjectArray Type="class" Name="dataFilesArray" Length="1">\n'
            '                                </ObjectArray>\n'
            '                            </Object>\n'
            '                            <Object Type="class" Class="CaretDataFile" Name="m_sceneAnnotationFile" Version="1"/>\n'
            '                            <Object Type="class" Class="DisplayPropertiesBorders" Name="displayPropertiesBorders" Version="1">\n'
            '                                <ObjectMap Type="boolean" Name="m_displayStatusInTab">\n'
            '                                    <Value Key="0">false</Value>\n'
            '                                    <Value Key="1">false</Value>\n'
            '                                </ObjectMap>\n'
            '                                <ObjectArray Type="boolean" Name="m_displayStatusInDisplayGroup" Length="2">\n'
            '                                    <Element Index="0">false</Element>\n'
            '                                    <Element Index="1">false</Element>\n'
            '                                </ObjectArray>\n'
            '                                <ObjectMap Type="enumeratedType" Name="m_coloringTypeInTab">\n'
            '                                    <Value Key="0">FEATURE_COLORING_TYPE_CUSTOM_COLOR</Value>\n'
            '                                    <Value Key="1">FEATURE_COLORING_TYPE_CUSTOM_COLOR</Value>\n'
            '                                </ObjectMap>\n'
            '                                <ObjectArray Type="enumeratedType" Name="m_coloringTypeInDisplayGroup" Length="2">\n'
            '                                    <Element Index="0">FEATURE_COLORING_TYPE_CUSTOM_COLOR</Element>\n'
            '                                    <Element Index="1">FEATURE_COLORING_TYPE_CUSTOM_COLOR</Element>\n'
            '                                </ObjectArray>\n'
            '                                <ObjectMap Type="enumeratedType" Name="m_drawingTypeInTab">\n'
            '                                    <Value Key="0">DRAW_AS_POINTS_SPHERES</Value>\n'
            '                                    <Value Key="1">DRAW_AS_POINTS_SPHERES</Value>\n'
            '                                </ObjectMap>\n'
            '                                <ObjectArray Type="enumeratedType" Name="m_drawingTypeInDisplayGroup" Length="2">\n'
            '                                    <Element Index="0">DRAW_AS_POINTS_SPHERES</Element>\n'
            '                                    <Element Index="1">DRAW_AS_POINTS_SPHERES</Element>\n'
            '                                </ObjectArray>\n'
            '                                <ObjectMap Type="enumeratedType" Name="m_standardColorTypeInTab">\n'
            '                                    <Value Key="0">RED</Value>\n'
            '                                    <Value Key="1">BLUE</Value>\n'
            '                                </ObjectMap>\n'
            '                                <ObjectArray Type="enumeratedType" Name="m_standardColorTypeInDisplayGroup" Length="2">\n'
            '                                    <Element Index="0">RED</Element>\n'
            '                                    <Element Index="1">BLUE</Element>\n'
            '                                </ObjectArray>\n'
            '                            </Object>\n'
        )
        out = add_review_scene_border_files(text, [("networks.L.border", "CORTEX_LEFT"), ("networks.R.border", "CORTEX_RIGHT")])
        self.assertIn('Name="allCaretDataFiles_V2" Length="3"', out)
        self.assertIn('Name="dataFilesArray" Length="3"', out)
        self.assertIn("<Object Type=\"enumeratedType\" Name=\"dataFileType\">BORDER</Object>", out)
        self.assertIn("<Object Type=\"enumeratedType\" Name=\"structure\">CORTEX_LEFT</Object>", out)
        self.assertIn("<Object Type=\"enumeratedType\" Name=\"structure\">CORTEX_RIGHT</Object>", out)
        self.assertNotIn(">false<", out)
        self.assertNotIn("FEATURE_COLORING_TYPE_CUSTOM_COLOR", out)
        self.assertNotIn("DRAW_AS_POINTS_SPHERES", out)
        self.assertGreaterEqual(out.count("FEATURE_COLORING_TYPE_STANDARD_COLOR"), 4)
        self.assertGreaterEqual(out.count("DRAW_AS_LINES"), 4)
        self.assertGreaterEqual(out.count(">BLACK<"), 4)


if __name__ == "__main__":
    unittest.main()
