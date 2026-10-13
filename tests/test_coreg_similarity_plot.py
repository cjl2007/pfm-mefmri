"""Keep run-similarity heatmaps readable without altering the matrix or scale."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from PIL import Image

ROOT = Path(os.environ.get('MEFMRI_TEST_ROOT', Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / 'lib'))
import mefmri_qc_coreg_sdc_report as report


class CoregSimilarityPlotTests(unittest.TestCase):
    def check_plot(self, n):
        sim = np.full((n, n), 0.98)
        np.fill_diagonal(sim, 1.0)
        sim[0, -1] = sim[-1, 0] = 0.965
        before = sim.copy()
        labels = [f'S{i // 2 + 1}R{i % 2 + 1}' for i in range(n)]
        figures_before = report.plt.get_fignums()
        fig, ax = report.plt.subplots(figsize=(5.2, 4.5), facecolor='white', constrained_layout=True)
        self.addCleanup(report.plt.close, fig)
        with tempfile.TemporaryDirectory(prefix='coreg-similarity-') as temp:
            output = Path(temp) / 'matrix.png'
            with patch.object(report.plt, 'subplots', return_value=(fig, ax)):
                report.write_run_similarity_matrix(sim, labels, output)
            self.assertTrue(output.is_file())
            with Image.open(output) as image:
                self.assertEqual(image.size, (780, 675))
            self.assertEqual(len(ax.texts), 0, 'Cells must not contain numeric text annotations')
            self.assertEqual([tick.get_text() for tick in ax.get_xticklabels()], labels)
            self.assertEqual([tick.get_text() for tick in ax.get_yticklabels()], labels)
            self.assertEqual(len(ax.images), 1)
            np.testing.assert_array_equal(ax.images[0].get_array(), before)
            np.testing.assert_array_equal(sim, before)
            self.assertEqual(ax.images[0].get_cmap().name, 'viridis')
            self.assertEqual(ax.images[0].get_clim(), (float(before.min()), 1.0))
            self.assertEqual(len(fig.axes), 2)
            self.assertEqual(fig.axes[1].get_ylabel(), 'Pearson r')
            self.assertEqual(report.plt.get_fignums(), figures_before)

    def test_two_run_heatmap_has_labels_and_colorbar_but_no_cell_text(self):
        self.check_plot(2)

    def test_twenty_six_run_heatmap_preserves_values_without_cell_text(self):
        self.check_plot(26)


if __name__ == '__main__':
    unittest.main()
