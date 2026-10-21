import tempfile
import unittest
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from mefmri_qc_motion_grayplots import RunInfo, load_physio


class MotionGrayplotPhysioTest(unittest.TestCase):
    def test_scanner_ticks_align_resp_and_pulse_to_frames(self):
        with tempfile.TemporaryDirectory() as folder:
            subject = Path(folder)
            physio = subject / 'physio/unprocessed/rest/session_2/run_3'
            physio.mkdir(parents=True)
            (physio / 'scan_Info.log').write_text(
                'VOLUME SLICE ACQ_START_TICS ACQ_FINISH_TICS ECHO\n'
                '0 0 100 120 0\n1 0 180 200 0\n'
            )
            for channel in ('RESP', 'PULS'):
                (physio / f'scan_{channel}.log').write_text(
                    'ACQ_TIME_TICS CHANNEL VALUE SIGNAL\n'
                    + '\n'.join(f'{tick} {channel} {tick % 17}' for tick in range(80, 221, 5))
                    + '\n'
                )
            run = RunInfo(2, 3, subject / 'func/rest/session_2/run_3')
            values = load_physio(subject, 'rest', run, 11)
            self.assertEqual(set(values), {'respiration', 'pulse'})
            for x, y in values.values():
                self.assertGreaterEqual(x.min(), 0)
                self.assertLessEqual(x.max(), 10)
                self.assertTrue(np.isfinite(y).all())
                self.assertAlmostEqual(float(np.std(y)), 1.0, places=6)


if __name__ == '__main__':
    unittest.main()
