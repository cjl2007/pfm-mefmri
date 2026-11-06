"""Validate confidence eligibility before expensive PFM preparation and fitting."""
import argparse
import nibabel as nib
import numpy as np
from pfm_ridge_confidence import load_sessions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--in-cifti', required=True)
    parser.add_argument('--scan-idx', required=True)
    parser.add_argument('--minimum-frames', type=int, default=30)
    args = parser.parse_args()
    if args.minimum_frames < 1:
        parser.error('--minimum-frames must be positive')
    try:
        sessions = load_sessions(args.scan_idx, nib.load(args.in_cifti).shape[0])
        ids, counts = np.unique(sessions, return_counts=True)
        eligible = ids[counts >= args.minimum_frames]
        if len(eligible) < 4:
            parser.error(f'PFM_RF_CONFIDENCE_ENABLE=1 requires at least four eligible sessions '
                         f'(two per half); found {len(eligible)} with at least {args.minimum_frames} '
                         'retained frames each. Set PFM_RF_CONFIDENCE_ENABLE=0 to run RidgePFM '
                         'without confidence analysis.')
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    print(f'[pfm] confidence preflight passed: {len(eligible)} eligible sessions')


if __name__ == '__main__':
    main()
