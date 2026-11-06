#!/usr/bin/env python3
"""Fail closed on invalid/mismatched FSL timing files. No input files are edited."""
import argparse
import gzip
import math
import struct
from pathlib import Path


def nifti_shape(path):
    path = Path(path)
    opener = gzip.open if path.name.endswith('.gz') else open
    with opener(path, 'rb') as f:
        h = f.read(540)
    if len(h) < 348:
        raise ValueError('Truncated NIfTI header: {}'.format(path))
    endian = next((e for e in ('<', '>') if struct.unpack(e+'i', h[:4])[0] in (348, 540)), None)
    if endian is None:
        raise ValueError('Invalid NIfTI header: {}'.format(path))
    version = struct.unpack(endian+'i', h[:4])[0]
    if version == 348:
        dims = struct.unpack(endian+'8h', h[40:56])
    else:
        if len(h) < 540:
            raise ValueError('Truncated NIfTI-2 header: {}'.format(path))
        dims = struct.unpack(endian+'8q', h[16:80])
    if not 3 <= dims[0] <= 7 or any(d <= 0 for d in dims[1:4]):
        raise ValueError('Invalid NIfTI dimensions: {}'.format(path))
    return tuple(dims[1:dims[0]+1])


def validate_metadata(data, n_slices, source='JSON'):
    times = data.get('SliceTiming')
    if not isinstance(times, list) or not times:
        raise ValueError('{}: missing/empty SliceTiming. Recover verified timing from DICOM; do not infer it from TR or MB alone.'.format(source))
    if len(times) != n_slices:
        raise ValueError('{}: {} SliceTiming entries for {} slices'.format(source, len(times), n_slices))
    if any(isinstance(t, bool) or not isinstance(t, (int, float)) or not math.isfinite(t) for t in times):
        raise ValueError('{}: SliceTiming must contain finite numbers'.format(source))
    tr = data.get('RepetitionTime')
    if isinstance(tr, bool) or not isinstance(tr, (int, float)) or not math.isfinite(tr) or tr <= 0:
        raise ValueError('{}: invalid RepetitionTime'.format(source))
    if any(t < 0 or t >= tr for t in times):
        raise ValueError('{}: SliceTiming must lie in [0, RepetitionTime) seconds'.format(source))
    # The existing headmotion code assumes z slicing and uses --direction=3.
    # Do not silently accept timing in another order/axis without remapping it.
    if data.get('SliceEncodingDirection', 'k') != 'k':
        raise ValueError('{}: SliceEncodingDirection={} needs explicit axis/order remapping; this pipeline currently requires k.'.format(source, data.get('SliceEncodingDirection')))
    return times


def validate_file(timing_path, image_path):
    lines = Path(timing_path).read_text().splitlines()
    n_slices = nifti_shape(image_path)[2]
    if len(lines) != n_slices:
        raise ValueError('{}: {} timing rows for {} slices in {}. Regenerate timing metadata before resuming headmotion; correction must not continue with this file.'.format(timing_path, len(lines), n_slices, image_path))
    if any(len(line.split()) != 1 for line in lines):
        raise ValueError('{}: timing file must have exactly one value per row'.format(timing_path))
    values = [float(line) for line in lines]
    if any(not math.isfinite(v) or abs(v) >= 1 for v in values):
        raise ValueError('{}: invalid fractional-TR timing shifts'.format(timing_path))
    return n_slices


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--timing', required=True)
    p.add_argument('--image', required=True)
    args = p.parse_args()
    try:
        count = validate_file(args.timing, args.image)
    except (ValueError, OSError) as exc:
        p.exit(2, '[ERROR] slice timing: {}\n'.format(exc))
    print('[slice timing] Validated {} per-slice shifts: {}'.format(count, args.timing), flush=True)
