#!/usr/bin/env python3
"""Append a zoomed coronal subcortical-label montage to a RidgePFM frame."""
from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to
from PIL import Image, ImageDraw, ImageFont


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frame', required=True)
    parser.add_argument('--labels-volume', required=True)
    parser.add_argument('--coverage-volume', required=True)
    parser.add_argument('--label-cifti', required=True)
    parser.add_argument('--anatomical', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--axial-mm', nargs='+', type=float, default=[-34, -14, 4, 14])
    args = parser.parse_args()

    label_img = nib.load(args.labels_volume)
    labels = np.asarray(label_img.dataobj)
    if labels.ndim == 4:
        labels = labels[..., 0]
    coverage = np.asarray(nib.load(args.coverage_volume).dataobj)
    if coverage.ndim == 4:
        coverage = coverage[..., 0]
    labels = np.where(coverage == 2, labels, 0).astype(int)
    anatomy = resample_from_to(
        nib.load(args.anatomical), (label_img.shape[:3], label_img.affine), order=1
    ).get_fdata(dtype=np.float32)
    positive = anatomy[anatomy > 0]
    lo, hi = np.percentile(positive, [2, 99])
    anatomy = np.clip((anatomy - lo) / max(hi - lo, 1e-6), 0, 1)

    table = nib.load(args.label_cifti).header.get_axis(0).label[0]
    rgba = np.zeros((max(table) + 1, 4), np.uint8)
    for key, (_, color) in table.items():
        rgba[key] = np.rint(np.asarray(color) * 255).astype(np.uint8)

    base = Image.open(args.frame).convert('RGB')
    insert_y = 805
    montage_h = 270
    canvas = Image.new('RGB', (base.width, base.height + montage_h), 'white')
    canvas.paste(base.crop((0, 0, base.width, insert_y)), (0, 0))
    canvas.paste(base.crop((0, insert_y, base.width, base.height)), (0, insert_y + montage_h))
    draw = ImageDraw.Draw(canvas)
    font_path = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
    draw.text((28, insert_y + 7), 'Subcortical assignments · axial slices',
              font=ImageFont.truetype(font_path, 18), fill='#334b61')

    xlo, xhi, ylo, yhi = 8, 83, 18, 93
    gap, panel_w, panel_h = 12, 276, 214
    start_x, top = 24, insert_y + 38
    inv = np.linalg.inv(label_img.affine)
    for panel, millimeters in enumerate(args.axial_mm):
        voxel_z = int(round(nib.affines.apply_affine(inv, [0, 0, millimeters])[2]))
        gray = anatomy[xlo:xhi, ylo:yhi, voxel_z].T[::-1, ::-1]
        lab = labels[xlo:xhi, ylo:yhi, voxel_z].T[::-1, ::-1]
        background = Image.fromarray(np.rint(gray * 255).astype(np.uint8)).resize(
            (panel_w, panel_h), Image.Resampling.BILINEAR
        )
        image = np.repeat(np.asarray(background)[..., None], 3, axis=2)
        lab = np.asarray(Image.fromarray(lab.astype(np.int16)).resize((panel_w, panel_h), Image.Resampling.NEAREST))
        selected = lab > 0
        colors = rgba[np.clip(lab, 0, len(rgba) - 1), :3]
        image[selected] = np.rint(0.80 * colors[selected] + 0.20 * image[selected]).astype(np.uint8)
        # Compute borders after enlargement so boundaries remain one pixel wide.
        edge = selected & (
            (np.roll(lab, 1, 0) != lab) | (np.roll(lab, -1, 0) != lab)
            | (np.roll(lab, 1, 1) != lab) | (np.roll(lab, -1, 1) != lab)
        )
        image[edge] = 0
        panel_image = Image.fromarray(image)
        x = start_x + panel * (panel_w + gap)
        canvas.paste(panel_image, (x, top))
        draw.rectangle((x, top, x + panel_w - 1, top + panel_h - 1), outline='#2f3944', width=1)
        draw.text((x + 8, top + 7), f'z = {millimeters:g} mm',
                  font=ImageFont.truetype(font_path, 15), fill='white', stroke_width=2, stroke_fill='black')

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output)


if __name__ == '__main__':
    main()
