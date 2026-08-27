# Fig4_NuclearEnvelopeIF_code

Python code for quantifying nuclear envelope protein intensity and nuclear
shape in fixed, immunostained single cells.

## Requirements

```
pip install numpy scikit-image scipy tifffile openpyxl
```

## Input

One cropped single-cell TIFF z-stack per file, multi-channel, in one folder.

## Running

```
python nuclear_IF_analysis.py
```

Set in the dialog (1-based channel indices):

- Channel for lamins
- Channel for LBR (0 if this dataset has no LBR stain)
- Channel for the nucleus mask (DAPI)
- Thresholding method for the DAPI mask
- Pixel size (µm/px)
- LBR shape percentile cutoff
- LBR shape edge smoothing (sigma)

## What it does

Channels are collapsed with an average-intensity projection. A nuclear mask
is built from the DAPI channel (threshold, fill holes, erode); background
(mean outside the mask) is subtracted from each channel, and the DAPI mask
is used to measure total and mean lamin, LBR, and DAPI intensity.

Nuclear shape (area, perimeter, circularity, roundness, solidity, aspect
ratio, compactness) is measured separately from the LBR channel, at its
best-focus Z-plane, thresholded by percentile within the DAPI mask.

If a dataset has no LBR channel (set to 0), protein/DNA intensities are
still quantified from the DAPI mask; LBR intensity and all shape fields are
left blank for those cells.

## Output

Per-cell TIFFs (projection, DAPI masks, masked intensity images, LBR shape
mask) plus a single `Analysis_Results.xlsx` with a Parameters sheet and an
All Data sheet (one row per cell).
