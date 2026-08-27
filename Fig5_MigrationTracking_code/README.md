# Fig5_MigrationTracking_code

Measures nuclear migration speed in PDMS microfluidic channels.

```
frames/ ──► 01_segment_and_track.py ──► tracked_features.csv
                                               │
                              ┌────────────────┴────────────────┐
                              ▼                                 ▼
                02_speed_whole_channel.py           03_speed_by_region.py
                  (straight channels)                 (constrictions)
```

Scripts 2 and 3 are alternatives, not sequential steps.

## Input

One TIFF per timepoint, single channel (DNA), in one folder per field of view.
Filenames must sort into acquisition order; natural sorting is used, so
zero-padding is not required.

## 1. Segment and track

Edit the `CONFIGURATION` block at the top of `01_segment_and_track.py` — paths,
`PIXEL_SIZE_UM`, `TIME_INTERVAL`, and the segmentation and linking cutoffs.
The script runs as six `#%%` cells.

Check the masks written to `NUC_MASK_DIR` and the labelled overlays written to
`BBOX_DIR` before continuing. A cell should keep the same label for as long as
it is in frame.

Output: one row per object per frame, with `cell_label`, `time_point`,
`x_c`, `y_c`, and nuclear shape and intensity features.

## 2a. Whole-channel speed

`02_speed_whole_channel.py` reports one speed per cell: total path length over
total elapsed time. A cell is treated as arrested once its centroid is
unchanged for `STOP_FRAME_COUNT` frames, and cells covering less than
`MIN_TOTAL_DISPLACEMENT_UM` are dropped.

## 2b. Region-specific speed

`03_speed_by_region.py` prompts for the migration axis (`x` or `y`), the
direction of travel along it, and the two coordinates bounding the
constriction. Direction sets which side of each boundary is Before and which is
After; all four orientations are supported.

Read the boundary coordinates off the Fiji status bar on one frame, and enter
them in the same unit as `x_c` / `y_c` in your CSV.

Speeds are reported only for cells observed in all three regions. Cells that
stall earlier appear in the categorisation summary printed to the console,
which is the basis for passage-rate analysis.

Output: `<input>_processedvelocity.csv`, one row per cell with before, during
and after speeds in µm/min.

## Troubleshooting

| Symptom | Fix |
| ------- | --- |
| Nuclei merge into one object | Lower `FG_DIST_FRAC` |
| One nucleus splits into fragments | Raise `FG_DIST_FRAC` |
| Labels swap between neighbours | Lower `DIST_CUTOFF_PX` |
| Labels change on the same cell | Raise `DIST_CUTOFF_PX` or `SIZE_CUTOFF` |
| Speeds off by a constant factor | Check `PIXEL_SIZE_UM` and `COORDS_IN_MICRONS` |
| Speeds off by exactly 2× | `TIME_INTERVAL` does not match acquisition |
| No cells traverse all three regions | Boundaries in the wrong unit, or direction set backwards |
