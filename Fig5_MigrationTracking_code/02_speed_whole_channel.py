"""
@author: ayesin

Step 2a — Average migration speed in straight channels (no constrictions).

Reads the tracked-feature CSV from 01_segment_and_track.py and computes one
average speed per cell across its whole track.

Two filters are applied:
  - A cell is considered to have stopped once its centroid is unchanged for
    STOP_FRAME_COUNT consecutive frames. Everything from the first such stop
    onward is discarded, so a cell that arrests partway through contributes
    only its motile portion.
  - Cells whose retained track covers less than MIN_TOTAL_DISPLACEMENT_UM of
    total path length are dropped entirely, which removes objects that were
    tracked but never meaningfully migrated.

Speed is total path length divided by total elapsed time, not the mean of the
per-frame instantaneous speeds.
"""

import pandas as pd
import numpy as np


# ============================ CONFIGURATION ==================================

INPUT_CSV  = "path/to/tracked_features.csv"
OUTPUT_CSV = "path/to/whole_channel_speeds.csv"

PIXEL_SIZE_UM = 0.65   # microns per pixel; ignored if COORDS_IN_MICRONS is True
TIME_INTERVAL_MIN = 0.5

# Set True if x_c / y_c in the input CSV are already in microns.
# 01_segment_and_track.py writes microns, so leave this consistent with
# whichever table you are feeding in.
COORDS_IN_MICRONS = False

STOP_THRESHOLD_PX = 0    # displacement at or below this counts as stationary
STOP_FRAME_COUNT = 5     # consecutive stationary frames that end a track
MIN_TOTAL_DISPLACEMENT_UM = 50


# ================================ ANALYSIS ===================================

scale = 1.0 if COORDS_IN_MICRONS else PIXEL_SIZE_UM

df = pd.read_csv(INPUT_CSV)
df.sort_values(by=["cell_label", "time_point"], inplace=True)

average_speeds = {}

for cell_label, group in df.groupby("cell_label"):
    x_coords = group["x_c"].values
    y_coords = group["y_c"].values

    displacements = np.sqrt(np.diff(x_coords) ** 2 + np.diff(y_coords) ** 2)
    displacements_um = displacements * scale

    # Locate the first run of STOP_FRAME_COUNT stationary frames
    stop_indices = [
        i for i in range(len(displacements) - STOP_FRAME_COUNT + 1)
        if np.all(displacements[i:i + STOP_FRAME_COUNT] <= STOP_THRESHOLD_PX)
    ]

    if stop_indices:
        first_stop = min(stop_indices)
        moving_displacements = displacements_um[:first_stop]
        moving_frames = first_stop
    else:
        moving_displacements = displacements_um
        moving_frames = len(displacements_um)

    total_movement = np.sum(moving_displacements)
    total_time = moving_frames * TIME_INTERVAL_MIN

    if total_movement >= MIN_TOTAL_DISPLACEMENT_UM and moving_frames > 0:
        average_speeds[cell_label] = total_movement / total_time

print("Cell Label : Average Speed (um/min)")
for cell_label, speed in average_speeds.items():
    print(f"{cell_label} : {speed:.2f}")

print(f"\nCells retained: {len(average_speeds)} of {df['cell_label'].nunique()}")

output_df = pd.DataFrame(list(average_speeds.items()),
                         columns=["cell_label", "average_speed_um_per_min"])
output_df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved speeds to: {OUTPUT_CSV}")
