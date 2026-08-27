"""
@author: ayesin

Step 2b — Migration speed before, during and after a constriction.

Reads the tracked-feature CSV from 01_segment_and_track.py, splits each track
into Before / During / After segments using two user-supplied boundary
coordinates, and reports an average speed per segment.

The script prompts for:
  - the axis along which cells migrate (x or y)
  - the direction of travel along that axis
  - the two boundary coordinates that delimit the constriction

Direction matters because it determines which side of a boundary counts as
Before and which as After. All four combinations are supported: left-to-right
and right-to-left along x, top-to-bottom and bottom-to-top along y.

Filtering behaviour:
  - Speeds are reported only for cells observed in all three regions. Cells
    that stall in Before, or in Before and During, are counted in the
    categorisation summary but excluded from the speed table.
  - A track is truncated if the cell re-enters the constriction after having
    reached the After region.
  - Within the After region, once the centroid is unchanged for
    STATIONARY_AFTER_FRAMES consecutive frames the cell is treated as arrested
    and subsequent frames are not accumulated.
"""

import pandas as pd
import numpy as np
from pathlib import Path


# ============================ CONFIGURATION ==================================

INPUT_CSV = "path/to/tracked_features.csv"

PIXEL_SIZE_UM = 0.65   # microns per pixel; ignored if COORDS_IN_MICRONS is True
TIME_INTERVAL_MIN = 0.5

# Set True if x_c / y_c in the input CSV are already in microns.
COORDS_IN_MICRONS = False

STATIONARY_AFTER_FRAMES = 5

# Boundary values are entered interactively, in the same unit as x_c / y_c
# in the input CSV.


# ============================== PROMPTS ======================================

def prompt_axis():
    while True:
        ax = input("Segment along which axis? Type 'x' or 'y': ").strip().lower()
        if ax in {"x", "y"}:
            return ax
        print("Please enter 'x' or 'y'.")


def prompt_direction(axis):
    if axis == "x":
        options = {"1": "left-to-right", "2": "right-to-left"}
        prompt = "Movement direction (x): [1] left-to-right, [2] right-to-left: "
    else:
        options = {"1": "top-to-bottom", "2": "bottom-to-top"}
        prompt = "Movement direction (y): [1] top-to-bottom, [2] bottom-to-top: "
    while True:
        choice = input(prompt).strip()
        if choice in options:
            return options[choice]
        print("Please choose 1 or 2.")


def prompt_float(msg):
    while True:
        try:
            return float(input(f"{msg}: ").strip())
        except ValueError:
            print("Please enter a number.")


# ================================ ANALYSIS ===================================

def main():
    csv_path = Path(INPUT_CSV)
    scale = 1.0 if COORDS_IN_MICRONS else PIXEL_SIZE_UM

    axis = prompt_axis()
    direction = prompt_direction(axis)

    print("\nEnter the region boundaries along the chosen axis, in the same "
          "unit as the x_c / y_c columns of your CSV.")
    before_val = prompt_float("Before-constriction boundary value")
    after_val = prompt_float("After-constriction boundary value")

    seg_key = "x_c" if axis == "x" else "y_c"

    # Which side of each boundary is Before and which is After depends on
    # the direction of travel.
    increasing = direction in {"left-to-right", "top-to-bottom"}
    if increasing:
        before_pred = lambda v: v < before_val
        after_pred = lambda v: v > after_val
        before_desc = f"{axis} < {before_val}"
        after_desc = f"{axis} > {after_val}"
    else:
        before_pred = lambda v: v > before_val
        after_pred = lambda v: v < after_val
        before_desc = f"{axis} > {before_val}"
        after_desc = f"{axis} < {after_val}"

    def region_of(val):
        if before_pred(val):
            return "before"
        if after_pred(val):
            return "after"
        return "during"

    print("\nSegmentation summary:")
    print(f"  Axis: {axis}")
    print(f"  Direction: {direction}")
    print(f"  BEFORE region: {before_desc}")
    print(f"  AFTER region:  {after_desc}")
    print("  DURING region: between the two boundaries.")
    print(f"  Scale: {scale} um per coordinate unit, dt = {TIME_INTERVAL_MIN} min/frame, "
          f"After stationary cutoff = {STATIONARY_AFTER_FRAMES} frames.\n")

    df = pd.read_csv(csv_path)
    required_cols = {"cell_label", "time_point", "x_c", "y_c"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df.sort_values(by=["cell_label", "time_point"], inplace=True)

    average_speeds = {}
    counts = {"before_only": 0, "before_during": 0, "all_three": 0, "total_cells": 0}
    grouped = df.groupby("cell_label")
    counts["total_cells"] = len(grouped)

    for cell_label, g in grouped:
        x = g["x_c"].to_numpy()
        y = g["y_c"].to_numpy()
        s = g[seg_key].to_numpy()

        in_before = in_during = in_after = False
        before_disps, during_disps, after_disps = [], [], []
        before_steps = during_steps = after_steps = 0

        after_stationary_frames = 0
        last_after_pos = None
        after_tracking_active = True

        for i in range(1, len(s)):
            disp = np.hypot(x[i] - x[i - 1], y[i] - y[i - 1]) * scale

            prev_region = region_of(s[i - 1])
            curr_region = region_of(s[i])

            # Stop tracking if the cell falls back into the constriction
            if prev_region == "after" and curr_region == "during":
                break

            if curr_region == "after":
                curr_pos = (x[i], y[i])
                if last_after_pos is not None and curr_pos == last_after_pos:
                    after_stationary_frames += 1
                    if after_stationary_frames >= STATIONARY_AFTER_FRAMES:
                        after_tracking_active = False
                else:
                    after_stationary_frames = 0
                    last_after_pos = curr_pos

            if curr_region == "after" and not after_tracking_active:
                continue

            if curr_region == "before":
                in_before = True
            elif curr_region == "during":
                in_during = True
            elif curr_region == "after":
                in_after = True

            # Assign each step to a region based on the transition it represents
            if prev_region == "before" and curr_region in {"before", "during"}:
                before_disps.append(disp); before_steps += 1
            elif prev_region == "during" and curr_region in {"during", "after"}:
                during_disps.append(disp); during_steps += 1
            elif prev_region == "after" and curr_region == "after":
                after_disps.append(disp); after_steps += 1

        if in_before and not in_during and not in_after:
            counts["before_only"] += 1
        elif in_before and in_during and not in_after:
            counts["before_during"] += 1
        elif in_before and in_during and in_after:
            counts["all_three"] += 1

        if in_before and in_during and in_after:
            def avg_speed(disps, steps):
                if not disps or steps == 0:
                    return np.nan
                total_t = steps * TIME_INTERVAL_MIN
                return float(np.sum(disps)) / total_t if total_t > 0 else np.nan

            average_speeds[cell_label] = {
                "before_constriction_speed_um_per_min": avg_speed(before_disps, before_steps),
                "during_constriction_speed_um_per_min": avg_speed(during_disps, during_steps),
                "after_constriction_speed_um_per_min": avg_speed(after_disps, after_steps),
            }

    total = counts["total_cells"]
    pct = lambda n: (n / total * 100.0) if total else 0.0

    print("\nCell categorisation:")
    print(f"  Total cells: {total}")
    print(f"  Before only: {counts['before_only']} ({pct(counts['before_only']):.2f}%)")
    print(f"  Before and during: {counts['before_during']} ({pct(counts['before_during']):.2f}%)")
    print(f"  All three regions: {counts['all_three']} ({pct(counts['all_three']):.2f}%)")

    if average_speeds:
        print("\nCell Label : Before : During : After  (um/min)")
        for cell_label, sp in average_speeds.items():
            print(f"{cell_label} : "
                  f"{sp['before_constriction_speed_um_per_min']:.2f} : "
                  f"{sp['during_constriction_speed_um_per_min']:.2f} : "
                  f"{sp['after_constriction_speed_um_per_min']:.2f}")
    else:
        print("\nNo cells traversed all three regions.")

    out_df = (pd.DataFrame.from_dict(average_speeds, orient="index")
              .reset_index().rename(columns={"index": "cell_label"}))
    out_path = csv_path.with_name(csv_path.stem + "_processedvelocity.csv")
    out_df.to_csv(out_path, index=False)
    print(f"\nSaved speeds to: {out_path}")


if __name__ == "__main__":
    main()
