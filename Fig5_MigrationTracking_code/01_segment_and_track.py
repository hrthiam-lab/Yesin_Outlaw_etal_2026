"""
@author: Minwoo Kang and Allen Yesin

Code adapted from:
M. Kang, et al., Classification of differentially activated groups of fibroblasts using morphodynamic and motile features. 
APL Bioeng. 9, 026116 (2025).

Step 1 — Nuclear segmentation and frame-to-frame tracking.

Takes a time-lapse sequence of single-channel fluorescence images (one TIFF per
frame, nuclei labelled with a DNA dye) and produces a per-cell, per-frame table
of morphological features with a persistent `cell_label` across frames.

Pipeline:
  1. Segment nuclei with a watershed on each frame and save binary masks.
  2. Extract region properties from the masks, using the fluorescence image
     for intensity features (background-subtracted).
  3. Link objects between consecutive frames by solving a linear assignment
     problem on centroid displacement penalised by area change.
  4. Optionally write labelled bounding-box overlays for visual QC.
  5. Add a `time_point` column and save the table as CSV.

The CSV produced here is the input to 02_speed_whole_channel.py and
03_speed_by_region.py.

"""

import os
import itertools
from collections import Counter
from operator import itemgetter
from itertools import groupby

import cv2
import numpy as np
import pandas as pd
from tqdm import tqdm
from skimage import measure, color
from skimage.segmentation import clear_border
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from natsort import natsort, natsorted, natsort_keygen
import mahotas


# ============================ CONFIGURATION ==================================

# --- Paths (edit these) ---
NUC_PATH        = "path/to/fluorescence_frames/"   # input TIFF sequence
NUC_MASK_DIR    = "path/to/masks/"                 # where binary masks are written
BBOX_DIR        = "path/to/bbox_overlays/"         # QC overlays (optional)
OUTPUT_CSV      = "path/to/tracked_features.csv"   # final table

# --- Acquisition parameters ---
PIXEL_SIZE_UM   = 0.65    # microns per pixel
TIME_INTERVAL   = 0.5     # minutes per frame

# --- Segmentation ---
AREA_THRESH_PX  = 50      # objects below this pixel area are discarded as debris
MIN_MARKER_PX   = 10      # watershed seeds below this pixel area are dropped
FG_DIST_FRAC    = 0.2     # distance-transform fraction defining sure foreground

# --- Tracking (Jaqaman-style linear assignment) ---
DIST_CUTOFF_PX  = 200     # maximum centroid displacement between frames, pixels
SIZE_CUTOFF     = 3.5     # maximum fold change in object area between frames
INT_EFF         = 1.4     # exponent applied to the area-change penalty

# --- Bookkeeping ---
SAMPLE_CONDITION = "condition_name"   # free-text label written to every row
WRITE_BBOX_QC    = True

# NOTE ON UNITS: this script writes x_c / y_c in MICRONS (pixels x PIXEL_SIZE_UM).
# The downstream speed scripts expect to know which unit they are given — set
# COORDS_IN_MICRONS accordingly there.


#%% ===================== 1. Segment nuclei, save masks =======================

nuc_image = natsort.natsorted(os.listdir(NUC_PATH))
os.makedirs(NUC_MASK_DIR, exist_ok=True)

for num, image in tqdm(enumerate(nuc_image), total=len(nuc_image)):
    if not image.lower().endswith(".tif"):
        continue

    img = cv2.imread(os.path.join(NUC_PATH, image), cv2.IMREAD_UNCHANGED)
    img_grey = cv2.convertScaleAbs(img, alpha=(255.0 / 65535.0))
    img_color = cv2.cvtColor(img_grey, cv2.COLOR_GRAY2BGR)

    # Threshold at the modal pixel value, i.e. the background peak
    blurred_img = cv2.GaussianBlur(img_grey, (5, 5), 0)
    binary_thresh = Counter(blurred_img.flatten()).most_common(1)[0][0]

    _, thresh = cv2.threshold(blurred_img, binary_thresh, 255, cv2.THRESH_BINARY)
    thresh = np.uint8(thresh)

    kernel = np.ones((3, 3), np.uint8)
    closing = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel, iterations=1)
    sure_bg = cv2.dilate(closing, kernel, iterations=3)

    dist_transform = cv2.distanceTransform(closing, cv2.DIST_L2, 5)
    _, sure_fg = cv2.threshold(dist_transform, FG_DIST_FRAC * dist_transform.max(), 255, 0)
    sure_fg = np.uint8(sure_fg)
    unknown = cv2.subtract(sure_bg, sure_fg)

    _, markers = cv2.connectedComponents(sure_fg)
    stats = cv2.connectedComponentsWithStats(sure_fg)[2]
    for i, stat in enumerate(stats):
        if stat[cv2.CC_STAT_AREA] < MIN_MARKER_PX:
            markers[markers == i] = 0

    markers = markers + 10
    markers[unknown == 255] = 0
    markers = cv2.watershed(img_color, markers)

    binary_mask = np.zeros(img_grey.shape, dtype=np.uint16)
    binary_mask[markers > 10] = 65535
    cv2.imwrite(os.path.join(NUC_MASK_DIR, f"nuc_mask_{num}.tif"), binary_mask)


#%% ===================== 2. Extract per-object features ======================

recon_image = natsort.natsorted(os.listdir(NUC_MASK_DIR))
inten_image = natsort.natsorted(os.listdir(NUC_PATH))
recon_image = [f for f in recon_image if f != ".DS_Store"]
inten_image = [f for f in inten_image if f != ".DS_Store"]

region_props_list = []
file_names = []
cell_label = []
zm_list = []
lst_contour_bbx = []

for num, (image, image_intensity) in tqdm(enumerate(zip(recon_image, inten_image)),
                                          total=len(recon_image)):
    if not image.lower().endswith(".tif"):
        continue

    img = cv2.imread(os.path.join(NUC_MASK_DIR, image))
    dim = (int(img.shape[1]), int(img.shape[0]))
    img_grey = cv2.resize(img[:, :, 0], dim)
    img_grey = clear_border(img_grey)

    img_intensity = cv2.imread(os.path.join(NUC_PATH, image_intensity), -1)

    # Background subtraction: mean intensity outside a dilated cell mask
    kernel = np.ones((3, 3), np.uint8)
    sure_bg_mask = cv2.dilate(img_grey, kernel, iterations=3)
    inverse_sure_bg_mask = np.uint16(1 - (sure_bg_mask / 255))
    background_intensity = np.uint16(img_intensity * inverse_sure_bg_mask)
    background_mean = np.full(background_intensity.shape, np.mean(background_intensity))

    img_intensity_corrected = img_intensity - background_mean
    img_intensity_corrected[img_intensity_corrected < 0] = 0

    # Relabel the mask, dropping objects below the area threshold
    img_8bit = cv2.imread(os.path.join(NUC_MASK_DIR, image))[:, :, 0]
    ret3, markers, box_stats, _ = cv2.connectedComponentsWithStats(img_8bit)

    img_without_smallarea = np.zeros(markers.shape, np.uint8)
    for k in range(1, ret3):
        if box_stats[k][4] > AREA_THRESH_PX:
            img_without_smallarea[markers == k] = 1

    nlabel, labels, stats, centroids = cv2.connectedComponentsWithStats(img_without_smallarea)
    img2 = color.label2rgb(labels, bg_label=0)

    props = measure.regionprops_table(
        labels, intensity_image=img_intensity_corrected,
        properties=['area', 'extent', 'centroid', 'perimeter',
                    'minor_axis_length', 'major_axis_length',
                    'equivalent_diameter', 'solidity', 'bbox',
                    'intensity_mean', 'intensity_max', 'intensity_min', 'intensity_std'])

    # Zernike moments per object, computed on the cropped binary mask
    contours, _ = cv2.findContours(img_without_smallarea.copy(),
                                   cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    for c in contours:
        empty_plane = np.zeros(img_without_smallarea.shape, dtype="uint8")
        cv2.drawContours(empty_plane, [c], -1, 255, -1)

        (xc, yc, wc, hc) = cv2.boundingRect(c)
        contour_bbx = np.array([xc, yc, xc + wc, yc + hc]).reshape(1, 4)

        cropped_cell_binary = empty_plane[yc:yc + hc, xc:xc + wc]
        _, radius = cv2.minEnclosingCircle(c)
        zernike_moments = mahotas.features.zernike_moments(cropped_cell_binary, radius, degree=9)

        df_zernike_moments = pd.DataFrame(zernike_moments).T
        for i in range(zernike_moments.shape[0]):
            df_zernike_moments = df_zernike_moments.rename(columns={i: f'zernike_moments_{i}'})

        zm_list.append(df_zernike_moments.to_numpy())
        lst_contour_bbx.append(contour_bbx)
        df_zm = pd.DataFrame.from_records(itertools.chain.from_iterable(zm_list))
        df_cont_bbx = pd.DataFrame.from_records(itertools.chain.from_iterable(lst_contour_bbx))

    region_props_list.append(pd.DataFrame(props).to_numpy())
    df_props = pd.DataFrame.from_records(itertools.chain.from_iterable(region_props_list))
    df_props.columns = ['area', 'extent', 'y_c', 'x_c', 'perimeter',
                        'minor_axis_length', 'major_axis_length',
                        'equivalent_diameter', 'solidity',
                        'bbox_y1', 'bbox_x1', 'bbox_y2', 'bbox_x2',
                        'intensity_mean', 'intensity_max', 'intensity_min', 'intensity_std']

    for i in range(1, nlabel):
        file_names.append(image)
        df0 = pd.DataFrame(file_names, columns=['file_name'])
        df = pd.concat([df0, df_props], axis=1)


#%% ===================== 3. Derived features and scaling =====================

df = df[['file_name', 'area', 'x_c', 'y_c', 'extent', 'perimeter',
         'minor_axis_length', 'major_axis_length', 'equivalent_diameter',
         'solidity', 'bbox_x1', 'bbox_y1', 'bbox_x2', 'bbox_y2',
         'intensity_mean', 'intensity_max', 'intensity_min', 'intensity_std']]

df['aspect_ratio'] = df['major_axis_length'] / df['minor_axis_length']
df['circularity'] = 4 * np.pi * df['area'] / (df['perimeter'] ** 2)

# Convert pixel-based features to microns
df['area'] = df['area'] * PIXEL_SIZE_UM * PIXEL_SIZE_UM
for col in ['x_c', 'y_c', 'perimeter', 'minor_axis_length', 'major_axis_length',
            'equivalent_diameter', 'bbox_x1', 'bbox_y1', 'bbox_x2', 'bbox_y2']:
    df[col] = df[col] * PIXEL_SIZE_UM

df['compactness'] = df['perimeter'] ** 2 / df['area']

df_cont_bbx.columns = ['bbox_x1', 'bbox_y1', 'bbox_x2', 'bbox_y2']
df_zm = pd.DataFrame(data=df_zm.values, columns=df_zernike_moments.columns)
df_temp = pd.concat([df_cont_bbx, df_zm], axis=1)
for col in ['bbox_x1', 'bbox_y1', 'bbox_x2', 'bbox_y2']:
    df_temp[col] = df_temp[col] * PIXEL_SIZE_UM

df = pd.merge(df, df_temp)
df.sort_values(by="file_name", key=natsort_keygen(), inplace=True, ignore_index=True)

fname = natsorted(list(df.loc[:, 'file_name'].unique()))

# Drop duplicate objects within a frame
df_ed_list = []
for tm in range(len(fname)):
    df_ed0 = df[df.loc[:, 'file_name'] == fname[tm]].drop_duplicates(
        subset=['area', 'x_c'], keep='first', ignore_index=True)
    df_ed_list.append(df_ed0.to_numpy())
    df_ed = pd.DataFrame.from_records(itertools.chain.from_iterable(df_ed_list))
    df_ed.columns = df.columns
df = df_ed

# Provisional per-frame labels, overwritten by the linking step below
for tmstp in range(len(fname)):
    lbl = len(df[df.loc[:, 'file_name'] == fname[tmstp]].values)
    for j in range(lbl):
        cell_label.append(j)
        df1 = pd.DataFrame(cell_label, columns=['cell_label'])

df = pd.concat([df, df1], axis=1)


#%% ===================== 4. Link objects across frames =======================
# Cost of linking object i (frame t) to object j (frame t+1) is the squared
# centroid displacement multiplied by an area-change penalty (r + 1/r)^INT_EFF.
# Links exceeding DIST_CUTOFF_PX or SIZE_CUTOFF are pushed above the cost of
# leaving an object unmatched, so they are never accepted.

dist_cutoff = DIST_CUTOFF_PX * PIXEL_SIZE_UM

for timestp in range(len(fname) - 1):
    try:
        ct00 = df[df.loc[:, 'file_name'] == fname[timestp]]
        ct01 = df[df.loc[:, 'file_name'] == fname[timestp + 1]]
        idx0, idx1 = list(ct00.index), list(ct01.index)

        ct0, ct1 = ct00.values, ct01.values
        cent0 = ct0[:, 2:4].astype(float)
        cent1 = ct1[:, 2:4].astype(float)
        area0 = ct0[:, 1].reshape((-1, 1))
        area1 = ct1[:, 1].reshape((-1, 1))

        # Area-change penalty
        int_dist_mat = area1.T / area0
        int_dist_mat = int_dist_mat + 1 / int_dist_mat
        int_dist_mat[int_dist_mat >= (SIZE_CUTOFF + 1 / SIZE_CUTOFF)] = 20.0
        int_dist_mat = int_dist_mat ** INT_EFF
        int_dist_baseline = np.percentile(int_dist_mat, 10)

        # Square cost matrix with no-match blocks on the diagonal
        cost_mat = np.ones((len(cent0) + len(cent1), len(cent0) + len(cent1))) \
            * (dist_cutoff ** 2 * 10) * int_dist_baseline
        dist_mat = cdist(cent0, cent1) ** 2
        dist_mat[dist_mat >= (dist_cutoff ** 2)] = dist_cutoff ** 2 * 10
        cost_mat[:len(cent0), :len(cent1)] = dist_mat * int_dist_mat

        for i in range(len(cent0)):
            cost_mat[i, i + len(cent1)] = 1.05 * (dist_cutoff ** 2) * int_dist_baseline
        for j in range(len(cent1)):
            cost_mat[len(cent0) + j, j] = 1.05 * (dist_cutoff ** 2) * int_dist_baseline
        cost_mat[len(cent0):, len(cent1):] = np.transpose(dist_mat)

        links = linear_sum_assignment(cost_mat)
        pairs = [p for p in zip(*links) if p[0] < len(cent0) and p[1] < len(cent1)]
        pairs.sort(key=lambda x: x[1])

        # Carry the frame-t label onto the matched object in frame t+1
        matching_idx = []
        for m in range(len(pairs)):
            target = ct01[ct01.loc[:, 'cell_label'] == ct1[pairs[m][1], -1]].index.astype(int)[0]
            matching_idx.append(target)
            df.loc[target, 'cell_label'] = ct0[pairs[m][0], -1]

        # Unmatched objects in frame t+1 become new tracks
        no_matching_idx = sorted(set(idx1) - set(matching_idx))
        lst_new_lbl = [max(df.loc[:, 'cell_label']) + 1]
        for i in range(1, len(no_matching_idx)):
            lst_new_lbl.append(lst_new_lbl[0] + i)
        for i, j in enumerate(no_matching_idx):
            df.loc[j, 'cell_label'] = lst_new_lbl[i]

        # Close gaps so labels stay consecutive
        all_lbl_list = list(df.loc[:, 'cell_label'].unique())
        group_consec_lbl = []
        for k, g in groupby(enumerate(all_lbl_list), lambda x: x[0] - x[1]):
            group = list(map(itemgetter(1), g))
            group_consec_lbl.append((group[0], group[-1]))

        common_diff = 0 if len(group_consec_lbl) == 1 else \
            group_consec_lbl[1][0] - (group_consec_lbl[0][1] + 1)

        for nl in range(len(df)):
            if df.loc[nl, 'cell_label'] > group_consec_lbl[0][1]:
                df.loc[nl, 'cell_label'] -= common_diff

    except Exception as e:
        print(f"Error during frame matching at timestamp {timestp}: {e}")

df['image_stack'] = NUC_MASK_DIR
df['condition'] = SAMPLE_CONDITION


#%% ===================== 5. Bounding-box overlays (QC) =======================
# Writes one PNG per frame with each tracked object boxed and labelled. Scroll
# through these to confirm labels stay on the same cell across frames.

if WRITE_BBOX_QC:
    os.makedirs(BBOX_DIR, exist_ok=True)

    for num, image in tqdm(enumerate(natsort.natsorted(os.listdir(NUC_MASK_DIR))),
                           total=len(recon_image)):
        if not image.lower().endswith(".tif"):
            continue

        img = cv2.imread(os.path.join(NUC_MASK_DIR, image))
        if img is None:
            print(f"Warning: could not read {image}, skipping.")
            continue

        img_b = cv2.resize(img.copy(), dim)
        sel = df.loc[:, 'file_name'] == fname[num]
        x1 = df.loc[sel, 'bbox_x1'].tolist()
        y1 = df.loc[sel, 'bbox_y1'].tolist()
        x2 = df.loc[sel, 'bbox_x2'].tolist()
        y2 = df.loc[sel, 'bbox_y2'].tolist()
        lbls = df.loc[sel, 'cell_label'].tolist()

        for k in range(len(lbls)):
            xs, ys = int(x1[k] / PIXEL_SIZE_UM), int(y1[k] / PIXEL_SIZE_UM)
            xe, ye = int(x2[k] / PIXEL_SIZE_UM), int(y2[k] / PIXEL_SIZE_UM)
            cv2.rectangle(img_b, (xs, ys), (xe, ye), (0, 255, 255), 2)
            cv2.putText(img_b, f'{lbls[k]}', (xs, ys - 10),
                        cv2.FONT_HERSHEY_SIMPLEX, 1.8, (36, 255, 12), 2)

        cv2.imwrite(os.path.join(BBOX_DIR, f"cyto_bbox_{num}.png"), img_b)


#%% ===================== 6. Add time and save ================================

df.sort_values(by=['cell_label', 'file_name'], inplace=True, ignore_index=True)

df['time_point'] = 0
for fn in range(len(fname)):
    for idx in df[df.loc[:, 'file_name'] == fname[fn]].index:
        df.loc[idx, 'time_point'] = fn * TIME_INTERVAL

df.sort_values(by=['cell_label', 'time_point'], key=natsort_keygen(),
               inplace=True, ignore_index=True)

df.to_csv(OUTPUT_CSV, index=True)
print(f"Saved tracked features to: {OUTPUT_CSV}")
