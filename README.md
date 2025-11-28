# MitoDiffTrack

MitoDiffTrack is a napari based pipeline for detecting and tracking mitochondrial motility in time lapse fluorescence movies. It combines CellProfiler style object detection, global LAP tracking, and motion based skeleton regions of interest to highlight fast moving mitochondria and export quantitative motility statistics in µm per second. :contentReference[oaicite:0]{index=0}

---

## Key features

- Load multi dimensional Zeiss CZI movies and perform Z projection
- Global LAP based tracking for all detected mitochondria
- Motion filament map from long term optical flow to reveal main “highways” of movement
- Fast track module that performs aggressive detection only inside motion skeleton ROIs
- Export per object, per track, and per step CSV tables including speed in µm per second
- Save and reload full GUI parameter sets as JSON to ensure reproducible analysis

---

## Requirements

Tested with Python 3.x and the following packages:

- `napari`
- `magicgui`
- `numpy`
- `pandas`
- `scikit-image`
- `scipy`
- `opencv-python` (cv2)
- `aicspylibczi`
- `qtpy`

Recommended workflow

Inside the napari window, follow the numbered docks from top to bottom:

📂 Load CZI

Load a .czi time lapse movie

Choose Z projection mode and time range preview if needed

⏱ Pixel / Time

Set pixel size in µm per pixel

Set frame interval in seconds

These values are required to compute speed_um_per_s for tracks

🔍 Segmentation preview (optional)

Preview object detection on the current time frame

Adjust threshold mode, smoothing, area filters and other parameters

① Detection (all frames)

Run intensity based object detection across a chosen time range

Generates a global per object table (centroid, area, shape, intensity)

② LAP tracking (global)

Perform global LAP based tracking on all detected mitochondria

Builds a preview track layer in napari

Computes displacement and speed in µm per second when pixel size and frame interval are set

③ Flow motion filaments

Build a motion based filament skeleton from long term optical flow patterns

The output layer motion_filament_skel highlights main paths of sustained mitochondrial movement

③b Fast tracks (skeleton ROI)

Dilate the motion skeleton to create a motion ROI

Run aggressive detection only inside this ROI

Apply LAP tracking with relaxed displacement limits to capture very fast movers

Outputs a separate fast tracks preview layer and per step table with speed_um_per_s

④ Export CSV / PNG

Export:

Global per object detection table

Global LAP per track summary

Optional global LAP per step table

Optional fast track per track and per step tables

Optional snapshot PNG of tracks in napari

All exports require valid pixel size and frame interval to ensure consistent speed units

💾 Save GUI config / 📂 Load GUI config

Save current GUI parameter settings into a JSON config file

Reload a config later to reproduce exactly the same analysis settings

Exported files

Depending on the options you select in ④ Export CSV / PNG, the pipeline can produce:

*_objects_YYYYMMDD_HHMMSS.csv
Global per object detection table for all frames

*_LAP_tracks_YYYYMMDD_HHMMSS.csv
Per track summary for global LAP tracking, including duration and speed statistics

*_LAP_steps_YYYYMMDD_HHMMSS.csv (optional)
Per step table for global LAP tracking

*_FAST_tracks_YYYYMMDD_HHMMSS.csv and *_FAST_steps_YYYYMMDD_HHMMSS.csv (optional)
Per track and per step tables for fast tracks inside the skeleton ROI

*_tracks_snapshot_YYYYMMDD_HHMMSS.png (optional)
PNG snapshot of the current napari view with tracks overlaid

All speed related columns are reported in µm per second when acquisition parameters are set correctly

You can create a virtual environment and install dependencies, for example:

```bash
pip install napari[all] magicgui numpy pandas scikit-image scipy opencv-python aicspylibczi qtpy
