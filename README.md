# MitoDiffTrack

MitoDiffTrack is a python based workflow for analysing mitochondrial motility in time lapse fluorescence imaging.  
It is implemented as a napari based tool with several dock widgets that guide the user through loading data, segmentation, tracking, motion analysis, and export.

Main panels:
- `load_czi`
- `calibration`
- `seg_params`
- `process_export (seg)`
- `track_mito`
- `batch_process_folder`

Advanced panels:
- `flow_motion_filaments`
- `fast_tracks_skeleton`
- `process_export (CSV)`
- `save_gui_config` / `load_gui_config`

---

## Installation

Create a Python environment and install dependencies:

```
pip install napari[all] magicgui numpy pandas scikit-image scipy opencv-python aicspylibczi qtpy
```

Clone this repository and run:

```
python MitoDiffTrack.py
```

This opens a napari viewer and loads all MitoDiffTrack dock widgets.


---

## 1. Launch the workflow in napari

### Run directly from Python

```
python MitoDiffTrack_V01.py
```

A napari viewer will open with all panels added automatically.

---

## 2. Load the time lapse data and set calibration  
### Panels: `load_czi` + `calibration`

### Load data

1. In the `load_czi` panel, click **Load CZI / VSI**.
2. Select a `.czi` (or supported) raw time lapse file.
3. Choose a Z-projection method: `max`, `median`, or `z0`.
4. Optionally set `Max preview T`.
5. Press **Load CZI / VSI**.

A new `CZI_movie` layer will appear in the viewer.

### Set calibration

1. In the `calibration` panel, enter:
   - **Pixel size (µm/px)**
   - **Frame interval (seconds)**
2. Click **Apply**.

All velocity measurements in exported tables will be in **µm per second**, based on these values.

---

## 3. Preview segmentation and tune parameters  
### Panel: `seg_params`

1. Select the detection **Channel**.
2. Adjust preprocessing (BG sigma, Gaussian sigma, optional CLAHE).
3. Choose a threshold method:
   - `otsu`
   - `mean+std`
   - `percentile`
   - `local`
4. Set object size filters:
   - Minimum area  
   - Maximum area  
   - Optional hole filling, opening, or watershed declumping
5. Press **Preview segmentation at current T**.

A `mito_labels_preview` layer is created.

**Tips**

- If faint mitochondria are missing → reduce threshold or min area.  
- If noise appears → increase threshold or min area.

---

## 4. Run full segmentation of all frames  
### Panel: `process_export (seg)`

1. Set **Channel**, **T start**, **T end**.  
2. Reuse segmentation parameters optimised in the preview.  
3. Press **Run detection on all frames**.

The workflow:

- Processes every frame in the selected time range  
- Adds a `mito_labels` layer  
- Creates a global per-object table (stored internally) containing:
  - Frame index  
  - Centroid  
  - Area  
  - Length, aspect ratio  
  - Intensity features  

---

## 5. LAP tracking and track preview  
### Panel: `track_mito`

Set linking parameters:

- **Max disp (px/frame)**
- **Fast disp (px/frame)** and **Fast penalty factor**
- **Area ratio** and **Intensity ratio** gating
- Display settings: `Color by` (ID or speed), Track tail length

Press **Run LAP tracking from detections**.

The workflow:

- Links objects across consecutive frames (LAP)
- Computes velocities (µm/s) if calibration is set
- Adds a `mito_tracks_preview` layer

Inspect several time points to ensure the tracks follow mitochondria correctly, including fast movers.

---

## 6. Export CSV results and optional snapshot  
### Panel: `process_export (CSV)`

1. Set **Output prefix**.
2. Optional:
   - Save PNG snapshot  
   - Export per-step global LAP table  
   - Export fast-track tables (if run)
3. Press **Export CSV / snapshot**.
4. Choose an output folder.

Exports include:

- Per-object table  
- LAP per-track table (duration, displacement, speed)  
- LAP per-step table  
- Optional fast-track per-track and per-step tables  
- Optional PNG overlay snapshot

All velocities are expressed in **µm per second**.

---

## 7. Motion filaments and fast tracks (optional)  
### Panels: `flow_motion_filaments` and `fast_tracks_skeleton`

### Motion-based filaments

1. Set **T start**, **T end**, **Channel**, and **Downscale**.  
2. Press **Flow motion filaments**.

A label layer `motion_filament_skel` is generated, highlighting main routes of sustained movement.

### Fast tracks in skeleton ROI

1. Choose **ROI dilate** to enlarge the skeleton region.  
2. Apply more aggressive segmentation and linking parameters.  
3. Press **Fast track in skeleton ROI**.

Outputs:

- `mito_fast_tracks_preview` layer  
- Fast-track CSV tables (exported via `process_export (CSV)`)

Useful for detecting very fast or transient motility events.

---

## 8. Batch processing (placeholder)  
### Panel: `batch_process_folder`

This version includes a placeholder batch-processing panel.  
Full automated batch segmentation + tracking will be added in a future release.

For now, process each file using the single-file workflow.

---

## 9. Save and load GUI configurations  
### Panels: `save_gui_config` / `load_gui_config`

To ensure reproducibility:

- **Save GUI config** writes all panel parameters to a JSON file  
- **Load GUI config** restores the settings later

After loading, press **Apply** in the `calibration` panel once to refresh pixel/temporal calibration.

---

## Notes

- Raw movies are not included in the repository.  
- This workflow is intended for mitochondrial motility research.  
- If you use MitoDiffTrack in your research, please cite or acknowledge the contribution of the Frederic Meunier group, Queensland Brain Institute, The University of Queensland.
