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

Then clone this repository and run:

```
python MitoDiffTrack.py
```

This opens a napari viewer and loads all MitoDiffTrack dock widgets.

If installed as a napari plugin (entry point: `MitoDiffTrack:mitochondria_motility_analyzer`), you can launch it from:

```
Plugins → Mitochondria Motility Analyzer (python based plugin)
```

---

## 1. Launch the workflow in napari

### Option A. Run directly from Python

```
python MitoDiffTrack.py
```

A napari viewer will open with all panels added automatically.

### Option B. Launch from napari plugins menu

If installed as a plugin:

```
Plugins → Mitochondria Motility Analyzer (python based plugin)
```

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

All velocity measurements in exported tables will be in **µm per second**, based 이
