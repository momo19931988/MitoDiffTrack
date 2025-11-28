# MitoDiffTrack

MitoDiffTrack is a napari based pipeline for detecting and tracking mitochondrial motility in time lapse fluorescence movies. It combines CellProfiler style object detection, global LAP tracking, and motion based skeleton regions of interest to highlight fast moving mitochondria and export quantitative motility statistics in µm per second. :contentReference[oaicite:0]{index=0}

---

## Key features

- Load multi dimensional Zeiss CZI movies and perform Z projection
- CellProfiler like IdentifyPrimaryObjects style segmentation
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

You can create a virtual environment and install dependencies, for example:

```bash
pip install napari[all] magicgui numpy pandas scikit-image scipy opencv-python aicspylibczi qtpy
