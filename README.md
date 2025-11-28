# MitoDiffTrack

MitoDiffTrack is a python based workflow for analysing mitochondrial motility in time lapse fluorescence imaging.  
It is implemented as a napari based tool with several dock widgets that guide the user from loading data, through segmentation and tracking, to export of quantitative results.

Main panels:

- `load_czi`  
- `calibration` (Pixel / Time)  
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

Create a python environment and install the required packages, for example:

```bash
pip install napari[all] magicgui numpy pandas scikit-image scipy opencv-python aicspylibczi qtpy

You can then clone this repository and run:

python MitoDiffTrack_V01.py


This will open a napari viewer and automatically add all panels on the right side.

If MitoDiffTrack is later installed as a napari plugin with the entry point
MitoDiffTrack:mitochondria_motility_analyzer, it can also be launched from:

Plugins → Mitochondria Motility Analyzer (python based plugin)

1. Launch the python based workflow in napari
Option A. Run directly from python

In your environment, run:

python MitoDiffTrack.py


A napari viewer will open.

The following dock widgets appear on the right:

load_czi

calibration

seg_params

process_export (seg)

track_mito

batch_process_folder

plus advanced panels for motion filaments, fast tracks, export and GUI config.

Option B. From the napari Plugins menu

If installed as a plugin, in napari:

Go to Plugins → Mitochondria Motility Analyzer (python based plugin).

The same dock widgets will appear on the right.

2. Load the time lapse data and set calibration
load_czi and calibration panels

Load data

In the load_czi panel, click the button Load CZI / VSI.

Select a .czi file (or other supported time lapse file).

Choose the Z projection mode:

max, median, or z0

Optionally limit the number of frames for preview with Max preview T.

Press Load CZI / VSI.

A new image layer named CZI_movie will appear in the viewer.
For multichannel data, the channel axis is handled automatically.

Set pixel and time calibration

In the calibration panel:

Set Pixel size (µm/px)

Set Frame interval (s)

Press Apply.

These values are used to convert displacements into velocities in µm per second, both for global tracks and for fast tracks.

Note: Correct pixel size and frame interval are essential, because all velocities in the exported tables are reported in µm per second.

3. Adjust segmentation parameters and preview mitochondrial detection
seg_params panel

In the seg_params panel:

Select the detection Channel.

Set Downscale if you want to preview on a lower resolution copy for speed.

Configure preprocessing:

BG sigma for background smoothing and subtraction

Gaussian sigma for noise suppression

Optional CLAHE clip for contrast-limited local enhancement

Choose threshold parameters:

Thr mode: otsu, mean+std, percentile, or local

Thr correction, k (mean+std), Percentile or Local block size / offset

Set object size constraints:

Min area (px) and Max area (px)

Optional Fill holes, Opening radius, and Declump (watershed)

Press Preview segmentation at current T.

A mito_labels_preview layer appears, overlaying the detected mitochondria on the current time frame.

Tips:

If small or dim mitochondria are missing, lower the threshold or decrease Min area (px).

If too many noisy spots are included, raise the threshold or increase Min area (px).

4. Run full segmentation and build a per object table
process_export (seg) panel

In the process_export (seg) panel:

Set the detection Channel.

Define T start and T end (or keep them at -1 to use the full range).

Copy or fine tune the threshold and morphology parameters you optimised in the preview.

Press Run detection on all frames.

The workflow:

Runs segmentation across all selected frames.

Creates a mito_labels layer for the last processed frame.

Builds a global per object table stored in memory (LAST_OBJECTS_DF), which includes:

Frame index

Local object label

Centroid coordinates

Area and shape related measurements

Mean intensity

A status message reports how many objects were detected, and over which frame range.

5. Set LAP tracking parameters and build a track preview
track_mito panel

Open the track_mito panel.

Set linking and gating parameters:

Max disp (px/frame) for typical frame to frame displacement

Fast disp (px/frame) and Fast penalty factor for handling rare large jumps

Max area ratio and Max intensity ratio to prevent linking objects whose size or brightness changes unrealistically

Choose how tracks are displayed:

Color by: track_id or speed (if speed is available)

Preview tail length: how many recent time points to show in the track tail

Press Run LAP tracking from detections.

The workflow uses the global detection table to construct tracks:

Track steps and velocities are computed when pixel size and frame interval are set.

A mito_tracks_preview layer appears, showing tracks overlaid on the movie.

Inspect several regions and time points to verify that the tracks follow the same mitochondrion over time, including fast and irregularly shaped mitochondria.

6. Export track data and optional snapshot
process_export (CSV) panel

When you are satisfied with the segmentation and tracking:

Open the process_export (CSV) panel.

Set Output prefix, for example mito.

Choose what to export:

Save track snapshot (PNG)

Export global LAP per step

Export fast tracks (skeleton ROI) if you have run fast track analysis

Press Export CSV / snapshot.

Select an output folder.

The workflow writes:

Global per object table:

PREFIX_objects_YYYYMMDD_HHMMSS.csv

Global LAP per track summary:

PREFIX_LAP_tracks_YYYYMMDD_HHMMSS.csv

Includes duration, total displacement and speed statistics

Optional global LAP per step table:

PREFIX_LAP_steps_YYYYMMDD_HHMMSS.csv

Optional fast track tables:

PREFIX_FAST_tracks_YYYYMMDD_HHMMSS.csv

PREFIX_FAST_steps_YYYYMMDD_HHMMSS.csv

Optional PNG snapshot of the viewer:

PREFIX_tracks_snapshot_YYYYMMDD_HHMMSS.png

All velocities are reported in µm per second, based on the pixel size and frame interval set in the calibration panel.

7. Optional: Motion filaments and fast tracks
flow_motion_filaments and fast_tracks_skeleton panels

These panels provide an advanced way to highlight regions of sustained mitochondrial movement and to focus tracking on those regions.

7.1 Motion filaments

Run flow_motion_filaments:

Select Channel and T start / T end.

Set Downscale and Motion percentile.

Press Flow motion filaments.

The panel generates a label layer motion_filament_skel that highlights main paths of sustained motion.

7.2 Fast tracks inside motion skeleton ROI

Run fast_tracks_skeleton:

Choose Channel and time range.

Set ROI dilate to grow the skeleton into a thicker region of interest.

Configure more aggressive segmentation and linking parameters.

Press Fast track in skeleton ROI.

The panel:

Detects mitochondria only inside the motion skeleton region.

Performs tracking with parameters tuned for very fast movers.

Creates a mito_fast_tracks_preview layer.

Stores fast track steps in memory so they can be exported from process_export (CSV).

Fast tracks are summarised separately from global tracks and can be used to focus on high speed events.

8. Optional batch analysis of multiple movies
batch_process_folder panel

The batch_process_folder panel is included as a placeholder for future development.

It currently opens a simple dialog but does not yet execute full batch segmentation and tracking.

For this version, it is recommended to:

Optimise parameters on a representative test file.

Apply the same parameters to additional files manually using the single file workflow.

Optionally save and reload the GUI configuration as described below.

9. Save and load GUI configurations
save_gui_config and load_gui_config panels

To make analyses reproducible and easier to share:

Use save_gui_config to write all current panel parameters into a JSON file.

Use load_gui_config to restore parameter values later.

After loading a configuration, you can click Apply in the calibration panel once, to make sure the internal acquisition parameters are refreshed.

Notes and citation

Large raw movies are not stored in the repository.
Only analysis code and configuration are included.

This workflow is intended for research use in mitochondrial motility analysis.

If you use MitoDiffTrack in your work, please cite this repository and the corresponding method description in your manuscrip




