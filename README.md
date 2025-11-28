# MitoDiffTrack
A python based workflow for analysing mitochondrial motility in time lapse fluorescence imaging, implemented as a napari plugin with multiple dock widgets.

---

## 1. Launch the python based workflow in napari

1. Open **napari**.
2. Go to **Plugins → Mitochondria Motility Analyzer (python based plugin)**.
3. This will open several dock widgets:
   - `load_czi`
   - `seg_params`
   - `process_export`
   - `track_mito`
   - `batch_process_folder`

These panels should be used from top to bottom during a typical analysis.

---

## 2. Load the time lapse data and set calibration  
### *(load_czi panel)*

1. In the **load_czi** panel, select the processing mode:
   - **single** for analysing one file  
   - **folder** for batch processing of multiple files
2. Click the **vsi_file** or **czi_file** field and select the raw time lapse image.
3. Set the **spatial calibration**:
   - Pixel size (µm per pixel)
4. Set the **temporal calibration**:
   - Frame interval (seconds)
5. Press **Load image**.

The image sequence will appear as a new layer (e.g. `CZI_movie`).

**Note**  
Correct pixel size and frame interval are essential, because all velocities are reported in **µm per second**.

---

## 3. Adjust segmentation parameters and preview detection  
### *(seg_params panel)*

1. In the **seg_params** panel, choose the appropriate detection channel.
2. Set preprocessing filters, for example:
   - Gaussian blur (sigma)  
   - Optional background subtraction
3. Choose the threshold method and parameters:
   - Global or adaptive thresholding  
   - Threshold correction factor (sensitivity)
4. Set the expected object size range:
   - Minimum area (pixels)  
   - Maximum area (pixels)
5. Press **Preview segmentation** or **Update preview**.

A layer named `mito_preview_mask` will appear, showing detected mitochondria overlaid on the raw image.

**Tips**
- If small or dim mitochondria are missing → lower the threshold or reduce minimum area.  
- If too many noisy spots appear → raise the threshold or increase minimum area.

---

## 4. Run full segmentation and export morphology data  
### *(process_export panel)*

1. Open the **process_export** panel.
2. Choose an output folder.
3. Select the items to export, such as:
   - Label image for each frame  
   - Object measurements (area, length, aspect ratio, circularity, intensity)
4. Press **Run segmentation and export**.

The workflow processes all frames and creates:
- A `mito_labels` layer in napari  
- One or more CSV files with per object morphology measurements

---

## 5. Set LAP tracking parameters and build preview tracks  
### *(track_mito panel)*

1. Open the **track_mito** panel.
2. Set tracking parameters:
   - Maximum linking distance (µm)  
   - Maximum gap closing frames  
   - Minimum track length (frames)
3. Press **Build tracks preview**.

A `mito_tracks` layer will appear, showing trajectories overlaid on the movie.

Inspect key regions to verify that tracks follow the same mitochondrion over time, especially for fast or irregularly shaped mitochondria.

---

## 6. Export track data including velocities  
### *(track_mito panel)*

1. Stay in the **track_mito** panel.
2. Choose an output folder for the tracking results.
3. Press **Export track table**.

The workflow generates a CSV file containing:
- Track ID and object ID  
- Frame number and time (seconds)  
- X and Y position (µm)  
- Instantaneous displacement and velocity (µm per second)  
- Track level summary statistics (e.g., mean speed)

---

## 7. Optional batch analysis for multiple movies  
### *(batch_process_folder panel)*

1. Open the **batch_process_folder** panel.
2. Select an **input folder** containing multiple movies.
3. Select an **output folder**.
4. Reuse or copy segmentation and tracking parameters optimised on a test file.
5. Press **Run batch**.

All files in the folder will be processed automatically, and segmentation and tracking results will be saved for each movie.

---

## Citation
If you use this workflow in your research, please cite this repository.

---

## Contact
For questions or suggestions, please open an issue in this repository.
