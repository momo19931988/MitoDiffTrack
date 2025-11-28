# ============================================
# Script name: MitoDiffTrack.py
#
# Purpose:
#   - Napari-based workflow "Mitochondria Motility Analyzer (python based plugin)"
#   - Main panels (dock widgets), matching the README:
#       * load_czi             – load time-lapse data and set calibration
#       * seg_params           – adjust segmentation parameters and preview
#       * process_export       – run full segmentation / tracking export
#       * track_mito           – LAP-based mitochondrial tracking
#       * batch_process_folder – (placeholder) batch processing of multiple movies
#   - Additional advanced panels:
#       * Flow motion filaments + fast tracks in skeleton ROI
#       * GUI config save / load
#
# Author: Mo
# Last modified: 2025-11-28 16:05
# ============================================

import os
os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")
os.environ.setdefault("QT_ENABLE_HIGHDPI_SCALING", "1")
os.environ.setdefault("QT_SCALE_FACTOR_ROUNDING_POLICY", "PassThrough")
os.environ.setdefault("NUMEXPR_MAX_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple
import json

import numpy as np
import pandas as pd

import napari
from napari.utils.notifications import show_info, show_error
from napari.utils import progress
from magicgui import magicgui
from qtpy.QtWidgets import QFileDialog
from qtpy import QtWidgets

from aicspylibczi import CziFile
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment

from skimage.filters import (
    threshold_otsu,
    gaussian,
    threshold_local,
)
from skimage.morphology import (
    remove_small_objects,
    remove_small_holes,
    binary_opening,
    disk,
    closing,
    skeletonize,
)
from skimage.segmentation import watershed
from skimage.measure import label, regionprops_table
from skimage.util import img_as_float32
from skimage import exposure

import cv2  # for dense optical flow


# --------------------------------------------
# Global caches
# --------------------------------------------
LAST_OBJECTS_DF: pd.DataFrame | None = None          # global per-object detection
LAST_TRACK_STEPS_DF: pd.DataFrame | None = None      # global LAP per-step
LAST_META: dict = {}                                 # pixel size, frame dt, downscale, t range etc

# Fast-track caches (skeleton ROI)
FAST_OBJECTS_DF: pd.DataFrame | None = None          # per-object inside skeleton ROI
FAST_TRACK_STEPS_DF: pd.DataFrame | None = None      # fast-track per-step


# ============================================
# CZI reading and Z projection
# ============================================

def _size_from_dim(v, default=1):
    """Helper to get size from CZI dims."""
    if isinstance(v, (tuple, list)):
        for x in v:
            try:
                x = int(x)
                if x > 0:
                    return x
            except Exception:
                continue
        return int(default)
    try:
        v = int(v)
        return v if v > 0 else int(default)
    except Exception:
        return int(default)


def _present_sizes(dims):
    """Parse dimension information from CziFile.get_dims_shape()."""
    d = dims[0] if isinstance(dims, list) and dims and isinstance(dims[0], dict) else dims
    if not isinstance(d, dict):
        return set(), 1, 1, 1, 1
    present = set(d.keys())
    S = _size_from_dim(d.get("S", 1), 1)
    T = _size_from_dim(d.get("T", 1), 1)
    C = _size_from_dim(d.get("C", 1), 1)
    Z = _size_from_dim(d.get("Z", 1), 1)
    return present, S, T, C, Z


def _read_plane(czi: CziFile, present: set, t=None, z=None, c=None, s=None):
    """Read a single plane from CZI with given indices."""
    kwargs = {}
    if "T" in present and t is not None:
        kwargs["T"] = int(t)
    if "Z" in present and z is not None:
        kwargs["Z"] = int(z)
    if "C" in present and c is not None:
        kwargs["C"] = int(c)
    if "S" in present and s is not None:
        kwargs["S"] = int(s)

    sb = czi.read_image(**kwargs)
    arr = np.asarray(sb[0])
    if arr.ndim == 3 and arr.shape[-1] in (3, 4):
        arr = arr[..., :3].mean(axis=-1)
    return arr


def _z_project_stack(path: str, z_method: str = "max"):
    """
    Read a CZI file and perform Z projection if needed.
    Return array of shape (T, C, Y, X).
    """
    czi = CziFile(path)
    dims = czi.get_dims_shape()
    present, S, T, C, Z = _present_sizes(dims)

    sample = _read_plane(
        czi,
        present,
        t=0 if "T" in present else None,
        z=0 if "Z" in present else None,
        c=0 if "C" in present else None,
        s=0 if "S" in present else None,
    )
    Y, X = sample.shape[-2], sample.shape[-1]
    dtype = sample.dtype

    data = np.zeros((max(T, 1), max(C, 1), Y, X), dtype=dtype)
    use_z = ("Z" in present and Z > 1 and z_method in ("max", "median"))

    for t in range(max(T, 1)):
        for c in range(max(C, 1)):
            if use_z:
                planes = []
                for z in range(Z):
                    planes.append(
                        _read_plane(
                            czi,
                            present,
                            t=t,
                            z=z,
                            c=c,
                            s=0 if "S" in present else None,
                        )
                    )
                stack = np.stack(planes, axis=0)
                if z_method == "median":
                    arr = np.median(stack, axis=0)
                else:
                    arr = np.max(stack, axis=0)
            else:
                zsel = 0 if ("Z" in present and Z > 0 and z_method == "z0") else None
                arr = _read_plane(
                    czi,
                    present,
                    t=t if "T" in present else None,
                    z=zsel,
                    c=c if "C" in present else None,
                    s=0 if "S" in present else None,
                )
            if arr.ndim != 2:
                arr = np.squeeze(arr)[-Y:, -X:]
            data[t, c] = arr

    info = {"T": int(T), "C": int(C), "Z": int(Z), "S": int(S), "present": sorted(list(present))}
    return data, info


# ============================================
# Helpers
# ============================================

def _get_layer_by_name(viewer: napari.Viewer, name: str):
    """Safely get a layer by name."""
    try:
        return viewer.layers[name]
    except KeyError:
        for lyr in viewer.layers:
            if getattr(lyr, "name", "") == name:
                return lyr
    return None


def _current_t(viewer: napari.Viewer) -> int:
    """Return current time index of viewer."""
    try:
        return int(viewer.dims.current_step[0])
    except Exception:
        return 0


def _get_or_create_labels_layer(viewer: napari.Viewer, name: str, shape, opacity=0.6):
    """Create or reuse a labels layer with given name."""
    try:
        lyr = viewer.layers[name]
    except KeyError:
        lyr = None

    if lyr is None:
        zeros = np.zeros(shape, dtype=np.int32)
        lyr = viewer.add_labels(zeros, name=name, opacity=opacity)
    try:
        lyr.visible = True
        idx = list(viewer.layers).index(lyr)
        viewer.layers.move(idx, len(viewer.layers) - 1)
    except Exception:
        pass
    return lyr


def _normalize_01(img: np.ndarray) -> np.ndarray:
    """Normalize image to [0, 1]."""
    img = img.astype(np.float32)
    vmin = float(img.min())
    vmax = float(img.max())
    if vmax > vmin:
        img = (img - vmin) / (vmax - vmin)
    else:
        img[:] = 0.0
    return img


# ============================================
# load_czi panel: Load CZI movie into Napari
# ============================================

@magicgui(
    call_button="Load CZI / VSI",
    z_projection={"label": "Z-projection", "choices": ["max", "median", "z0"], "value": "max"},
    single_channel_only={"label": "Load only current channel", "widget_type": "CheckBox", "value": True},
    max_preview_T={"label": "Max preview T (0=all)", "min": 0, "max": 5000, "step": 1, "value": 0},
)
def load_czi(viewer: napari.Viewer,
             z_projection: str = "max",
             single_channel_only: bool = True,
             max_preview_T: int = 0):
    """
    load_czi panel (README step: "Load the time lapse data and set calibration"):

      - Select .czi (or .vsi via CZI-compatible reader)
      - Load as CZI_movie layer
      - Combined with acquisition_params panel for pixel size & frame interval
    """
    try:
        path, _ = QFileDialog.getOpenFileName(
            None,
            "Select .czi",
            "",
            "Zeiss CZI (*.czi);;All files (*.*)",
        )
    except Exception as e:
        show_error(f"Open dialog failed: {e}")
        return
    if not path:
        show_info("No file selected.")
        return

    try:
        data, info = _z_project_stack(path, z_method=z_projection)
    except Exception as e:
        show_error(f"Failed to read CZI: {e}")
        return

    T_all = data.shape[0]
    if single_channel_only and data.ndim == 4:
        data = data[:, 0:1, ...]

    if max_preview_T and max_preview_T > 0 and T_all > max_preview_T:
        data = data[:max_preview_T, ...]

    layer = viewer.add_image(
        data,
        name="CZI_movie",
        channel_axis=1 if data.ndim == 4 else None,
        rgb=False,
        blending="additive",
    )
    try:
        if hasattr(layer, "metadata") and isinstance(layer.metadata, dict):
            layer.metadata.update(info)
            layer.metadata["source_path"] = path
    except Exception:
        pass

    viewer.dims.axis_labels = ["T", "C", "Y", "X"] if data.ndim == 4 else ["T", "Y", "X"]
    show_info(f"Loaded: {Path(path).name} - T={info.get('T')} C={info.get('C')}")


# ============================================
# seg_params: Preprocess and identify primary objects
# ============================================

def _preprocess_intensity(img: np.ndarray,
                          bg_sigma: float = 0.0,
                          smooth_sigma: float = 1.0,
                          clahe_clip: float = 0.0) -> np.ndarray:
    """
    Mito-friendly intensity preprocessing.
    """
    a = img_as_float32(img)

    if bg_sigma and bg_sigma > 0:
        bg = gaussian(a, sigma=float(bg_sigma), preserve_range=True)
        a = a - bg

    a = a - a.min()
    max_val = a.max()
    if max_val > 0:
        a /= max_val

    if clahe_clip and clahe_clip > 0:
        a = exposure.equalize_adapthist(a, clip_limit=float(clahe_clip))

    if smooth_sigma and smooth_sigma > 0:
        a = gaussian(a, sigma=float(smooth_sigma), preserve_range=True)

    return a


def _identify_primary_objects(img: np.ndarray,
                              thr_mode: str = "otsu",
                              thr_corr: float = 1.0,
                              thr_k: float = 1.0,
                              thr_percentile: float = 95.0,
                              smooth_sigma: float = 1.0,
                              bg_sigma: float = 0.0,
                              clahe_clip: float = 0.0,
                              local_block_size: int = 51,
                              local_offset: float = 0.0,
                              min_area_px: int = 10,
                              max_area_px: int = 999999,
                              fill_holes_flag: bool = True,
                              opening_radius: int = 0,
                              declump_watershed_flag: bool = True,
                              peak_min_distance: int = 3,
                              filament_mode: bool = False,
                              close_radius: int = 0) -> np.ndarray:
    """
    Identify primary objects (mitochondria) with flexible thresholding.
    """
    a = _preprocess_intensity(
        img,
        bg_sigma=bg_sigma,
        smooth_sigma=smooth_sigma,
        clahe_clip=clahe_clip,
    )

    thr_mode = str(thr_mode).lower()
    if thr_mode == "local":
        bs = max(3, int(local_block_size))
        if bs % 2 == 0:
            bs += 1
        th_map = threshold_local(a, block_size=bs, offset=float(local_offset))
        th_map = th_map * float(thr_corr)
        mask = a > th_map
    else:
        if thr_mode == "mean+std":
            m = float(a.mean())
            s = float(a.std())
            thr_base = m + float(thr_k) * s
        elif thr_mode == "percentile":
            thr_base = float(np.percentile(a, float(thr_percentile)))
        else:
            thr_base = float(threshold_otsu(a))
        thr = thr_base * float(thr_corr)
        mask = a > thr

    if opening_radius and opening_radius > 0:
        mask = binary_opening(mask, footprint=disk(int(opening_radius)))

    if fill_holes_flag:
        mask = remove_small_holes(mask, area_threshold=int(max(1, min_area_px)))

    mask = remove_small_objects(mask, min_size=int(min_area_px))
    if max_area_px < 1e12:
        lbl_tmp = label(mask, connectivity=1)
        props = regionprops_table(lbl_tmp, properties=["label", "area"])
        keep = set(
            int(l)
            for l, area in zip(props["label"], props["area"])
            if area <= max_area_px
        )
        mask = np.isin(lbl_tmp, list(keep))

    if filament_mode and close_radius and close_radius > 0:
        mask = closing(mask, footprint=disk(int(close_radius)))

    if filament_mode:
        labels_out = label(mask, connectivity=1)
    else:
        if declump_watershed_flag:
            dist = ndi.distance_transform_edt(mask)
            from skimage.feature import peak_local_max
            coords = peak_local_max(
                dist,
                min_distance=max(1, int(peak_min_distance)),
                labels=mask,
            )
            markers = np.zeros_like(mask, dtype=int)
            for i, (r, c) in enumerate(coords, start=1):
                markers[r, c] = i
            if markers.max() == 0:
                labels_out = label(mask, connectivity=1)
            else:
                labels_out = watershed(-dist, markers, mask=mask)
        else:
            labels_out = label(mask, connectivity=1)

    return labels_out.astype(np.int32, copy=False)


def _centroids_from_labels(lbl: np.ndarray, intensity_img: np.ndarray | None = None) -> pd.DataFrame:
    """
    Compute regionprops for labels.
    """
    if lbl.max() == 0:
        cols = [
            "label", "cy", "cx", "area",
            "major_axis_length", "minor_axis_length",
            "eccentricity", "perimeter", "mean_intensity",
        ]
        return pd.DataFrame(columns=cols)

    props = regionprops_table(
        lbl,
        intensity_image=intensity_img,
        properties=(
            "label",
            "centroid",
            "area",
            "major_axis_length",
            "minor_axis_length",
            "eccentricity",
            "perimeter",
            "mean_intensity",
        ),
    )
    df = pd.DataFrame(props).rename(columns={"centroid-0": "cy", "centroid-1": "cx"})
    return df


# seg_params panel = segmentation preview
@magicgui(
    call_button="Preview segmentation at current T",
    channel={"label": "Channel", "min": 0, "max": 63, "step": 1, "value": 0},
    downscale={"label": "Downscale", "choices": [1, 2, 4, 8], "value": 1},

    thr_mode={"label": "Thr mode",
              "choices": ["otsu", "mean+std", "percentile", "local"],
              "value": "otsu"},
    thr_corr={"label": "Thr correction", "min": 0.1, "max": 3.0, "step": 0.05, "value": 1.0},
    thr_k={"label": "k (mean+std)", "min": 0.0, "max": 5.0, "step": 0.1, "value": 1.0},
    thr_percentile={"label": "Percentile", "min": 50.0, "max": 99.9, "step": 0.1, "value": 95.0},
    local_block_size={"label": "Local block size", "min": 3, "max": 201, "step": 2, "value": 51},
    local_offset={"label": "Local offset", "min": -0.5, "max": 0.5, "step": 0.01, "value": 0.0},

    bg_sigma={"label": "BG sigma", "min": 0.0, "max": 50.0, "step": 0.5, "value": 8.0},
    clahe_clip={"label": "CLAHE clip (0=off)", "min": 0.0, "max": 0.1, "step": 0.005, "value": 0.015},
    smooth_sigma={"label": "Gaussian sigma", "min": 0.0, "max": 5.0, "step": 0.2, "value": 1.0},

    min_area_px={"label": "Min area (px)", "min": 1, "max": 50000, "step": 1, "value": 20},
    max_area_px={"label": "Max area (px)", "min": 10, "max": 1000000, "step": 10, "value": 5000},
    fill_holes_flag={"label": "Fill holes", "widget_type": "CheckBox", "value": True},
    opening_radius={"label": "Opening radius (px)", "min": 0, "max": 10, "step": 1, "value": 0},
    declump_watershed_flag={"label": "Declump (watershed)", "widget_type": "CheckBox", "value": True},
    peak_min_distance={"label": "Peak min distance", "min": 1, "max": 20, "step": 1, "value": 3},
    filament_mode={"label": "Filament mode", "widget_type": "CheckBox", "value": False},
    close_radius={"label": "Closing radius (px)", "min": 0, "max": 10, "step": 1, "value": 0},
)
def identify_preview(viewer: napari.Viewer,
                     channel: int = 0,
                     downscale: int = 1,
                     thr_mode: str = "otsu",
                     thr_corr: float = 1.0,
                     thr_k: float = 1.0,
                     thr_percentile: float = 95.0,
                     local_block_size: int = 51,
                     local_offset: float = 0.0,
                     bg_sigma: float = 8.0,
                     clahe_clip: float = 0.015,
                     smooth_sigma: float = 1.0,
                     min_area_px: int = 20,
                     max_area_px: int = 5000,
                     fill_holes_flag: bool = True,
                     opening_radius: int = 0,
                     declump_watershed_flag: bool = True,
                     peak_min_distance: int = 3,
                     filament_mode: bool = False,
                     close_radius: int = 0):
    """
    seg_params panel (README step: "Adjust segmentation parameters and preview"):

      - Adjust preprocessing and threshold parameters
      - Preview segmentation on the current time frame
      - A mito_labels_preview layer will appear
    """
    base = _get_layer_by_name(viewer, "CZI_movie")
    if base is None or getattr(base, "data", None) is None:
        show_error("No base image. Please load CZI first.")
        return
    data = np.asarray(base.data)
    t = _current_t(viewer)

    if data.ndim == 3:
        frame = data[t]
    else:
        T, C, Y, X = data.shape
        if channel >= C:
            channel = 0
        frame = data[t, channel]

    if frame.ndim == 3:
        frame = np.max(frame, axis=0)
    ds = max(1, int(downscale))
    img = frame[::ds, ::ds]

    labels = _identify_primary_objects(
        img,
        thr_mode=thr_mode,
        thr_corr=thr_corr,
        thr_k=thr_k,
        thr_percentile=thr_percentile,
        smooth_sigma=smooth_sigma,
        bg_sigma=bg_sigma,
        clahe_clip=clahe_clip,
        local_block_size=local_block_size,
        local_offset=local_offset,
        min_area_px=min_area_px,
        max_area_px=max_area_px,
        fill_holes_flag=fill_holes_flag,
        opening_radius=opening_radius,
        declump_watershed_flag=declump_watershed_flag,
        peak_min_distance=peak_min_distance,
        filament_mode=filament_mode,
        close_radius=close_radius,
    )

    lyr = _get_or_create_labels_layer(viewer, "mito_labels_preview", labels.shape, opacity=0.6)
    lyr.data = labels


# ============================================
# Global acquisition parameters (calibration)
# ============================================

@magicgui(
    call_button="Apply",
    pixel_size_um={"label": "Pixel size (µm/px)", "min": 0.0, "max": 10.0, "step": 0.001, "value": 0.11},
    frame_dt_s={"label": "Frame interval (s)", "min": 0.0, "max": 3600.0, "step": 0.001, "value": 5.0},
)
def acquisition_params(viewer: napari.Viewer,
                       pixel_size_um: float = 0.11,
                       frame_dt_s: float = 5.0):
    """
    Calibration panel (README: part of load_czi step):

      - Set pixel size and frame interval
      - These values are required to compute speed in µm/s
    """
    global LAST_META
    LAST_META["pixel_size_um"] = float(pixel_size_um) if pixel_size_um > 0 else None
    LAST_META["frame_dt_s"] = float(frame_dt_s) if frame_dt_s > 0 else None
    show_info(f"Acquisition params updated: pixel={pixel_size_um} µm/px, dt={frame_dt_s} s")


# ============================================
# Stage 1 - global detection for all frames
# (used by process_export step)
# ============================================

@magicgui(
    call_button="Run detection on all frames",
    channel={"label": "Channel", "min": 0, "max": 63, "step": 1, "value": 0},
    t_start={"label": "T start (-1=auto)", "min": -1, "max": 999999, "step": 1, "value": -1},
    t_end={"label": "T end (-1=auto)", "min": -1, "max": 999999, "step": 1, "value": -1},
    downscale={"label": "Downscale", "choices": [1, 2, 4, 8], "value": 1},

    thr_mode={"label": "Thr mode",
              "choices": ["otsu", "mean+std", "percentile", "local"],
              "value": "otsu"},
    thr_corr={"label": "Thr correction", "min": 0.1, "max": 3.0, "step": 0.05, "value": 1.15},
    thr_k={"label": "k (mean+std)", "min": 0.0, "max": 5.0, "step": 0.1, "value": 1.0},
    thr_percentile={"label": "Percentile", "min": 50.0, "max": 99.9, "step": 0.1, "value": 95.0},
    local_block_size={"label": "Local block size", "min": 3, "max": 201, "step": 2, "value": 51},
    local_offset={"label": "Local offset", "min": -0.5, "max": 0.5, "step": 0.01, "value": 0.0},

    bg_sigma={"label": "BG sigma", "min": 0.0, "max": 50.0, "step": 0.5, "value": 8.0},
    clahe_clip={"label": "CLAHE clip (0=off)", "min": 0.0, "max": 0.1, "step": 0.005, "value": 0.015},
    smooth_sigma={"label": "Gaussian sigma", "min": 0.0, "max": 5.0, "step": 0.2, "value": 1.0},

    min_area_px={"label": "Min area (px)", "min": 1, "max": 50000, "step": 1, "value": 20},
    max_area_px={"label": "Max area (px)", "min": 10, "max": 1000000, "step": 10, "value": 5000},
    fill_holes_flag={"label": "Fill holes", "widget_type": "CheckBox", "value": True},
    opening_radius={"label": "Opening radius (px)", "min": 0, "max": 10, "step": 1, "value": 0},
    declump_watershed_flag={"label": "Declump (watershed)", "widget_type": "CheckBox", "value": True},
    peak_min_distance={"label": "Peak min distance", "min": 1, "max": 20, "step": 1, "value": 3},
    filament_mode={"label": "Filament mode", "widget_type": "CheckBox", "value": False},
    close_radius={"label": "Closing radius (px)", "min": 0, "max": 10, "step": 1, "value": 0},
)
def detect_all_frames(viewer: napari.Viewer,
                      channel: int = 0,
                      t_start: int = -1,
                      t_end: int = -1,
                      downscale: int = 1,
                      thr_mode: str = "otsu",
                      thr_corr: float = 1.15,
                      thr_k: float = 1.0,
                      thr_percentile: float = 95.0,
                      local_block_size: int = 51,
                      local_offset: float = 0.0,
                      bg_sigma: float = 8.0,
                      clahe_clip: float = 0.015,
                      smooth_sigma: float = 1.0,
                      min_area_px: int = 20,
                      max_area_px: int = 5000,
                      fill_holes_flag: bool = True,
                      opening_radius: int = 0,
                      declump_watershed_flag: bool = True,
                      peak_min_distance: int = 3,
                      filament_mode: bool = False,
                      close_radius: int = 0):
    """
    process_export (segmentation part):

      - Run intensity-based object detection on all frames
      - Fill LAST_OBJECTS_DF (per-object table) used for morphology export and tracking
    """
    global LAST_OBJECTS_DF, LAST_META

    base = _get_layer_by_name(viewer, "CZI_movie")
    if base is None or getattr(base, "data", None) is None:
        show_error("No base image. Please load CZI first.")
        return
    data = np.asarray(base.data)

    if data.ndim == 3:
        T, Y, X = data.shape

        def get_frame(tt):
            return data[tt]
    else:
        T, C, Y, X = data.shape
        if channel >= C:
            show_error(f"Channel {channel} out of range: C={C}")
            return

        def get_frame(tt):
            return data[tt, channel]

    if t_start < 0:
        t0 = 0
    else:
        t0 = max(0, min(t_start, T - 1))
    if t_end < 0:
        t1 = T - 1
    else:
        t1 = max(0, min(t_end, T - 1))

    if t1 <= t0:
        show_error("Need at least 2 frames (t_end > t_start).")
        return

    ds = max(1, int(downscale))
    obj_rows: List[Dict] = []

    lyr_preview = _get_layer_by_name(viewer, "mito_labels")
    if lyr_preview is not None:
        try:
            viewer.layers.remove(lyr_preview)
        except Exception:
            pass

    for t in progress(range(t0, t1 + 1), desc="Detect objects (global)"):
        frame = get_frame(t)
        if frame.ndim == 3:
            frame = np.max(frame, axis=0)
        img = frame[::ds, ::ds]

        labels = _identify_primary_objects(
            img,
            thr_mode=thr_mode,
            thr_corr=thr_corr,
            thr_k=thr_k,
            thr_percentile=thr_percentile,
            smooth_sigma=smooth_sigma,
            bg_sigma=bg_sigma,
            clahe_clip=clahe_clip,
            local_block_size=local_block_size,
            local_offset=local_offset,
            min_area_px=min_area_px,
            max_area_px=max_area_px,
            fill_holes_flag=fill_holes_flag,
            opening_radius=opening_radius,
            declump_watershed_flag=declump_watershed_flag,
            peak_min_distance=peak_min_distance,
            filament_mode=filament_mode,
            close_radius=close_radius,
        )
        props_df = _centroids_from_labels(labels, intensity_img=img)
        if props_df.empty:
            continue
        props_df = props_df.reset_index(drop=True)
        for _, r in props_df.iterrows():
            obj_rows.append(
                {
                    "t": int(t),
                    "obj_local": int(r["label"]),
                    "cy": float(r["cy"]),
                    "cx": float(r["cx"]),
                    "area_px": float(r["area"]),


                    "major_axis": float(r["major_axis_length"]),
                    "minor_axis": float(r["minor_axis_length"]),
                    "eccentricity": float(r["eccentricity"]),
                    "perimeter": float(r["perimeter"]),
                    "mean_intensity": float(r["mean_intensity"]),
                }
            )

        if t == t1:
            lyr = _get_or_create_labels_layer(viewer, "mito_labels", labels.shape, opacity=0.5)
            lyr.data = labels

    if not obj_rows:
        show_error("No objects detected in the given time range.")
    else:
        LAST_OBJECTS_DF = pd.DataFrame(obj_rows)
        LAST_OBJECTS_DF = LAST_OBJECTS_DF.sort_values(["t", "cy", "cx"]).reset_index(drop=True)
        LAST_OBJECTS_DF["obj_id"] = LAST_OBJECTS_DF.index.astype(int)

        LAST_META["downscale"] = ds
        LAST_META["t_start"] = int(t0)
        LAST_META["t_end"] = int(t1)

        show_info(f"Detection finished. Objects: {len(LAST_OBJECTS_DF)} from t={t0} to t={t1}.")


# ============================================
# Stage 2 - LAP tracking (global)
# track_mito panel
# ============================================

def _build_tracks_by_LAP(objects_df: pd.DataFrame,
                         max_disp_px: float,
                         pixel_size_um: float | None,
                         frame_dt_s: float | None,
                         area_ratio_max: float = 3.0,
                         intensity_ratio_max: float = 3.0,
                         fast_disp_px: float | None = None,
                         fast_penalty_factor: float = 5.0) -> Tuple[pd.DataFrame, List[List[float]], List[Dict]]:
    """
    Build tracks using a simplified LAP between consecutive frames.
    """
    df = objects_df.copy()
    if "obj_id" not in df.columns:
        df["obj_id"] = df.index.astype(int)

    t_all = sorted(df["t"].unique().tolist())
    if len(t_all) < 2:
        raise RuntimeError("Need at least 2 frames for tracking.")

    frames = {int(t): sub_df.reset_index(drop=True) for t, sub_df in df.groupby("t")}

    assignment: Dict[Tuple[int, int], int] = {}
    next_track_id = 1

    steps_rows: List[Dict] = []
    tracks_points: List[List[float]] = []
    props_rows: List[Dict] = []

    big_cost = 1e6
    max_disp_px = float(max_disp_px) if max_disp_px > 0 else float("inf")
    if fast_disp_px is None or fast_disp_px <= 0 or fast_disp_px <= max_disp_px:
        fast_disp_px_val = max_disp_px
    else:
        fast_disp_px_val = float(fast_disp_px)
    fast_penalty_factor = max(1.0, float(fast_penalty_factor))
    eps = 1e-6

    for t in progress(range(min(t_all), max(t_all)), desc="LAP tracking"):
        if t not in frames or (t + 1) not in frames:
            continue
        df0 = frames[t]
        df1 = frames[t + 1]
        if df0.empty or df1.empty:
            continue

        pts0 = df0[["cy", "cx"]].to_numpy(float)
        pts1 = df1[["cy", "cx"]].to_numpy(float)
        n0, n1 = pts0.shape[0], pts1.shape[0]
        if n0 == 0 or n1 == 0:
            continue

        diff = pts0[:, None, :] - pts1[None, :, :]
        dist = np.sqrt(np.sum(diff ** 2, axis=-1))

        cost = dist.copy()
        fast_mask = np.zeros_like(dist, dtype=bool)

        if np.isfinite(max_disp_px) or np.isfinite(fast_disp_px_val):
            mask_far = dist > fast_disp_px_val
            cost[mask_far] = big_cost

            if fast_disp_px_val > max_disp_px:
                mask_fast_band = (dist > max_disp_px) & (dist <= fast_disp_px_val)
                cost[mask_fast_band] = dist[mask_fast_band] * fast_penalty_factor
                fast_mask[mask_fast_band] = True

        if "area_px" in df0.columns and "area_px" in df1.columns:
            a0 = df0["area_px"].to_numpy(float)[:, None]
            a1 = df1["area_px"].to_numpy(float)[None, :]
            ratio_a = (a0 + eps) / (a1 + eps)
            invalid_area = (ratio_a > float(area_ratio_max)) | (ratio_a < 1.0 / float(area_ratio_max))
            cost[invalid_area] = big_cost
            fast_mask[invalid_area] = False

        if "mean_intensity" in df0.columns and "mean_intensity" in df1.columns:
            i0 = df0["mean_intensity"].to_numpy(float)[:, None]
            i1 = df1["mean_intensity"].to_numpy(float)[None, :]
            ratio_i = (i0 + eps) / (i1 + eps)
            invalid_int = (ratio_i > float(intensity_ratio_max)) | (ratio_i < 1.0 / float(intensity_ratio_max))
            cost[invalid_int] = big_cost
            fast_mask[invalid_int] = False

        row_ind, col_ind = linear_sum_assignment(cost)
        for i, j in zip(row_ind, col_ind):
            c_ij = cost[i, j]
            if c_ij >= big_cost:
                continue

            r0 = df0.loc[i]
            r1 = df1.loc[j]

            key0 = (int(t), int(r0["obj_id"]))
            key1 = (int(t + 1), int(r1["obj_id"]))

            if key0 in assignment:
                tid = assignment[key0]
            else:
                tid = next_track_id
                next_track_id += 1
            assignment[key0] = tid
            assignment[key1] = tid

            cy0, cx0 = float(r0["cy"]), float(r0["cx"])
            cy1, cx1 = float(r1["cy"]), float(r1["cx"])
            disp_px = float(dist[i, j])
            is_fast = bool(fast_mask[i, j])

            step = {
                "track_id": tid,
                "t0": int(t),
                "t1": int(t + 1),
                "cy0": cy0,
                "cx0": cx0,
                "cy1": cy1,
                "cx1": cx1,
                "disp_px": disp_px,
                "fast_jump": is_fast,
            }

            if pixel_size_um is not None and frame_dt_s is not None and frame_dt_s > 0:
                disp_um = disp_px * float(pixel_size_um)
                speed = disp_um / float(frame_dt_s)
                step["disp_um"] = disp_um
                step["speed_um_per_s"] = speed

            steps_rows.append(step)

            tracks_points.append([tid, int(t), cy0, cx0])
            tracks_points.append([tid, int(t + 1), cy1, cx1])

            props_rows.append({"track_id": tid})
            props_rows.append({"track_id": tid})

    if not steps_rows:
        raise RuntimeError("No track steps found. Try increasing max_disp_px or fast_disp_px, or relax gating.")

    steps_df = pd.DataFrame(steps_rows)

    if "speed_um_per_s" in steps_df.columns:
        speed_map = steps_df.groupby("track_id")["speed_um_per_s"].median().to_dict()
        new_props = []
        for p in props_rows:
            tid = int(p["track_id"])
            new_props.append(
                {
                    "track_id": tid,
                    "speed_um_per_s": float(speed_map.get(tid, np.nan)),
                }
            )
        props_rows = new_props

    return steps_df, tracks_points, props_rows


def _update_tracks_preview(viewer: napari.Viewer,
                           tracks_points: List[List[float]],
                           props_rows: List[Dict],
                           color_by: str,
                           tail_len: int,
                           layer_name: str):
    """
    Update or create preview tracks layer.
    """
    if not tracks_points:
        return

    tracks_arr = np.asarray(tracks_points, dtype=float)
    props_df = pd.DataFrame(props_rows) if props_rows else pd.DataFrame(index=range(len(tracks_arr)))
    if "track_id" not in props_df.columns:
        props_df["track_id"] = tracks_arr[:, 0].astype(int)

    lyr = _get_layer_by_name(viewer, layer_name)
    if lyr is None:
        lyr = viewer.add_tracks(
            data=tracks_arr,
            properties=props_df.to_dict(orient="list"),
            name=layer_name,
            tail_length=int(tail_len),
            head_length=0,
        )
    else:
        lyr.data = tracks_arr
        lyr.properties = props_df.to_dict(orient="list")
        try:
            lyr.tail_length = int(tail_len)
        except Exception:
            pass

    try:
        if color_by == "speed" and "speed_um_per_s" in props_df.columns:
            lyr.color_by = "speed_um_per_s"
        else:
            lyr.color_by = "track_id"
    except Exception:
        pass


@magicgui(
    call_button="Run LAP tracking from detections",
    max_disp_px={"label": "Max disp (px/frame)", "min": 0.0, "max": 200.0, "step": 0.5, "value": 5.0},
    fast_disp_px={"label": "Fast disp (px/frame)", "min": 0.0, "max": 400.0, "step": 0.5, "value": 5.0},
    fast_penalty_factor={"label": "Fast penalty factor", "min": 1.0, "max": 50.0, "step": 0.5, "value": 5.0},
    area_ratio_max={"label": "Max area ratio", "min": 1.0, "max": 10.0, "step": 0.1, "value": 3.0},
    intensity_ratio_max={"label": "Max intensity ratio", "min": 1.0, "max": 10.0, "step": 0.1, "value": 3.0},
    color_by={"label": "Color by", "choices": ["track_id", "speed"], "value": "track_id"},
    preview_tail={"label": "Preview tail length", "min": 0, "max": 500, "step": 10, "value": 100},
)
def lap_tracking(viewer: napari.Viewer,
                 max_disp_px: float = 5.0,
                 fast_disp_px: float = 5.0,
                 fast_penalty_factor: float = 5.0,
                 area_ratio_max: float = 3.0,
                 intensity_ratio_max: float = 3.0,
                 color_by: str = "track_id",
                 preview_tail: int = 100):
    """
    track_mito panel (README step: "Set LAP tracking parameters and build a track preview"):

      - Use LAST_OBJECTS_DF from detection
      - Run LAP tracking
      - Create mito_tracks_preview layer
    """
    global LAST_OBJECTS_DF, LAST_TRACK_STEPS_DF, LAST_META

    if LAST_OBJECTS_DF is None or LAST_OBJECTS_DF.empty:
        show_error("No detection result. Please run 'Run detection on all frames' first.")
        return

    p_um = LAST_META.get("pixel_size_um", None)
    dt_s = LAST_META.get("frame_dt_s", None)
    ds = LAST_META.get("downscale", 1)
    if ds is None or ds < 1:
        ds = 1

    if p_um is None or dt_s is None:
        show_info("Pixel size or frame interval not set. µm/s cannot be computed.")
        p_eff = None
    else:
        p_eff = float(p_um) * float(ds)

    try:
        steps_df, tracks_points, props_rows = _build_tracks_by_LAP(
            LAST_OBJECTS_DF,
            max_disp_px=max_disp_px,
            pixel_size_um=p_eff,
            frame_dt_s=dt_s,
            area_ratio_max=area_ratio_max,
            intensity_ratio_max=intensity_ratio_max,
            fast_disp_px=fast_disp_px,
            fast_penalty_factor=fast_penalty_factor,
        )
    except Exception as e:
        show_error(f"LAP tracking error: {e}")
        return

    LAST_TRACK_STEPS_DF = steps_df
    LAST_META["max_disp_px"] = float(max_disp_px)
    LAST_META["fast_disp_px"] = float(fast_disp_px)
    LAST_META["fast_penalty_factor"] = float(fast_penalty_factor)
    LAST_META["area_ratio_max"] = float(area_ratio_max)
    LAST_META["intensity_ratio_max"] = float(intensity_ratio_max)

    _update_tracks_preview(
        viewer,
        tracks_points,
        props_rows,
        color_by=color_by,
        tail_len=int(preview_tail),
        layer_name="mito_tracks_preview",
    )

    n_steps = len(LAST_TRACK_STEPS_DF)
    n_tracks = LAST_TRACK_STEPS_DF["track_id"].nunique()
    n_fast = int(LAST_TRACK_STEPS_DF.get("fast_jump", pd.Series([], dtype=bool)).sum())
    msg = f"LAP tracking finished. Tracks: {n_tracks}, steps: {n_steps}, fast jumps: {n_fast}."
    if "speed_um_per_s" in LAST_TRACK_STEPS_DF.columns:
        msg += " (speed in µm/s computed)"
    show_info(msg)


# ============================================
# Stage 3 - Flow motion filaments (advanced)
# ============================================

def _flow_motion_filament_core(data: np.ndarray,
                               channel: int,
                               t0: int,
                               t1: int,
                               downscale: int,
                               thr_percentile: float,
                               min_area_px: int) -> np.ndarray:
    """
    Flow-based filament map:

      - Dense optical flow between consecutive frames
      - Accumulate motion magnitude
      - Smooth, threshold, skeletonize
    """
    if data.ndim == 3:
        T, Y, X = data.shape

        def get_frame(tt):
            return data[tt]
    else:
        T, C, Y, X = data.shape
        if channel >= C:
            channel = 0

        def get_frame(tt):
            return data[tt, channel]

    ds = max(1, int(downscale))

    frame0 = get_frame(t0)
    if frame0.ndim == 3:
        frame0 = np.max(frame0, axis=0)
    img0 = frame0[::ds, ::ds]
    H, W = img0.shape

    motion_accum = np.zeros((H, W), dtype=np.float32)

    for t in progress(range(t0, t1), desc="Flow motion accumulation"):
        frame0 = get_frame(t)
        frame1 = get_frame(t + 1)

        if frame0.ndim == 3:
            frame0 = np.max(frame0, axis=0)
        if frame1.ndim == 3:
            frame1 = np.max(frame1, axis=0)

        img0 = frame0[::ds, ::ds]
        img1 = frame1[::ds, ::ds]

        h = min(H, img0.shape[0], img1.shape[0])
        w = min(W, img0.shape[1], img1.shape[1])
        img0 = img0[:h, :w]
        img1 = img1[:h, :w]
        if motion_accum.shape != (h, w):
            motion_accum = motion_accum[:h, :w]

        i0 = _normalize_01(img0)
        i1 = _normalize_01(img1)
        prev = (i0 * 255).astype(np.uint8)
        nxt = (i1 * 255).astype(np.uint8)

        flow = cv2.calcOpticalFlowFarneback(
            prev, nxt, None,
            pyr_scale=0.5,
            levels=3,
            winsize=15,
            iterations=3,
            poly_n=5,
            poly_sigma=1.2,
            flags=0,
        )
        fx = flow[..., 0]
        fy = flow[..., 1]
        mag = np.sqrt(fx ** 2 + fy ** 2)

        motion_accum += mag.astype(np.float32)

    max_val = float(motion_accum.max())
    if max_val <= 0:
        raise RuntimeError("Motion accumulation is zero. Check input data or T range.")

    motion_smooth = ndi.gaussian_filter(motion_accum, sigma=2.0)

    thr_hi = np.percentile(motion_smooth, float(thr_percentile))
    thr_lo_percentile = max(40.0, float(thr_percentile) - 35.0)
    thr_lo = np.percentile(motion_smooth, thr_lo_percentile)

    mask_lo = motion_smooth >= thr_lo

    mask = ndi.binary_dilation(mask_lo, iterations=2)
    mask = remove_small_objects(mask, min_size=int(max(1, min_area_px)))
    mask = closing(mask, footprint=disk(2))

    skel = skeletonize(mask.astype(bool))
    labels_skel = label(skel, connectivity=2).astype(np.int32)

    return labels_skel


@magicgui(
    call_button="Flow motion filaments",
    channel={"label": "Channel", "min": 0, "max": 63, "step": 1, "value": 0},
    t_start={"label": "T start (-1=auto)", "min": -1, "max": 999999, "step": 1, "value": -1},
    t_end={"label": "T end (-1=auto)", "min": -1, "max": 999999, "step": 1, "value": -1},
    downscale={"label": "Downscale", "choices": [1, 2, 4, 8], "value": 2},
    thr_percentile={"label": "Motion percentile", "min": 70.0, "max": 99.9, "step": 0.5, "value": 95.0},
    min_area_px={"label": "Min branch area (px)", "min": 1, "max": 50000, "step": 1, "value": 10},
)
def flow_motion_filaments(viewer: napari.Viewer,
                          channel: int = 0,
                          t_start: int = -1,
                          t_end: int = -1,
                          downscale: int = 2,
                          thr_percentile: float = 95.0,
                          min_area_px: int = 10):
    """
    Advanced: build motion-based filament skeleton from long-term optical flow pattern.

    Output layer:
      - motion_filament_skel (labels)
    """
    global LAST_META

    base = _get_layer_by_name(viewer, "CZI_movie")
    if base is None or getattr(base, "data", None) is None:
        show_error("No base image 'CZI_movie' found.")
        return
    data = np.asarray(base.data)

    T = data.shape[0]
    if T < 2:
        show_error("Need at least 2 frames for flow motion filaments.")
        return

    if t_start < 0:
        t0 = 0
    else:
        t0 = max(0, min(int(t_start), T - 2))
    if t_end < 0:
        t1 = T - 1
    else:
        t1 = max(1, min(int(t_end), T - 1))

    if t1 <= t0:
        show_error("Flow motion: need t_end > t_start and at least 2 frames.")
        return

    try:
        labels_skel = _flow_motion_filament_core(
            data=data,
            channel=channel,
            t0=t0,
            t1=t1,
            downscale=downscale,
            thr_percentile=thr_percentile,
            min_area_px=min_area_px,
        )
    except Exception as e:
        show_error(f"Flow motion filaments error: {e}")
        return

    lyr = _get_or_create_labels_layer(viewer, "motion_filament_skel", labels_skel.shape, opacity=0.7)
    lyr.data = labels_skel

    n_branches = int(labels_skel.max())

    LAST_META["flow_downscale"] = int(downscale)
    LAST_META["flow_t_start"] = int(t0)
    LAST_META["flow_t_end"] = int(t1)

    show_info(
        f"Flow motion filaments built.\n"
        f"T range: {t0} to {t1}\n"
        f"Downscale: {downscale}, percentile: {thr_percentile}\n"
        f"Branches: {n_branches}"
    )


# ============================================
# Stage 3b - Fast-track inside skeleton ROI
# ============================================

@magicgui(
    call_button="Fast track in skeleton ROI",
    channel={"label": "Channel", "min": 0, "max": 63, "step": 1, "value": 0},
    t_start={"label": "T start (-1=use flow T)", "min": -1, "max": 999999, "step": 1, "value": -1},
    t_end={"label": "T end (-1=use flow T)", "min": -1, "max": 999999, "step": 1, "value": -1},
    roi_dilate_px={"label": "ROI dilate (px, flow scale)", "min": 0, "max": 20, "step": 1, "value": 2},

    thr_mode={"label": "Thr mode",
              "choices": ["otsu", "mean+std", "percentile", "local"],
              "value": "percentile"},
    thr_corr={"label": "Thr correction", "min": 0.1, "max": 3.0, "step": 0.05, "value": 1.0},
    thr_k={"label": "k (mean+std)", "min": 0.0, "max": 5.0, "step": 0.1, "value": 1.0},
    thr_percentile={"label": "Percentile", "min": 50.0, "max": 99.9, "step": 0.1, "value": 85.0},
    local_block_size={"label": "Local block size", "min": 3, "max": 201, "step": 2, "value": 51},
    local_offset={"label": "Local offset", "min": -0.5, "max": 0.5, "step": 0.01, "value": 0.0},

    bg_sigma={"label": "BG sigma", "min": 0.0, "max": 50.0, "step": 0.5, "value": 4.0},
    clahe_clip={"label": "CLAHE clip (0=off)", "min": 0.0, "max": 0.1, "step": 0.005, "value": 0.02},
    smooth_sigma={"label": "Gaussian sigma", "min": 0.0, "max": 5.0, "step": 0.2, "value": 0.8},

    min_area_px={"label": "Min area (px)", "min": 1, "max": 50000, "step": 1, "value": 10},
    max_area_px={"label": "Max area (px)", "min": 10, "max": 1000000, "step": 10, "value": 4000},
    fill_holes_flag={"label": "Fill holes", "widget_type": "CheckBox", "value": True},
    opening_radius={"label": "Opening radius (px)", "min": 0, "max": 10, "step": 1, "value": 0},
    declump_watershed_flag={"label": "Declump (watershed)", "widget_type": "CheckBox", "value": True},
    peak_min_distance={"label": "Peak min distance", "min": 1, "max": 20, "step": 1, "value": 2},

    max_disp_px={"label": "Max disp (px/frame)", "min": 0.0, "max": 400.0, "step": 0.5, "value": 10.0},
    fast_disp_px={"label": "Fast disp (px/frame)", "min": 0.0, "max": 800.0, "step": 0.5, "value": 20.0},
    fast_penalty_factor={"label": "Fast penalty factor", "min": 1.0, "max": 50.0, "step": 0.5, "value": 2.0},
    area_ratio_max={"label": "Max area ratio", "min": 1.0, "max": 10.0, "step": 0.1, "value": 4.0},
    intensity_ratio_max={"label": "Max intensity ratio", "min": 1.0, "max": 10.0, "step": 0.1, "value": 4.0},

    color_by={"label": "Color by", "choices": ["track_id", "speed"], "value": "speed"},
    preview_tail={"label": "Preview tail length", "min": 0, "max": 500, "step": 10, "value": 120},
)
def fast_tracks_skeleton(viewer: napari.Viewer,
                         channel: int = 0,
                         t_start: int = -1,
                         t_end: int = -1,
                         roi_dilate_px: int = 2,
                         thr_mode: str = "percentile",
                         thr_corr: float = 1.0,
                         thr_k: float = 1.0,
                         thr_percentile: float = 85.0,
                         local_block_size: int = 51,
                         local_offset: float = 0.0,
                         bg_sigma: float = 4.0,
                         clahe_clip: float = 0.02,
                         smooth_sigma: float = 0.8,
                         min_area_px: int = 10,
                         max_area_px: int = 4000,
                         fill_holes_flag: bool = True,
                         opening_radius: int = 0,
                         declump_watershed_flag: bool = True,
                         peak_min_distance: int = 2,
                         max_disp_px: float = 10.0,
                         fast_disp_px: float = 20.0,
                         fast_penalty_factor: float = 2.0,
                         area_ratio_max: float = 4.0,
                         intensity_ratio_max: float = 4.0,
                         color_by: str = "speed",
                         preview_tail: int = 120):
    """
    Advanced fast-track pipeline inside skeleton ROI.

    Requires motion_filament_skel first.
    """
    global FAST_OBJECTS_DF, FAST_TRACK_STEPS_DF, LAST_META

    base = _get_layer_by_name(viewer, "CZI_movie")
    if base is None or getattr(base, "data", None) is None:
        show_error("No base image 'CZI_movie' found.")
        return
    data = np.asarray(base.data)

    skel_layer = _get_layer_by_name(viewer, "motion_filament_skel")
    if skel_layer is None:
        show_error("No skeleton layer 'motion_filament_skel'. Please run Flow motion filaments first.")
        return

    skel = np.asarray(skel_layer.data)
    roi = skel > 0
    if roi_dilate_px and roi_dilate_px > 0:
        roi = ndi.binary_dilation(roi, iterations=int(roi_dilate_px), structure=disk(1))

    Hf, Wf = roi.shape

    ds_flow = int(LAST_META.get("flow_downscale", 1))
    if ds_flow <= 0:
        ds_flow = 1

    T = data.shape[0]
    flow_t0 = int(LAST_META.get("flow_t_start", 0))
    flow_t1 = int(LAST_META.get("flow_t_end", T - 1))

    if t_start < 0:
        t0 = flow_t0
    else:
        t0 = max(0, min(int(t_start), T - 1))
    if t_end < 0:
        t1 = flow_t1
    else:
        t1 = max(0, min(int(t_end), T - 1))

    if t1 <= t0:
        show_error("Fast tracks: need t_end > t_start.")
        return

    if data.ndim == 3:
        def get_frame(tt):
            return data[tt]
    else:
        T_, C, Y, X = data.shape
        if channel >= C:
            channel = 0

        def get_frame(tt):
            return data[tt, channel]

    fast_rows: List[Dict] = []

    for t in progress(range(t0, t1 + 1), desc="Fast detect inside skeleton ROI"):
        frame = get_frame(t)
        if frame.ndim == 3:
            frame = np.max(frame, axis=0)

        img_full = frame[::ds_flow, ::ds_flow]
        img = img_full[:Hf, :Wf]

        labels = _identify_primary_objects(
            img,
            thr_mode=thr_mode,
            thr_corr=thr_corr,
            thr_k=thr_k,
            thr_percentile=thr_percentile,
            smooth_sigma=smooth_sigma,
            bg_sigma=bg_sigma,
            clahe_clip=clahe_clip,
            local_block_size=local_block_size,
            local_offset=local_offset,
            min_area_px=min_area_px,
            max_area_px=max_area_px,
            fill_holes_flag=fill_holes_flag,
            opening_radius=opening_radius,
            declump_watershed_flag=declump_watershed_flag,
            peak_min_distance=peak_min_distance,
            filament_mode=False,
            close_radius=0,
        )

        labels = labels * roi.astype(labels.dtype)

        props_df = _centroids_from_labels(labels, intensity_img=img)
        if props_df.empty:
            continue
        props_df = props_df.reset_index(drop=True)
        for _, r in props_df.iterrows():
            fast_rows.append(
                {
                    "t": int(t),
                    "obj_local": int(r["label"]),
                    "cy": float(r["cy"]),
                    "cx": float(r["cx"]),
                    "area_px": float(r["area"]),
                    "major_axis": float(r["major_axis_length"]),
                    "minor_axis": float(r["minor_axis_length"]),
                    "eccentricity": float(r["eccentricity"]),
                    "perimeter": float(r["perimeter"]),
                    "mean_intensity": float(r["mean_intensity"]),
                }
            )

    if not fast_rows:
        show_error("Fast tracks: no objects detected inside skeleton ROI. Try relaxing thresholds.")
        FAST_OBJECTS_DF = None
        FAST_TRACK_STEPS_DF = None
        return

    FAST_OBJECTS_DF = pd.DataFrame(fast_rows)
    FAST_OBJECTS_DF = FAST_OBJECTS_DF.sort_values(["t", "cy", "cx"]).reset_index(drop=True)
    FAST_OBJECTS_DF["obj_id"] = FAST_OBJECTS_DF.index.astype(int)

    pixel_size_um = LAST_META.get("pixel_size_um", None)
    frame_dt_s = LAST_META.get("frame_dt_s", None)
    if pixel_size_um is None or frame_dt_s is None or frame_dt_s <= 0:
        show_error("Pixel size or frame interval not set. Please set them in calibration panel first.")
        FAST_TRACK_STEPS_DF = None
        return

    p_eff_fast = float(pixel_size_um) * float(ds_flow)

    try:
        steps_df, tracks_points, props_rows = _build_tracks_by_LAP(
            FAST_OBJECTS_DF,
            max_disp_px=max_disp_px,
            pixel_size_um=p_eff_fast,
            frame_dt_s=frame_dt_s,
            area_ratio_max=area_ratio_max,
            intensity_ratio_max=intensity_ratio_max,
            fast_disp_px=fast_disp_px,
            fast_penalty_factor=fast_penalty_factor,
        )
    except Exception as e:
        show_error(f"Fast LAP tracking error: {e}")
        FAST_TRACK_STEPS_DF = None
        return

    steps_df["track_type"] = "FAST_skelROI"
    FAST_TRACK_STEPS_DF = steps_df

    _update_tracks_preview(
        viewer,
        tracks_points,
        props_rows,
        color_by=color_by,
        tail_len=int(preview_tail),
        layer_name="mito_fast_tracks_preview",
    )

    n_steps = len(FAST_TRACK_STEPS_DF)
    n_tracks = FAST_TRACK_STEPS_DF["track_id"].nunique()
    msg = f"Fast tracks finished. Tracks: {n_tracks}, steps: {n_steps}."
    if "speed_um_per_s" in FAST_TRACK_STEPS_DF.columns:
        msg += " (speed_um_per_s available)"
    show_info(msg)


# ============================================
# process_export panel: export detection / tracks
# ============================================

@magicgui(
    call_button="Export CSV / snapshot",
    out_prefix={"label": "Output prefix", "value": "mito"},
    save_snapshot={"label": "Save track snapshot (PNG)", "widget_type": "CheckBox", "value": True},
    export_steps={"label": "Export global LAP per-step", "widget_type": "CheckBox", "value": True},
    export_fast={"label": "Export fast tracks (skeleton ROI)", "widget_type": "CheckBox", "value": True},
)
def export_tracks(viewer: napari.Viewer,
                  out_prefix: str = "mito",
                  save_snapshot: bool = True,
                  export_steps: bool = True,
                  export_fast: bool = True):
    """
    process_export panel (README step):

      - Export global per-object detection table
      - Export LAP track tables (per-track, per-step)
      - Optional fast-track tables and snapshot PNG
    """
    global LAST_OBJECTS_DF, LAST_TRACK_STEPS_DF, FAST_TRACK_STEPS_DF, LAST_META

    has_det = LAST_OBJECTS_DF is not None and not LAST_OBJECTS_DF.empty
    has_lap = LAST_TRACK_STEPS_DF is not None and not LAST_TRACK_STEPS_DF.empty
    has_fast = FAST_TRACK_STEPS_DF is not None and not FAST_TRACK_STEPS_DF.empty

    if not has_det and not has_lap and not has_fast:
        show_error("No detection or tracking results to export.")
        return

    if not has_det:
        show_error("No global detection result to export. Please run detection first.")
        return

    if not has_lap:
        show_error("No global LAP tracking results to export. Please run LAP tracking first.")
        return

    frame_dt_s = LAST_META.get("frame_dt_s", None)
    pixel_size_um = LAST_META.get("pixel_size_um", None)

    if pixel_size_um is None or frame_dt_s is None or frame_dt_s <= 0:
        show_error(
            "Pixel size or frame interval not set.\n"
            "To guarantee speed_um_per_s in all track files, export is blocked.\n"
            "Please set them in calibration panel and click Apply."
        )
        return

    out_dir = QFileDialog.getExistingDirectory(None, "Select output folder", "")
    if not out_dir:
        show_info("Canceled.")
        return

    out_dir = Path(out_dir)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    msg_lines = []

    p_obj = out_dir / f"{out_prefix}_objects_{ts}.csv"
    LAST_OBJECTS_DF.to_csv(p_obj, index=False)
    msg_lines.append(f"Exported global per-object CSV:\n{p_obj}")

    ds_global = LAST_META.get("downscale", 1)
    if ds_global is None or ds_global < 1:
        ds_global = 1

    steps = LAST_TRACK_STEPS_DF.copy()
    grp = steps.groupby("track_id", as_index=False)

    summary = grp.agg(
        t_start=("t0", "min"),
        t_end=("t1", "max"),
        n_steps=("track_id", "size"),
        total_disp_px=("disp_px", "sum"),
        max_step_disp_px=("disp_px", "max"),
    )

    summary["duration_frames"] = summary["t_end"] - summary["t_start"] + 1
    summary["duration_s"] = summary["duration_frames"] * float(frame_dt_s)

    p_eff = float(pixel_size_um) * float(ds_global)
    summary["total_disp_um"] = summary["total_disp_px"] * p_eff
    summary["max_step_disp_um"] = summary["max_step_disp_px"] * p_eff

    if "speed_um_per_s" in steps.columns:
        speed_stats = grp["speed_um_per_s"].agg(
            mean_speed_um_s="mean",
            median_speed_um_s="median",
            max_speed_um_s="max",
        )
        summary = summary.merge(speed_stats, on="track_id", how="left")

    p_tracks = out_dir / f"{out_prefix}_LAP_tracks_{ts}.csv"
    summary.to_csv(p_tracks, index=False)
    msg_lines.append(f"Exported global LAP per-track CSV:\n{p_tracks}")

    if export_steps:
        steps_lap = steps.copy()
        steps_lap["track_type"] = steps_lap.get("track_type", "LAP_global")
        p_steps_lap = out_dir / f"{out_prefix}_LAP_steps_{ts}.csv"
        steps_lap.to_csv(p_steps_lap, index=False)
        msg_lines.append(f"Exported LAP per-step CSV:\n{p_steps_lap}")
    else:
        msg_lines.append("Global per-step CSV export disabled by user.")

    if export_fast and has_fast:
        steps_fast = FAST_TRACK_STEPS_DF.copy()
        grp_fast = steps_fast.groupby("track_id", as_index=False)

        summary_fast = grp_fast.agg(
            t_start=("t0", "min"),
            t_end=("t1", "max"),
            n_steps=("track_id", "size"),
            total_disp_px=("disp_px", "sum"),
            max_step_disp_px=("disp_px", "max"),
        )

        summary_fast["duration_frames"] = summary_fast["t_end"] - summary_fast["t_start"] + 1
        summary_fast["duration_s"] = summary_fast["duration_frames"] * float(frame_dt_s)

        ds_flow = int(LAST_META.get("flow_downscale", 1))
        if ds_flow <= 0:
            ds_flow = 1
        p_eff_fast = float(pixel_size_um) * float(ds_flow)

        summary_fast["total_disp_um"] = summary_fast["total_disp_px"] * p_eff_fast
        summary_fast["max_step_disp_um"] = summary_fast["max_step_disp_px"] * p_eff_fast

        if "speed_um_per_s" in steps_fast.columns:
            speed_stats_fast = grp_fast["speed_um_per_s"].agg(
                mean_speed_um_s="mean",
                median_speed_um_s="median",
                max_speed_um_s="max",
            )
            summary_fast = summary_fast.merge(speed_stats_fast, on="track_id", how="left")

        p_fast_tracks = out_dir / f"{out_prefix}_FAST_tracks_{ts}.csv"
        summary_fast.to_csv(p_fast_tracks, index=False)
        msg_lines.append(f"Exported FAST per-track CSV (skeleton ROI):\n{p_fast_tracks}")

        p_fast_steps = out_dir / f"{out_prefix}_FAST_steps_{ts}.csv"
        steps_fast.to_csv(p_fast_steps, index=False)
        msg_lines.append(f"Exported FAST per-step CSV:\n{p_fast_steps}")
    elif export_fast:
        msg_lines.append("No FAST tracks to export.")

    if save_snapshot:
        p_png = out_dir / f"{out_prefix}_tracks_snapshot_{ts}.png"
        try:
            viewer.screenshot(path=p_png, canvas_only=True)
            msg_lines.append(f"Snapshot saved:\n{p_png}")
        except Exception as e:
            show_error(f"Failed to save snapshot PNG: {e}")

    show_info("\n".join(msg_lines))


# ============================================
# GUI parameter save / load (advanced)
# ============================================

GUI_WIDGET_PARAM_MAP: Dict[str, List[str]] = {
    "acquisition_params": [
        "pixel_size_um",
        "frame_dt_s",
    ],
    "identify_preview": [
        "channel",
        "downscale",
        "thr_mode",
        "thr_corr",
        "thr_k",
        "thr_percentile",
        "local_block_size",
        "local_offset",
        "bg_sigma",
        "clahe_clip",
        "smooth_sigma",
        "min_area_px",
        "max_area_px",
        "fill_holes_flag",
        "opening_radius",
        "declump_watershed_flag",
        "peak_min_distance",
        "filament_mode",
        "close_radius",
    ],
    "detect_all_frames": [
        "channel",
        "t_start",
        "t_end",
        "downscale",
        "thr_mode",
        "thr_corr",
        "thr_k",
        "thr_percentile",
        "local_block_size",
        "local_offset",
        "bg_sigma",
        "clahe_clip",
        "smooth_sigma",
        "min_area_px",
        "max_area_px",
        "fill_holes_flag",
        "opening_radius",
        "declump_watershed_flag",
        "peak_min_distance",
        "filament_mode",
        "close_radius",
    ],
    "lap_tracking": [
        "max_disp_px",
        "fast_disp_px",
        "fast_penalty_factor",
        "area_ratio_max",
        "intensity_ratio_max",
        "color_by",
        "preview_tail",
    ],
    "flow_motion_filaments": [
        "channel",
        "t_start",
        "t_end",
        "downscale",
        "thr_percentile",
        "min_area_px",
    ],
    "fast_tracks_skeleton": [
        "channel",
        "t_start",
        "t_end",
        "roi_dilate_px",
        "thr_mode",
        "thr_corr",
        "thr_k",
        "thr_percentile",
        "local_block_size",
        "local_offset",
        "bg_sigma",
        "clahe_clip",
        "smooth_sigma",
        "min_area_px",
        "max_area_px",
        "fill_holes_flag",
        "opening_radius",
        "declump_watershed_flag",
        "peak_min_distance",
        "max_disp_px",
        "fast_disp_px",
        "fast_penalty_factor",
        "area_ratio_max",
        "intensity_ratio_max",
        "color_by",
        "preview_tail",
    ],
    "export_tracks": [
        "out_prefix",
        "save_snapshot",
        "export_steps",
        "export_fast",
    ],
}


def _get_widget_by_name(widget_name: str):
    """
    Return magicgui widget by its function name.
    """
    return globals().get(widget_name, None)


def _collect_gui_params() -> Dict:
    """
    Collect GUI parameters into a dict for JSON saving.
    """
    config: Dict[str, Dict[str, object]] = {}

    for widget_name, param_list in GUI_WIDGET_PARAM_MAP.items():
        w = _get_widget_by_name(widget_name)
        if w is None:
            continue

        widget_cfg: Dict[str, object] = {}
        for p in param_list:
            if not hasattr(w, p):
                continue
            try:
                val = getattr(w, p).value
            except Exception:
                try:
                    val = w[p].value
                except Exception:
                    continue
            widget_cfg[p] = val

        if widget_cfg:
            config[widget_name] = widget_cfg

    meta = {
        "version": "MitoDiffTrack",
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    return {
        "meta": meta,
        "widgets": config,
    }


def _apply_gui_params(config: Dict):
    """
    Apply JSON config back to GUI widgets.
    Note: for Pixel / Time, user can click Apply again to refresh internal state.
    """
    widgets_cfg = config.get("widgets", {})
    if not isinstance(widgets_cfg, dict):
        return

    for widget_name, params in widgets_cfg.items():
        w = _get_widget_by_name(widget_name)
        if w is None or not isinstance(params, dict):
            continue

        for p_name, p_val in params.items():
            if not hasattr(w, p_name):
                continue
            try:
                child = getattr(w, p_name)
                child.value = p_val
            except Exception:
                try:
                    w[p_name].value = p_val
                except Exception:
                    continue


@magicgui(
    call_button="Save GUI config",
    config_label={"label": "Config label", "value": "default"},
)
def save_gui_config(viewer: napari.Viewer,
                    config_label: str = "default"):
    """
    Save GUI parameters to a JSON file.
    """
    try:
        default_name = f"mito_gui_config_{config_label}.json"
        path, _ = QFileDialog.getSaveFileName(
            None,
            "Save GUI config",
            default_name,
            "JSON files (*.json);;All files (*.*)",
        )
    except Exception as e:
        show_error(f"Open save dialog failed: {e}")
        return

    if not path:
        show_info("Save canceled.")
        return

    cfg = _collect_gui_params()
    cfg["meta"]["label"] = config_label

    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=2, ensure_ascii=False)
    except Exception as e:
        show_error(f"Failed to save config: {e}")
        return

    show_info(f"GUI config saved:\n{path}")


@magicgui(
    call_button="Load GUI config",
)
def load_gui_config(viewer: napari.Viewer):
    """
    Load GUI parameters from a JSON file and update widgets.
    """
    try:
        path, _ = QFileDialog.getOpenFileName(
            None,
            "Load GUI config",
            "",
            "JSON files (*.json);;All files (*.*)",
        )
    except Exception as e:
        show_error(f"Open load dialog failed: {e}")
        return

    if not path:
        show_info("Load canceled.")
        return

    try:
        with open(path, "r", encoding="utf-8") as f:
            cfg = json.load(f)
    except Exception as e:
        show_error(f"Failed to read config: {e}")
        return

    if not isinstance(cfg, dict):
        show_error("Invalid config format.")
        return

    _apply_gui_params(cfg)

    meta = cfg.get("meta", {})
    label = meta.get("label", "")
    msg = f"GUI config loaded from:\n{path}"
    if label:
        msg += f"\nLabel: {label}"
    msg += "\nNote: for Pixel / Time, it is recommended to click Apply once after loading."
    show_info(msg)


# ============================================
# batch_process_folder panel (placeholder)
# ============================================

@magicgui(
    call_button="Run batch",
    in_folder={"label": "Input folder", "mode": "d"},
    out_folder={"label": "Output folder", "mode": "d"},
)
def batch_process_folder(viewer: napari.Viewer,
                         in_folder: Path = Path(""),
                         out_folder: Path = Path("")):
    """
    batch_process_folder panel (README step: "Optional batch analysis of multiple movies").

    Currently implemented as a placeholder.
    In this version, please process single files using load_czi / seg_params / process_export / track_mito.
    """
    show_info(
        "Batch processing placeholder.\n"
        "This version does not yet run full batch.\n"
        "Please process single files following the README steps."
    )


# ============================================
# Window and main + plugin entry
# ============================================

def _center_and_resize(viewer, width_ratio=0.9, height_ratio=0.9, min_w=900, min_h=600):
    """Center and resize napari window."""
    app = QtWidgets.QApplication.instance()
    if app is None:
        return
    screen = app.primaryScreen()
    if screen is None:
        return
    geo = screen.availableGeometry()
    w = max(int(geo.width() * width_ratio), min_w)
    h = max(int(geo.height() * height_ratio), min_h)
    win = viewer.window._qt_window
    win.setMinimumSize(min_w, min_h)
    win.resize(w, h)
    center = geo.center()
    win.move(center.x() - w // 2, center.y() - h // 2)


def mitochondria_motility_analyzer(viewer: napari.Viewer):
    """
    Entry point for napari plugin:

      Plugins → Mitochondria Motility Analyzer (python based plugin)

    It will open the main panels described in the README.
    """
    # Main panels, following README order
    viewer.window.add_dock_widget(load_czi,              name="load_czi",              area="right")
    viewer.window.add_dock_widget(acquisition_params,    name="calibration",           area="right")
    viewer.window.add_dock_widget(identify_preview,      name="seg_params",            area="right")
    viewer.window.add_dock_widget(detect_all_frames,     name="process_export (seg)",  area="right")
    viewer.window.add_dock_widget(lap_tracking,          name="track_mito",            area="right")
    viewer.window.add_dock_widget(batch_process_folder,  name="batch_process_folder",  area="right")

    # Advanced panels
    viewer.window.add_dock_widget(flow_motion_filaments, name="flow_motion_filaments", area="right")
    viewer.window.add_dock_widget(fast_tracks_skeleton,  name="fast_tracks_skeleton",  area="right")
    viewer.window.add_dock_widget(export_tracks,         name="process_export (CSV)",  area="right")
    viewer.window.add_dock_widget(save_gui_config,       name="save_gui_config",       area="right")
    viewer.window.add_dock_widget(load_gui_config,       name="load_gui_config",       area="right")

    _center_and_resize(viewer)
    show_info(
        "Mitochondria Motility Analyzer loaded.\n"
        "Main panels:\n"
        " 1) load_czi + calibration\n"
        " 2) seg_params\n"
        " 3) process_export (seg)\n"
        " 4) track_mito\n"
        " 5) batch_process_folder (placeholder)\n"
        "Advanced:\n"
        " - flow_motion_filaments, fast_tracks_skeleton, process_export (CSV), GUI config\n"
    )
    return viewer


def main():
    viewer = napari.Viewer()
    mitochondria_motility_analyzer(viewer)
    napari.run()


if __name__ == "__main__":
    main()
