"""
Unified Cellpose-based segmentation pipeline for single- or dual-population images.

The script can:
  - segment nuclei alone, or segment cells plus nuclei when a cell/cytoplasmic
    signal is available;
  - work with a single nuclear channel or with two population-specific nuclear
    channels;
  - optionally quantify a protein/intensity channel in the nucleus, whole cell,
    cytoplasm, Voronoi territory, or pseudo-cytoplasm;
  - optionally generate Voronoi territories from nuclear centroids as a
    pseudo-cell/pseudo-cytoplasm compartment when no actin/cell signal is
    available;
  - export label masks, QC overlays and per-object measurements.

Expected image layout:
  - single-plane multichannel images: (C, Y, X)
  - z-stack multichannel images: (Z, C, Y, X)
For z-stacks, each requested channel is extracted along axis 1 and projected along Z axis 0.
Channel indices in the YAML always refer to the C axis.

This version uses Cellpose 4 / Cellpose-SAM via models.CellposeModel. The legacy
Cellpose 3 models.Cellpose class and eval(channels=...) API are not used.

Run:
  python scripts/unified_segmentation.py path/to/config.yml
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yaml
from scipy import ndimage as ndi
from skimage import filters, measure, morphology
from skimage.segmentation import watershed
from tifffile import imread, imwrite
from tqdm import tqdm

import warnings

warnings.filterwarnings(
    "ignore",
    message=".*Sparse invariant checks are implicitly disabled.*",
    category=UserWarning,
)

mpl.rcParams["figure.dpi"] = 300

Mask = np.ndarray
Image = np.ndarray
MISSING = ""

def require_key(mapping: Mapping, key: str):
    """Return a required configuration value with a clear error message."""
    if key not in mapping or mapping[key] is None:
        raise ValueError(f"Missing required configuration key: {key}")
    return mapping[key]


def get_param(params: Mapping, key: str, default=None):
    """Read a configuration value with defaults."""
    return params.get(key, default)


def get_nuclear_mode(params: Mapping) -> str:
    """Return a normalized nuclear mode.

    The legacy value "mixed" is accepted for backward compatibility but is
    normalized to "dual_population". New configuration files should use
    "dual_population".
    """
    nuclear_mode = get_param(params, "nuclear_mode", "single")
    if nuclear_mode == "mixed":
        print(
            "Warning: nuclear_mode: 'mixed' is deprecated. "
            "Use nuclear_mode: 'dual_population' instead.",
            flush=True,
        )
        return "dual_population"
    if nuclear_mode not in {"single", "dual_population"}:
        raise ValueError(
            "nuclear_mode must be either 'single' or 'dual_population'. "
            "The legacy alias 'mixed' is deprecated."
        )
    return nuclear_mode

def get_cellpose_device(params: Mapping):

    import torch

    use_gpu = bool(get_param(params, "use_gpu", False))
    requested = str(get_param(params, "gpu_device", "auto")).lower()

    if requested not in {"auto", "cuda", "mps", "cpu"}:
        raise ValueError("gpu_device must be one of: 'auto', 'cuda', 'mps', 'cpu'.")

    if not use_gpu or requested == "cpu":
        return torch.device("cpu")

    if requested == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError(
                "gpu_device: 'cuda' was requested, but CUDA is not available. "
                "Set use_gpu: false, gpu_device: 'cpu', or install a CUDA-compatible PyTorch build."
            )
        return torch.device("cuda")

    if requested == "mps":
        if not hasattr(torch.backends, "mps") or not torch.backends.mps.is_available():
            raise RuntimeError(
                "gpu_device: 'mps' was requested, but PyTorch MPS is not available. "
                "Set use_gpu: false, gpu_device: 'cpu', or install a PyTorch build with MPS support."
            )
        return torch.device("mps")

    # auto mode
    if torch.cuda.is_available():
        return torch.device("cuda")

    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")

    return torch.device("cpu")

def load_raw_image(filename: str | Path) -> Image:
    """Load a TIFF as a NumPy array."""
    return np.asarray(imread(str(filename)))

def project_image_channel(raw: Image, channel_index: int, projection: str = "max") -> Image:
    """Extract one 2D channel using the pipeline's fixed dimension convention.

    Expected input formats:
      - 3D single-plane multichannel image: (C, Y, X)
      - 4D z-stack multichannel image: (Z, C, Y, X)

    For 4D images, the requested channel is extracted along axis 1 and then
    projected along Z axis 0.
    """
    arr = np.asarray(raw)

    if arr.ndim == 2:
        if channel_index not in (0, None):
            raise IndexError(
                "A 2D image has no channel axis; use channel index 0 or omit the channel."
            )
        return arr.astype(np.float32, copy=False)

    if arr.ndim == 3:
        # Expected layout: (C, Y, X)
        if channel_index < 0 or channel_index >= arr.shape[0]:
            raise IndexError(
                f"Channel index {channel_index} is out of bounds for image with "
                f"{arr.shape[0]} channels. Expected image layout: (C, Y, X)."
            )
        return arr[int(channel_index)].astype(np.float32, copy=False)

    if arr.ndim == 4:
        if channel_index < 0 or channel_index >= arr.shape[1]:
            raise IndexError(
                f"Channel index {channel_index} is out of bounds for image with "
                f"{arr.shape[1]} channels. Expected image layout: (Z, C, Y, X)."
            )

        channel_stack = arr[:, int(channel_index), :, :]

        if projection == "max":
            channel = np.max(channel_stack, axis=0)
        elif projection == "mean":
            channel = np.mean(channel_stack, axis=0)
        elif projection == "sum":
            channel = np.sum(channel_stack, axis=0)
        else:
            raise ValueError(
                f"Unsupported projection: '{projection}'. Use 'max', 'mean', or 'sum'."
            )

        return channel.astype(np.float32, copy=False)

    raise ValueError(
        f"Unsupported image shape {arr.shape}. "
        "Expected a 2D image, a 3D image with shape (C, Y, X), "
        "or a 4D z-stack with shape (Z, C, Y, X)."
    )


def extract_channels(raw: Image, params: Mapping) -> Tuple[Optional[Image], Dict[str, Image], Optional[Image]]:
    """Extract optional cell, nuclear and optional protein channels.

    Returned nuclear keys are output/analysis labels rather than channel names:
      - single mode: ``all_cells``
      - dual-population mode: ``all_cells``, ``population_1``, ``population_2``

    The YAML still uses explicit channel names such as ``nucleus``,
    ``population_1_nucleus`` and ``population_2_nucleus`` because these describe
    the input channels.
    """
    channels = require_key(params, "channels")
    nuclear_mode = get_nuclear_mode(params)
    segmentation_mode = get_param(params, "segmentation_mode", "cell_nucleus")
    quantify_protein = bool(get_param(params, "quantify_protein", False))
    default_projection = get_param(params, "z_projection", "max")
    protein_projection = get_param(params, "protein_z_projection", default_projection)

    print(f"Raw image shape: {raw.shape}", flush=True)
    if raw.ndim == 3:
        print("Interpreting image as (C, Y, X).", flush=True)
    elif raw.ndim == 4:
        print("Interpreting image as (Z, C, Y, X).", flush=True)

    def channel(name: str, *, required: bool = True, projection: Optional[str] = None) -> Optional[Image]:
        if name not in channels or channels[name] is None:
            if required:
                raise ValueError(f"Missing channel mapping for '{name}'.")
            return None
        return project_image_channel(raw, int(channels[name]), projection or default_projection)
    

    cell = None
    if segmentation_mode == "cell_nucleus":
        cell = channel("cell", required=True)
    elif segmentation_mode == "nucleus_only":
        cell = channel("cell", required=False)
    else:
        raise ValueError("segmentation_mode must be 'cell_nucleus' or 'nucleus_only'.")

    if nuclear_mode == "single":
        nucleus_key = "nucleus" if "nucleus" in channels else "nuclei"
        nuclei = {"all_cells": channel(nucleus_key, required=True)}
    elif nuclear_mode == "dual_population":
        population_1_key = "population_1_nucleus" if "population_1_nucleus" in channels else "red_nuclei"
        population_2_key = "population_2_nucleus" if "population_2_nucleus" in channels else "green_nuclei"
        population_1 = channel(population_1_key, required=True)
        population_2 = channel(population_2_key, required=True)
        if "total_nuclei" in channels and channels["total_nuclei"] is not None:
            total = channel("total_nuclei", required=True)
        else:
            total = population_1 + population_2
        nuclei = {
            "all_cells": total,
            "population_1": population_1,
            "population_2": population_2,
        }
    else:
        raise ValueError("nuclear_mode must be either 'single' or 'dual_population'.")

    protein = channel("protein", required=True, projection=protein_projection) if quantify_protein else None
    return cell, nuclei, protein

def preprocess_channel(raw: Image, gamma: float, sigma: float) -> Image:
    """Apply gamma correction and optional Gaussian filtering."""
    if gamma <= 0:
        raise ValueError(f"gamma must be > 0, got {gamma}.")
    raw = raw.astype(np.float32, copy=False)
    raw = np.clip(raw, a_min=0, a_max=None)
    processed = np.power(raw, gamma)
    if sigma and sigma > 0:
        processed = filters.gaussian(processed, sigma=sigma, preserve_range=True)
    return processed.astype(np.float32, copy=False)

def initialize_cellpose_models(params: Mapping) -> Tuple[Optional[object], object]:

    from cellpose import models

    device = get_cellpose_device(params)
    use_gpu = device.type != "cpu"

    print(f"Using Cellpose device: {device}", flush=True)

    nucleus_model = models.CellposeModel(
        gpu=use_gpu,
        device=device,
    )

    cell_model = None
    if get_param(params, "segmentation_mode", "cell_nucleus") == "cell_nucleus":
        cell_model = models.CellposeModel(
            gpu=use_gpu,
            device=device,
        )

    return cell_model, nucleus_model

def run_cellpose_cell_segmentation(cell_raw: Image, nuclei_for_cell: Optional[Image], params: Mapping, model: object) -> Tuple[Mask, object]:
    """Segment cell/cytoplasmic signal using Cellpose 4 / Cellpose-SAM."""
    cell_pre = preprocess_channel(
        cell_raw,
        gamma=float(get_param(params, "gamma_cell", 1.0)),
        sigma=float(get_param(params, "gaussian_cell", 0.0)),
    )
    if bool(get_param(params, "nu_in_cell_seg", False)) and nuclei_for_cell is not None:
        nuclei_pre = filters.gaussian(
            nuclei_for_cell.astype(np.float32, copy=False),
            sigma=float(get_param(params, "gaussian_cell", 0.0)),
            preserve_range=True,
        )
        image_for_cellpose = np.stack((cell_pre, nuclei_pre), axis=0)
        print("Segmenting cells with cytoplasmic + nuclear signal...", flush=True)
    else:
        image_for_cellpose = cell_pre
        print("Segmenting cells with cytoplasmic signal only...", flush=True)

    result = model.eval(
        image_for_cellpose,
        diameter=get_param(params, "estimated_cell_diameter", None),
        flow_threshold=float(get_param(params, "flow_threshold", 0.4)),
        cellprob_threshold=float(get_param(params, "cellprob_threshold", 0.0)),
    )
    masks, flows = result[0], result[1]
    return masks.astype(np.uint16, copy=False), flows


def run_cellpose_nuclear_segmentation(nuclei_raw: Image, params: Mapping, model: object) -> Tuple[Mask, object]:
    """Segment nuclear signal using Cellpose 4 / Cellpose-SAM."""
    nuclei_pre = preprocess_channel(
        nuclei_raw,
        gamma=float(get_param(params, "gamma_nu", 1.0)),
        sigma=float(get_param(params, "gaussian_nu", 0.0)),
    )

    result = model.eval(
        nuclei_pre,
        diameter=get_param(params, "estimated_nucleus_diameter", None),
        flow_threshold=float(get_param(params, "flow_threshold", 0.4)),
        cellprob_threshold=float(get_param(params, "cellprob_threshold", 0.0)),
    )
    masks, flows = result[0], result[1]
    return masks.astype(np.uint16, copy=False), flows


def filter_by_area(mask: Mask, intensity_image: Image, min_area: float, object_name: str = "objects", report_lines: Optional[List[str]] = None) -> Mask:
    out = mask.copy()
    props = measure.regionprops(mask, intensity_image)
    before = len(props)
    for prop in props:
        if prop.area < min_area:
            out[out == prop.label] = 0
    after = len(measure.regionprops(out, intensity_image))
    message = (
        f"Filtered {before - after} {object_name} by area "
        f"({before} -> {after}; min_area={min_area})."
    )
    print(message, flush=True)
    if report_lines is not None:
        report_lines.append(message)

    return out.astype(np.uint16, copy=False)


def filter_nuclei(mask: Mask, nuclei_raw: Image, params: Mapping, object_name: str = "nuclei", report_lines: Optional[List[str]] = None) -> Mask:
    out = mask.copy()
    min_area = float(get_param(params, "filter_area_nu", 0))
    max_ecc = float(get_param(params, "filter_ecc", 1.0))
    props = measure.regionprops(mask, nuclei_raw)
    before = len(props)
    removed_area = 0
    removed_eccentricity = 0
    for prop in props:
        if prop.area < min_area:
            out[out == prop.label] = 0
            removed_area += 1
        elif prop.eccentricity > max_ecc:
            out[out == prop.label] = 0
            removed_eccentricity += 1
    after = len(measure.regionprops(out, nuclei_raw))
    message = (
        f"Filtered {before - after} {object_name} "
        f"({before} -> {after}; area={removed_area}, eccentricity={removed_eccentricity})."
    )
    print(message, flush=True)
    if report_lines is not None:
        report_lines.append(message)
    return out.astype(np.uint16, copy=False)


def make_foreground_mask(image: Image, params: Mapping) -> np.ndarray:
    """Create a permissive foreground mask, mostly to limit Voronoi territories."""
    sigma = float(get_param(params, "foreground_sigma", 2.0))
    percentile = float(get_param(params, "foreground_percentile", 5.0))
    min_size = int(get_param(params, "foreground_min_size", 0))
    smooth = filters.gaussian(image.astype(np.float32), sigma=sigma, preserve_range=True)
    mask = smooth > np.percentile(smooth, percentile)
    if min_size > 0:
        mask = morphology.remove_small_objects(mask, min_size=min_size)
    return mask.astype(bool)


def voronoi_from_centroids(nuclei_labels: Mask, mask: Optional[np.ndarray] = None) -> Mask:

    nuclei_labels = nuclei_labels.astype(np.int32, copy=False)
    props = measure.regionprops(nuclei_labels)
    if not props:
        return np.zeros_like(nuclei_labels, dtype=np.uint16)

    markers = np.zeros_like(nuclei_labels, dtype=np.int32)
    for prop in props:
        row, col = prop.centroid
        rr = int(np.clip(round(row), 0, markers.shape[0] - 1))
        cc = int(np.clip(round(col), 0, markers.shape[1] - 1))
        markers[rr, cc] = int(prop.label)

    distance_to_markers = ndi.distance_transform_edt(markers == 0)
    if mask is None:
        mask = np.ones_like(nuclei_labels, dtype=bool)
    territories = watershed(distance_to_markers, markers=markers, mask=mask)
    return territories.astype(np.uint16, copy=False)

def voronoi_from_nuclear_masks(nuclei_labels: Mask, mask: Optional[np.ndarray] = None) -> Mask:

    nuclei_labels = nuclei_labels.astype(np.int32, copy=False)

    if nuclei_labels.max() == 0:
        return np.zeros_like(nuclei_labels, dtype=np.uint16)

    if mask is None:
        mask = np.ones_like(nuclei_labels, dtype=bool)
    else:
        mask = mask.astype(bool)

    # Make sure nuclear seed regions are included in the foreground.
    mask = np.logical_or(mask, nuclei_labels > 0)

    # Expand labelled nuclei into the foreground mask.
    distance = ndi.distance_transform_edt(nuclei_labels == 0)
    territories = watershed(distance, markers=nuclei_labels, mask=mask)

    return territories.astype(np.uint16, copy=False)

def make_voronoi_territories(nuclei_labels: Mask, mask: Optional[np.ndarray], params: Mapping) -> Mask:

    seed_mode = get_param(params, "voronoi_seed_mode", "nuclear_mask")

    if seed_mode == "nuclear_mask":
        return voronoi_from_nuclear_masks(nuclei_labels, mask=mask)

    if seed_mode == "centroid":
        return voronoi_from_centroids(nuclei_labels, mask=mask)

    raise ValueError(
        f"Unsupported voronoi_seed_mode: {seed_mode}. "
        "Use 'nuclear_mask' or 'centroid'."
    )

def make_pseudo_cytoplasm(territory_mask: Mask, nucleus_mask: Mask, params: Mapping) -> Mask:
    """Create a pseudo-cytoplasm mask from Voronoi territories minus nuclei.

    If voronoi_ring_width_px is set to a positive value, the pseudo-cytoplasm is
    restricted to a perinuclear ring within the territory. Otherwise the full
    Voronoi territory outside the nucleus is used.
    """
    pseudo = np.where(nucleus_mask == 0, territory_mask, 0).astype(np.uint16)
    ring_width = int(get_param(params, "voronoi_ring_width_px", 0) or 0)
    if ring_width <= 0:
        return pseudo

    restricted = np.zeros_like(pseudo, dtype=np.uint16)
    for prop in measure.regionprops(nucleus_mask):
        label = int(prop.label)
        nuc = nucleus_mask == label
        territory = territory_mask == label
        dilated = morphology.binary_dilation(nuc, morphology.disk(ring_width))
        restricted[np.logical_and.reduce((territory, dilated, ~nuc))] = label
    return restricted.astype(np.uint16, copy=False)

def associate_nuclei_to_cells(nuclei_mask: Mask, cell_mask: Mask, cell_raw: Image, params: Mapping, object_name: str = "nuclei",container_name: str = "cells/territories", report_lines: Optional[List[str]] = None) -> Tuple[Mask, Mask, Mask]:
    """Filter nuclei by cell overlap and return final cell/nucleus/cytoplasm masks.

    The nuclear mask is converted to cell labels, so the final nuclear label id matches
    the corresponding cell label id.
    """
    min_nucleus_area = float(get_param(params, "filter_area_nu", 0))
    ambiguous_ratio = float(get_param(params, "ambiguous_overlap_ratio", 0.8))

    cell_props = measure.regionprops(cell_mask, cell_raw)
    cell_area = {p.label: p.area for p in cell_props}
    filtered_nuclei = nuclei_mask.copy()

    n_nuclei_before = len(measure.regionprops(nuclei_mask))
    n_containers_before = len(cell_props)

    removed_no_container = 0
    removed_too_many_containers = 0
    removed_ambiguous_overlap = 0
    removed_small_overlap = 0
    removed_overlap_larger_than_container = 0
    ignored_tiny_overlaps = 0

    for prop in measure.regionprops(nuclei_mask):
        nucleus_region = nuclei_mask == prop.label
        overlapping_cells, counts = np.unique(cell_mask[nucleus_region], return_counts=True)
        valid = overlapping_cells != 0
        overlapping_cells = overlapping_cells[valid]
        counts = counts[valid]

        min_overlap_fraction = float(get_param(params, "min_container_overlap_fraction", 0.02))

        if len(counts) > 0:
            nucleus_area = np.sum(nucleus_region)

            keep = counts / nucleus_area >= min_overlap_fraction
            ignored_tiny_overlaps += int(np.sum(~keep))

            overlapping_cells = overlapping_cells[keep]
            counts = counts[keep]

        if len(overlapping_cells) == 0:
            filtered_nuclei[nucleus_region] = 0
            removed_no_container += 1
            continue

        if len(overlapping_cells) > 2:
            filtered_nuclei[nucleus_region] = 0
            removed_too_many_containers += 1
            continue

        if len(overlapping_cells) == 2:
            order = np.argsort(counts)[::-1]
            major_cell = int(overlapping_cells[order[0]])
            minor_cell = int(overlapping_cells[order[1]])
            major_count = int(counts[order[0]])
            minor_count = int(counts[order[1]])

            if min(major_count, minor_count) / max(major_count, minor_count) > ambiguous_ratio:
                filtered_nuclei[nucleus_region] = 0
                removed_ambiguous_overlap += 1
                continue

            if major_count >= cell_area.get(major_cell, np.inf):
                filtered_nuclei[nucleus_region] = 0
                removed_overlap_larger_than_container += 1
                continue

            filtered_nuclei[np.logical_and(nucleus_region, cell_mask != major_cell)] = 0
            continue

        cell_label = int(overlapping_cells[0])
        overlap_area = int(counts[0])

        if overlap_area < min_nucleus_area:
            filtered_nuclei[nucleus_region] = 0
            removed_small_overlap += 1
        elif overlap_area >= cell_area.get(cell_label, np.inf):
            filtered_nuclei[nucleus_region] = 0
            removed_overlap_larger_than_container += 1
        else:
            # Keep only the portion of this nucleus inside the assigned cell.
            # This removes tiny ignored overlaps with adjacent cells.
            filtered_nuclei[np.logical_and(nucleus_region, cell_mask != cell_label)] = 0

    final_cell_mask = cell_mask.copy()
    removed_empty_containers = 0
    removed_multi_nucleus_containers = 0
    removed_nuclei_in_multi_nucleus_containers = 0
    kept_largest_nucleus_in_multi_nucleus_containers = 0

    multi_nucleus_policy = get_param(params, "multi_nucleus_cell_policy", "remove_cell")
    if multi_nucleus_policy not in {"remove_cell", "keep_largest", "keep_all"}:
        raise ValueError(
            "multi_nucleus_cell_policy must be 'remove_cell', 'keep_largest', or 'keep_all'."
        )

    for prop in cell_props:
        container_region = cell_mask == prop.label
        nuclear_labels = np.unique(filtered_nuclei[container_region])
        nuclear_labels = nuclear_labels[nuclear_labels != 0]

        if len(nuclear_labels) == 0:
            final_cell_mask[final_cell_mask == prop.label] = 0
            removed_empty_containers += 1
            continue

        if len(nuclear_labels) > 1:
            if multi_nucleus_policy == "remove_cell":
                final_cell_mask[final_cell_mask == prop.label] = 0
                for nuc_label in nuclear_labels:
                    filtered_nuclei[filtered_nuclei == nuc_label] = 0
                removed_multi_nucleus_containers += 1
                removed_nuclei_in_multi_nucleus_containers += len(nuclear_labels)

            elif multi_nucleus_policy == "keep_largest":
                areas = {
                    int(nuc_label): int(np.sum(filtered_nuclei[container_region] == nuc_label))
                    for nuc_label in nuclear_labels
                }
                keep_label = max(areas, key=areas.get)

                for nuc_label in nuclear_labels:
                    if int(nuc_label) != int(keep_label):
                        filtered_nuclei[filtered_nuclei == nuc_label] = 0

                kept_largest_nucleus_in_multi_nucleus_containers += 1

            elif multi_nucleus_policy == "keep_all":
                pass

    n_nuclei_after = len(measure.regionprops(filtered_nuclei))
    n_containers_after = len(measure.regionprops(final_cell_mask))

    message = (
        f"Association filtering for {object_name} in {container_name}: "
        f"nuclei {n_nuclei_before} -> {n_nuclei_after} "
        f"(removed: no container={removed_no_container}, "
        f">2 containers={removed_too_many_containers}, "
        f"ambiguous overlap={removed_ambiguous_overlap}, "
        f"small overlap={removed_small_overlap}, "
        f"overlap >= container={removed_overlap_larger_than_container}, "
        f"ignored tiny overlaps={ignored_tiny_overlaps}, "
        f"multi-nucleus containers removed={removed_multi_nucleus_containers}, "
        f"nuclei removed in multi-nucleus containers={removed_nuclei_in_multi_nucleus_containers}, "
        f"multi-nucleus containers kept with largest nucleus={kept_largest_nucleus_in_multi_nucleus_containers}); "
        f"{container_name} {n_containers_before} -> {n_containers_after} "
        f"(removed empty={removed_empty_containers})."
    )
    print(message, flush=True)
    if report_lines is not None:
        report_lines.append(message)

    nucleus_with_cell_labels = np.where(filtered_nuclei != 0, final_cell_mask, 0).astype(np.uint16)
    cytoplasm_mask = np.where(nucleus_with_cell_labels == 0, final_cell_mask, 0).astype(np.uint16)
    return final_cell_mask.astype(np.uint16), nucleus_with_cell_labels, cytoplasm_mask


def save_mask(mask: Mask, path: str | Path):
    """Save a label mask as uint16 TIFF."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    imwrite(str(path), mask.astype(np.uint16, copy=False))


def random_label_cmap(max_labels: int = 65536):
    rng = np.random.default_rng(42)
    colors = rng.random((max_labels, 3))
    colors[0] = 0
    return mpl.colors.ListedColormap(colors)


def save_overlay(raw: Image, mask: Mask, out_path: str | Path, title: str):
    """Save a QC overlay of a label mask on a raw image."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmap = random_label_cmap()
    plt.figure(figsize=(6, 6))
    plt.imshow(raw, cmap="gray", vmin=np.percentile(raw, 0.1), vmax=np.percentile(raw, 99.7))
    plt.imshow(cmap(mask), alpha=0.45)
    plt.title(title, fontsize=14)
    plt.axis("off")
    plt.savefig(str(out_path), bbox_inches="tight", pad_inches=0, transparent=True)
    plt.close()

def save_outline_overlay(raw: Image, masks: Mapping[str, Mask], out_path: str | Path, title: str):

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(6, 6))
    plt.imshow(
        raw,
        cmap="gray",
        vmin=np.percentile(raw, 0.1),
        vmax=np.percentile(raw, 99.7),
    )

    # Use fixed colors for readability.
    outline_colors = {
        "cell": "lime",
        "nucleus": "magenta",
        "cytoplasm": "cyan",
        "territory": "yellow",
        "pseudo_cytoplasm": "orange",
    }

    for name, mask in masks.items():
        outlines = label_boundaries(mask)
        color = outline_colors.get(name, "red")
        plt.contour(outlines.astype(float), levels=[0.5], colors=color, linewidths=0.5)

    plt.title(title, fontsize=14)
    plt.axis("off")
    plt.savefig(str(out_path), bbox_inches="tight", pad_inches=0, dpi=300)
    plt.close()

def label_boundaries(mask: Mask) -> np.ndarray:
    """Return a boolean image with object boundaries from a label mask."""
    mask = mask.astype(np.uint16, copy=False)
    boundaries = np.zeros(mask.shape, dtype=bool)

    boundaries[:-1, :] |= mask[:-1, :] != mask[1:, :]
    boundaries[1:, :] |= mask[1:, :] != mask[:-1, :]
    boundaries[:, :-1] |= mask[:, :-1] != mask[:, 1:]
    boundaries[:, 1:] |= mask[:, 1:] != mask[:, :-1]

    boundaries &= mask > 0
    return boundaries

def prop_dict(mask: Mask, intensity: Image) -> Dict[int, measure._regionprops.RegionProperties]:
    """Return regionprops indexed by label."""
    return {int(p.label): p for p in measure.regionprops(mask, intensity)}


def write_measurements_csv(
    csv_path: str | Path,
    nucleus_mask: Mask,
    nuclei_raw: Image,
    protein_raw: Optional[Image],
    pixel_size: float,
    population: str,
    source_file: str,
    cell_mask: Optional[Mask] = None,
    cytoplasm_mask: Optional[Mask] = None,
    territory_mask: Optional[Mask] = None,
    pseudo_cytoplasm_mask: Optional[Mask] = None,
):
    """Write per-object measurements to CSV.

    Nucleus measurements are always exported. Cell, cytoplasm, Voronoi territory
    and pseudo-cytoplasm measurements are exported only when the corresponding
    masks are available. Protein intensity columns are exported only when a
    protein channel is provided.
    """
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    intensity_for_shape = protein_raw if protein_raw is not None else nuclei_raw
    nuc_props = prop_dict(nucleus_mask, intensity_for_shape)
    cell_props = prop_dict(cell_mask, intensity_for_shape) if cell_mask is not None else {}
    cyto_props = prop_dict(cytoplasm_mask, intensity_for_shape) if cytoplasm_mask is not None else {}
    territory_props = prop_dict(territory_mask, intensity_for_shape) if territory_mask is not None else {}
    pseudo_props = prop_dict(pseudo_cytoplasm_mask, intensity_for_shape) if pseudo_cytoplasm_mask is not None else {}

    header = [
        "source_file",
        "image_id",
        "population",
        "object_id",
        "centroid_nucleus_y_px",
        "centroid_nucleus_x_px",
        "area_nucleus_um2",
        "nucleus_major_axis_px",
        "nucleus_minor_axis_px",
        "nucleus_aspect_ratio",
    ]
    if cell_mask is not None:
        header += ["area_cell_um2", "cell_major_axis_px", "cell_minor_axis_px", "cell_aspect_ratio"]
    if cytoplasm_mask is not None:
        header += ["area_cytoplasm_um2"]
    if territory_mask is not None:
        header += ["area_voronoi_territory_um2"]
    if pseudo_cytoplasm_mask is not None:
        header += ["area_pseudo_cytoplasm_um2"]
    if protein_raw is not None:
        header += ["mean_intensity_protein_nucleus"]
        if cell_mask is not None:
            header += ["mean_intensity_protein_cell"]
        if cytoplasm_mask is not None:
            header += ["mean_intensity_protein_cytoplasm"]
        if territory_mask is not None:
            header += ["mean_intensity_protein_voronoi_territory"]
        if pseudo_cytoplasm_mask is not None:
            header += ["mean_intensity_protein_pseudo_cytoplasm"]

    labels = sorted(nuc_props)
    with open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for label in labels:
            nuc = nuc_props[label]
            nuc_aspect = nuc.axis_major_length / nuc.axis_minor_length if nuc.axis_minor_length else MISSING
            row = [
                source_file,
                Path(source_file).stem,
                population,
                label,
                nuc.centroid[0],
                nuc.centroid[1],
                nuc.area * pixel_size * pixel_size,
                nuc.axis_major_length,
                nuc.axis_minor_length,
                nuc_aspect,
            ]
            if cell_mask is not None:
                cell = cell_props.get(label)
                if cell is not None:
                    cell_aspect = cell.axis_major_length / cell.axis_minor_length if cell.axis_minor_length else MISSING
                    row += [cell.area * pixel_size * pixel_size, cell.axis_major_length, cell.axis_minor_length, cell_aspect]
                else:
                    row += [MISSING, MISSING, MISSING, MISSING]
            if cytoplasm_mask is not None:
                cyto = cyto_props.get(label)
                row += [cyto.area * pixel_size * pixel_size if cyto is not None else MISSING]
            if territory_mask is not None:
                territory = territory_props.get(label)
                row += [territory.area * pixel_size * pixel_size if territory is not None else MISSING]
            if pseudo_cytoplasm_mask is not None:
                pseudo = pseudo_props.get(label)
                row += [pseudo.area * pixel_size * pixel_size if pseudo is not None else MISSING]
            if protein_raw is not None:
                row += [nuc.intensity_mean]
                if cell_mask is not None:
                    row += [cell_props[label].intensity_mean if label in cell_props else MISSING]
                if cytoplasm_mask is not None:
                    row += [cyto_props[label].intensity_mean if label in cyto_props else MISSING]
                if territory_mask is not None:
                    row += [territory_props[label].intensity_mean if label in territory_props else MISSING]
                if pseudo_cytoplasm_mask is not None:
                    row += [pseudo_props[label].intensity_mean if label in pseudo_props else MISSING]
            writer.writerow(row)


def process_image(filename: str | Path, params: Mapping, cell_model: Optional[object], nucleus_model: object) -> List[Path]:
    filename = Path(filename)
    print(f"\nProcessing {filename}", flush=True)

    raw = load_raw_image(filename)
    cell_raw, nuclei_channels, protein_raw = extract_channels(raw, params)

    segmentation_mode = get_param(params, "segmentation_mode", "cell_nucleus")
    use_voronoi = bool(get_param(params, "use_voronoi", False))
    pixel_size = float(get_param(params, "pixel_size", 1.0))
    output_root = Path(require_key(params, "output_directory"))
    out_dir = output_root / f"{filename.stem}_output"
    masks_dir = out_dir / "masks"
    overlays_dir = out_dir / "qc_overlays"
    out_dir.mkdir(parents=True, exist_ok=True)
    measurement_files: List[Path] = []

    filter_report: List[str] = [
        f"Filtering report for {filename.name}",
        f"Input file: {filename}",
        "",
    ]

    total_nuclei_for_cell = nuclei_channels.get("all_cells", next(iter(nuclei_channels.values())))
    cell_mask = None
    global_voronoi = None

    if segmentation_mode == "cell_nucleus":
        if cell_raw is None:
            raise ValueError("segmentation_mode is 'cell_nucleus', but no cell channel was provided.")
        if cell_model is None:
            raise ValueError("Cellpose cell model was not initialized, but cell segmentation was requested.")

        cell_mask = run_cellpose_cell_segmentation(cell_raw,total_nuclei_for_cell,params,cell_model)[0]    
        cell_mask = filter_by_area(cell_mask, cell_raw, float(get_param(params, "filter_area_cell", 0)),object_name="cells", report_lines=filter_report)
        save_mask(cell_mask, masks_dir / "cell_mask_initial.tif")

    for population, nuclei_raw in nuclei_channels.items():
        print(f"Segmenting {population}...", flush=True)
        nuclei_mask = run_cellpose_nuclear_segmentation(nuclei_raw,params,nucleus_model)[0]
        nuclei_mask = filter_nuclei(nuclei_mask, nuclei_raw, params, object_name=f"{population} nuclei", report_lines=filter_report)
        save_mask(nuclei_mask, masks_dir / f"{population}_nucleus_mask_initial.tif")

        final_cell = None
        final_nucleus = nuclei_mask
        cytoplasm = None
        territory = None
        pseudo_cytoplasm = None

        if segmentation_mode == "cell_nucleus":
            final_cell, final_nucleus, cytoplasm = associate_nuclei_to_cells(nuclei_mask, cell_mask, cell_raw, params,object_name=f"{population} nuclei",container_name="cells",report_lines=filter_report)
            save_mask(final_cell, masks_dir / f"{population}_cell_mask_final.tif")
            save_mask(final_nucleus, masks_dir / f"{population}_nucleus_mask_final.tif")
            save_mask(cytoplasm, masks_dir / f"{population}_cytoplasm_mask_final.tif")
            save_overlay(cell_raw + nuclei_raw, final_cell, overlays_dir / f"{population}_cell_overlay.png", f"{population}: cell masks")
            save_overlay(nuclei_raw, final_nucleus, overlays_dir / f"{population}_nucleus_overlay.png", f"{population}: nuclear masks")
            save_outline_overlay(nuclei_raw,{"nucleus": final_nucleus,},overlays_dir / f"{population}_nucleus_outlines.png",f"{population}: nuclear outlines",)
            save_outline_overlay(cell_raw,{"cell": final_cell},overlays_dir / f"{population}_cell_outlines.png",f"{population}: cell outlines",)
            save_outline_overlay(cell_raw + nuclei_raw,{"cell": final_cell,"nucleus": final_nucleus,},overlays_dir / f"{population}_cell_nucleus_outlines.png",f"{population}: cell and nuclear outlines",)

            if protein_raw is not None:
                save_overlay(protein_raw, cytoplasm, overlays_dir / f"{population}_cytoplasm_on_protein.png", f"{population}: cytoplasm on protein")
                save_outline_overlay(protein_raw,{"cell": final_cell,"nucleus": final_nucleus,},overlays_dir / f"{population}_cell_nucleus_outlines_on_protein.png",f"{population}: cell and nuclear outlines on protein",)
        else:
            save_mask(final_nucleus, masks_dir / f"{population}_nucleus_mask_final.tif")
            save_overlay(nuclei_raw, final_nucleus, overlays_dir / f"{population}_nucleus_overlay.png", f"{population}: nuclear masks")
            save_outline_overlay(nuclei_raw,{"nucleus": final_nucleus},overlays_dir / f"{population}_nucleus_outlines.png",f"{population}: nuclear outlines")
            if protein_raw is not None:
                save_overlay(protein_raw, final_nucleus, overlays_dir / f"{population}_nucleus_on_protein.png", f"{population}: nucleus on protein")
                save_outline_overlay(protein_raw,{"nucleus": final_nucleus},overlays_dir / f"{population}_nucleus_outlines_on_protein.png",f"{population}: nuclear outlines on protein")

        if use_voronoi:
            fg_source = protein_raw if protein_raw is not None else total_nuclei_for_cell

            # In dual-population mode, use one global Voronoi generated from all nuclei.
            if get_nuclear_mode(params) == "dual_population":
                if global_voronoi is None:
                    if population != "all_cells":
                        raise RuntimeError(
                            "Expected 'all_cells' to be processed before population-specific masks."
                        )

                    if bool(get_param(params, "use_foreground_mask", False)):
                        fg_mask = make_foreground_mask(fg_source, params)
                        save_mask(fg_mask.astype(np.uint16), masks_dir / "all_cells_foreground_mask.tif")
                    else:
                        fg_mask = None

                    global_voronoi = make_voronoi_territories(final_nucleus, fg_mask, params)

                if population == "all_cells":
                    territory = global_voronoi
                    pseudo_cytoplasm = make_pseudo_cytoplasm(territory, final_nucleus, params)
                else:
                    # Use global Voronoi as pseudo-cell territories and assign this
                    # population's nuclei to those territories.
                    territory, final_nucleus, pseudo_cytoplasm = associate_nuclei_to_cells(
                        final_nucleus,
                        global_voronoi,
                        nuclei_raw,
                        params,
                        object_name=f"{population} nuclei",
                        container_name="Voronoi territories",
                        report_lines=filter_report
                    )

                save_mask(territory, masks_dir / f"{population}_voronoi_territory_mask.tif")
                save_mask(pseudo_cytoplasm, masks_dir / f"{population}_pseudo_cytoplasm_mask.tif")
                save_overlay(
                    fg_source,
                    territory,
                    overlays_dir / f"{population}_voronoi_territory_overlay.png",
                    f"{population}: Voronoi territories",
                )
                save_outline_overlay(
                    fg_source,
                    {
                        "territory": territory,
                        "nucleus": final_nucleus,
                    },
                    overlays_dir / f"{population}_voronoi_nucleus_outlines.png",
                    f"{population}: Voronoi and nuclear outlines",
                )
                if protein_raw is not None:
                    save_overlay(
                        protein_raw,
                        pseudo_cytoplasm,
                        overlays_dir / f"{population}_pseudo_cytoplasm_on_protein.png",
                        f"{population}: pseudo-cytoplasm on protein",
                    )

            else:
                # Single-population case: Voronoi from the current nuclear mask.
                if bool(get_param(params, "use_foreground_mask", False)):
                    fg_mask = make_foreground_mask(fg_source, params)
                    save_mask(fg_mask.astype(np.uint16), masks_dir / f"{population}_foreground_mask.tif")
                else:
                    fg_mask = None

                territory = make_voronoi_territories(final_nucleus, fg_mask, params)
                pseudo_cytoplasm = make_pseudo_cytoplasm(territory, final_nucleus, params)

                save_mask(territory, masks_dir / f"{population}_voronoi_territory_mask.tif")
                save_mask(pseudo_cytoplasm, masks_dir / f"{population}_pseudo_cytoplasm_mask.tif")
                save_overlay(
                    fg_source,
                    territory,
                    overlays_dir / f"{population}_voronoi_territory_overlay.png",
                    f"{population}: Voronoi territories",
                )
                save_outline_overlay(
                    fg_source,
                    {
                        "territory": territory,
                        "nucleus": final_nucleus,
                    },
                    overlays_dir / f"{population}_voronoi_nucleus_outlines.png",
                    f"{population}: Voronoi and nuclear outlines",
                )
                if protein_raw is not None:
                    save_overlay(
                        protein_raw,
                        pseudo_cytoplasm,
                        overlays_dir / f"{population}_pseudo_cytoplasm_on_protein.png",
                        f"{population}: pseudo-cytoplasm on protein",
                    )

        measurement_csv = out_dir / f"{population}_measurements.csv"
        write_measurements_csv(
            measurement_csv,
            final_nucleus,
            nuclei_raw,
            protein_raw,
            pixel_size,
            population,
            source_file=filename.name,
            cell_mask=final_cell,
            cytoplasm_mask=cytoplasm,
            territory_mask=territory,
            pseudo_cytoplasm_mask=pseudo_cytoplasm,
        )
        measurement_files.append(measurement_csv)
        print(f"Saved outputs for {population} in {out_dir}", flush=True)

    write_filter_report(filter_report, out_dir / "filtering_report.txt")
    return measurement_files


def combine_measurement_tables(csv_paths: Iterable[str | Path], out_path: str | Path) -> Optional[Path]:
    """Concatenate per-image/per-population measurement CSVs into one table.

    The combined table is written to the main output directory. The function uses
    the union of all columns, so it can safely combine tables generated with or
    without protein quantification or Voronoi measurements, although a single run
    normally produces a consistent set of columns.
    """
    rows: List[Dict[str, str]] = []
    fieldnames: List[str] = []
    for csv_path in csv_paths:
        csv_path = Path(csv_path)
        if not csv_path.exists():
            continue
        with open(csv_path, "r", newline="") as handle:
            reader = csv.DictReader(handle)
            for name in reader.fieldnames or []:
                if name not in fieldnames:
                    fieldnames.append(name)
            rows.extend(dict(row) for row in reader)

    if not rows:
        return None

    preferred = ["source_file", "image_id", "population", "object_id"]
    ordered_fields = [name for name in preferred if name in fieldnames]
    ordered_fields += [name for name in fieldnames if name not in ordered_fields]

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered_fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, MISSING) for name in ordered_fields})
    return out_path


def iter_input_images(folder: str | Path) -> Iterable[Path]:
    folder = Path(folder)
    paths = [
        path
        for path in folder.iterdir()
        if path.is_file() and path.suffix.lower() in {".tif", ".tiff"}
    ]
    return sorted(paths)

def save_run_parameters(params: Mapping, out_path: str | Path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as handle:
        yaml.safe_dump(dict(params), handle, sort_keys=False)

def write_filter_report(report_lines: List[str], path: str | Path):
    """Write filtering/report messages to a text file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as handle:
        for line in report_lines:
            handle.write(str(line) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Run unified segmentation pipeline.")
    parser.add_argument("config", help="Path to YAML configuration file.")
    parser.add_argument(
        "--test",
        type=int,
        default=None,
        help="Run the pipeline only on the first N images for parameter testing.",
    )
    parser.add_argument(
        "--test-output-suffix",
        default="_test",
        help="Suffix added to the output directory when running in test mode.",
    )
    args = parser.parse_args()

    with open(args.config, "r") as handle:
        params = yaml.safe_load(handle)

    folder = Path(require_key(params, "folder_path"))
    images = list(iter_input_images(folder))
    if not images:
        raise FileNotFoundError(f"No TIFF images found in {folder}")
    
    if args.test is not None:
        if args.test <= 0:
            raise ValueError("--test must be a positive integer.")

        images = images[: args.test]

        params = dict(params)
        original_output = Path(require_key(params, "output_directory"))
        params["output_directory"] = str(original_output.with_name(original_output.name + args.test_output_suffix))

        print(
            f"Test mode enabled: processing only {len(images)} image(s). "
            f"Outputs will be saved to {params['output_directory']}",
            flush=True,
        )
    output_root = Path(require_key(params, "output_directory"))
    save_run_parameters(params, output_root / "run_parameters.yml")

    cell_model, nucleus_model = initialize_cellpose_models(params)
    all_measurement_files: List[Path] = []
    for image_path in tqdm(images, desc="Images"):
        all_measurement_files.extend(process_image(image_path, params, cell_model, nucleus_model))

    if bool(get_param(params, "write_combined_table", True)):
        output_root = Path(require_key(params, "output_directory"))
        combined = combine_measurement_tables(all_measurement_files, output_root / "combined_measurements.csv")
        if combined is not None:
            print(f"Combined measurement table saved to {combined}", flush=True)
        else:
            print("No measurement rows found; combined table was not written.", flush=True)


if __name__ == "__main__":
    main()
