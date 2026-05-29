# Unified Cellpose-SAM segmentation pipeline

This repository contains custom analysis scripts used to segment fluorescence microscopy images and quantify per-object protein intensity. The code is provided to document the computational workflow used for the manuscript and to support reproducibility. It is not intended as a general-purpose image-analysis software package.

This version is adapted to segmenting cells using **Cellpose-SAM**. 

Pachitariu, M., Rariden, M., & Stringer, C. (2025). Cellpose-SAM: superhuman generalization for cellular segmentation. <em>bioRxiv</em>. [Link to the original paper](https://www.biorxiv.org/content/10.1101/2025.04.28.651001v1)

## Overview

The main script is:

```bash
scripts/unified_segmentation.py
```

The pipeline supports:

| Mode | Description |
|---|---|
| `cell_nucleus` + `single` | Segment cells and one nuclear channel |
| `cell_nucleus` + `dual_population` | Segment cells and two population-specific nuclear channels |
| `nucleus_only` + `single` | Segment one nuclear channel without a cell marker |
| `nucleus_only` + `dual_population` | Segment two population-specific nuclear channels without a cell marker |

If `quantify_protein: true`, the protein/intensity channel is quantified in the available compartments: nucleus, whole cell, cytoplasm, nuclear territory and/or pseudo-cytoplasm, depending on the selected mode.

When no cell/cytoplasmic marker is available, the pipeline can generate watershed-derived nuclear territories as pseudo-cell compartments.

## Installation

Create and activate the base conda environment:

```bash
conda env create -f environment.yml
conda activate segmentation-pipeline-cellposeSAM
```

The `environment.yml` file creates a minimal Python environment. After activating it, install the PyTorch backend appropriate for your hardware **before installing Cellpose**, especially if GPU acceleration is required.

The pipeline dependencies, including Cellpose, are installed with:

```bash
python -m pip install -r requirements.txt
```

### CPU installation

For CPU-only use, install the requirements directly:

```bash
python -m pip install -r requirements.txt
```

Use this in the YAML:

```yaml
use_gpu: false
gpu_device: "auto"
```

This is the safest option on laptops or machines without a supported GPU.

### NVIDIA CUDA GPU

For NVIDIA GPU acceleration, install a PyTorch build compatible with your GPU and driver **before** installing Cellpose.

For recent NVIDIA Blackwell GPUs with compute capability `sm_120`, a CUDA 12.8+ PyTorch build may be required. One possible starting point is:

```bash
python -m pip install --pre torch torchvision --index-url https://download.pytorch.org/whl/nightly/cu128
python -m pip install -r requirements.txt
```

Then use:

```yaml
use_gpu: true
gpu_device: "cuda"
```

For other NVIDIA GPUs, choose the appropriate PyTorch command from the official PyTorch installation selector.

### Apple Silicon Mac: MPS backend

On Apple Silicon Macs, Cellpose-SAM can use the PyTorch Metal/MPS backend instead of CUDA. Install the requirements normally:

```bash
python -m pip install -r requirements.txt
```

Then verify MPS availability:

```bash
python - <<'PY'
import torch
print("MPS available:", torch.backends.mps.is_available())
print("MPS built:", torch.backends.mps.is_built())
PY
```

To test Apple Silicon acceleration, use:

```yaml
use_gpu: true
gpu_device: "mps"
```

If MPS execution fails, is unstable, or gives unexpected results, use CPU mode instead:

```yaml
use_gpu: false
gpu_device: "auto"
```

### AMD ROCm GPU

For AMD GPUs on Linux, install the ROCm drivers and a ROCm-compatible PyTorch build first, then install Cellpose and the pipeline requirements:

```bash
python -m pip install -r requirements.txt
```

ROCm support is less mature than CUDA and may require additional troubleshooting.

## Verify the installation

After installation, check the environment:

```bash
python - <<'PY'
import torch, cellpose

print("torch:", torch.__version__)
print("torch CUDA:", torch.version.cuda)
print("CUDA available:", torch.cuda.is_available())

if torch.cuda.is_available():
    print("GPU:", torch.cuda.get_device_name(0))
    print("capability:", torch.cuda.get_device_capability(0))

if hasattr(torch.backends, "mps"):
    print("MPS available:", torch.backends.mps.is_available())
    print("MPS built:", torch.backends.mps.is_built())
else:
    print("MPS available:", False)
    print("MPS built:", False)

print("cellpose:", getattr(cellpose, "__version__", getattr(cellpose, "version_str", "unknown")))
PY
```

If GPU execution is not available or not reliable, set:

```yaml
use_gpu: false
```

## Input image layout

The pipeline assumes channel-first TIFF files:

```text
single-plane multichannel image: (C, Y, X)
z-stack multichannel image:      (Z, C, Y, X)
```

For z-stacks, each requested channel is extracted along the C axis and projected along Z.

Channel indices in the YAML always refer to the C axis.

```yaml
z_projection: "max"           # used for cell and nuclear channels
protein_z_projection: "mean"  # used for the protein/intensity channel
```

## Running the pipeline

Before processing a full batch, test the configuration on one or a few images:

```bash
python scripts/unified_segmentation.py config/params.yml --test 1
```

Run the full batch:

```bash
python scripts/unified_segmentation.py config/params.yml
```

In test mode, outputs are written to a separate folder with a `_test` suffix, so the final output directory is not overwritten.

## Key YAML options

```yaml
segmentation_mode: "cell_nucleus"    # "cell_nucleus" or "nucleus_only"
nuclear_mode: "dual_population"      # "single" or "dual_population"
quantify_protein: true               # true or false
use_voronoi: false                   # true or false

# Hardware acceleration
use_gpu: false                       # false forces CPU execution
gpu_device: "auto"                  # "auto", "cuda", "mps", or "cpu"
```

Example: cell + two nuclear populations + protein quantification

```yaml
segmentation_mode: "cell_nucleus"
nuclear_mode: "dual_population"
quantify_protein: true

channels:
  protein: 0
  population_1_nucleus: 1
  population_2_nucleus: 2
  cell: 3
```

Example: nucleus-only image with one nuclear channel and protein quantification

```yaml
segmentation_mode: "nucleus_only"
nuclear_mode: "single"
quantify_protein: true
use_voronoi: false

channels:
  nucleus: 0
  protein: 1
```

Example: nucleus-only dual-population image with nuclear territories

```yaml
segmentation_mode: "nucleus_only"
nuclear_mode: "dual_population"
quantify_protein: true
use_voronoi: true

channels:
  total_nuclei: 0
  population_1_nucleus: 1
  population_2_nucleus: 2
  protein: 3
```

If `total_nuclei` is not provided in dual-population mode, the total nuclear signal is generated by summing the two population-specific nuclear channels.

## Nuclear territories and pseudo-cytoplasm

When `use_voronoi: true`, the pipeline generates nuclear territories to approximate pseudo-cell regions in images lacking a cell/cytoplasmic marker.

Two seed modes are supported:

```yaml
voronoi_seed_mode: "nuclear_mask"  # "nuclear_mask" or "centroid"
```

- `nuclear_mask`: territories are generated by watershed expansion from labelled nuclear masks.
- `centroid`: territories are generated from nuclear centroids.

The pseudo-cytoplasm is defined as:

```text
pseudo-cytoplasm = nuclear territory - nuclear mask
```

If `voronoi_ring_width_px` is set to a positive value, the pseudo-cytoplasm is restricted to a perinuclear ring. If it is `0`, the full territory outside the nucleus is used.

For dual-population analyses, a global nuclear territory mask is generated from all nuclei and then assigned to each population by spatial association with population-specific nuclei.

## Filtering and object association

After segmentation, nuclei are associated with cell masks or nuclear territories. The association step removes likely segmentation errors, including nuclei that do not overlap any cell/territory, nuclei shared ambiguously between masks, and cells/territories without a valid nucleus.

Small marginal overlaps between nuclei and adjacent cells can be ignored before assignment:

```yaml
min_container_overlap_fraction: 0.02
```

Cells or territories containing more than one valid nucleus are handled according to:

```yaml
multi_nucleus_cell_policy: "remove_cell"  # "remove_cell", "keep_largest", or "keep_all"
```

For conservative single-cell measurements, `remove_cell` is recommended.

A filtering summary is saved for each image as:

```text
filtering_report.txt
```

## Outputs

For each input TIFF, outputs are written to:

```text
<output_directory>/<image_name>_output/
```

Main outputs include:

```text
masks/*_nucleus_mask_initial.tif
masks/*_nucleus_mask_final.tif
masks/*_cell_mask_final.tif                 # cell_nucleus mode
masks/*_cytoplasm_mask_final.tif            # cell_nucleus mode
masks/*_voronoi_territory_mask.tif          # if use_voronoi: true
masks/*_pseudo_cytoplasm_mask.tif           # if use_voronoi: true
qc_overlays/*.png
*_measurements.csv
filtering_report.txt
```

For each batch, the script also writes:

```text
<output_directory>/combined_measurements.csv
<output_directory>/run_parameters.yml
<output_directory>/environment_report.txt
```

`combined_measurements.csv` concatenates all per-image and per-population measurement tables and includes `source_file` and `image_id`.

To disable the combined table:

```yaml
write_combined_table: false
```
Credits to Leonardo Barzaghi
