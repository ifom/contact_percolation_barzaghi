# Cell Motility Analysis Pipeline

## Overview
This folder contains analysis pipelines for studying cell motility in mixed cancer cell populations, including trajectory tracking, velocity analysis, and statistical characterization of cell movement.

## Folder Contents

- **MixedPop_Code4Sharing/**: Main analysis pipeline with MATLAB and Python scripts
  - `DYNAMIX_*.m`: Trajectory building, velocity calculations, MSD analysis
  - `STATIX_*.m`: Statistical analysis and clustering
  - `PIV_*.m`: Particle Image Velocimetry analysis
  - [`Functions/`](./MixedPop_Code4Sharing/Functions/): Helper functions for data processing
  - `Macro_Fiji_*.ijm`: ImageJ/FIJI macros for segmentation

## Quick Start

1. **Trajectory Analysis**: Execute MATLAB scripts in numerical order (01-07) in `MixedPop_Code4Sharing/`
2. **Statistical Analysis**: Run `STATIX_*.m` scripts for clustering and statistics