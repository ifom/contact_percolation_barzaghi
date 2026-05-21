# Cell Motility Analysis Pipeline

## Overview
This folder contains analysis pipelines for studying cell motility in mixed cancer cell populations, including trajectory tracking, velocity analysis, and statistical characterization of cell movement.

## Folder Contents

- **MixedPop_Code4Sharing/**: Main analysis pipeline with MATLAB and Python scripts
  - `DYNAMIX_*.m`: Trajectory building, velocity calculations, MSD analysis
  - `STATIX_*.m`: Statistical analysis and clustering
  - `PIV_*.m`: Particle Image Velocimetry analysis
  - `Functions/`: Helper functions for data processing
  - `EdU_quantification_*.py`: Automated EdU staining quantification
  - `Macro_Fiji_*.ijm`: ImageJ/FIJI macros for segmentation

- **Image_analisys/**: Image segmentation pipeline (Python)
  - Cell segmentation using neural networks
  
- **RNA_seq/**: RNA sequencing data and analysis

## Quick Start

1. **Image Segmentation**: Run `Image_analisys/segmentation.py` to segment cell regions
2. **Trajectory Analysis**: Execute MATLAB scripts in numerical order (01-07) in `MixedPop_Code4Sharing/`
3. **Statistical Analysis**: Run `STATIX_*.m` scripts for clustering and statistics

## Notes
- Code is provided as-is for the original study; parameters may need adjustment for new datasets
- Requires MATLAB with Image Processing Toolbox, Python dependencies in `RNA_seq_env.yml`
- See individual script headers for specific requirements
