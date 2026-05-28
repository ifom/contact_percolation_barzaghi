# Data dictionary

Each image produces one or more per-image measurement tables named `<population>_measurements.csv` inside:

```text
<output_directory>/<image_name>_output/
```

At the end of the run, these tables are concatenated into:

```text
<output_directory>/combined_measurements.csv
```

The combined table contains all rows from all segmented images and includes the original file name for traceability.

Missing values are written as empty fields.

## Populations

| Value | Meaning |
|---|---|
| `all_cells` | All segmented nuclei/cells. In single-nuclear-channel mode, this corresponds to the single nuclear channel. In dual-population mode, this corresponds to the total nuclear channel if provided, or to the sum of the two population-specific nuclear channels. |
| `population_1` | Objects associated with `population_1_nucleus`. |
| `population_2` | Objects associated with `population_2_nucleus`. |

## Always present

| Column | Meaning |
|---|---|
| `source_file` | Name of the input TIFF file that generated the measurement row. |
| `image_id` | Input file name without extension. Useful for grouping fields of view. |
| `population` | Analysed population: `all_cells`, `population_1`, or `population_2`. |
| `object_id` | Label ID of the measured object. In `cell_nucleus` mode, final nuclear labels are matched to their associated cell labels. |
| `centroid_nucleus_y_px` | Nuclear centroid y-coordinate in pixels. |
| `centroid_nucleus_x_px` | Nuclear centroid x-coordinate in pixels. |
| `area_nucleus_um2` | Nuclear area in square micrometers. |
| `nucleus_major_axis_px` | Major axis length of the nuclear mask in pixels. |
| `nucleus_minor_axis_px` | Minor axis length of the nuclear mask in pixels. |
| `nucleus_aspect_ratio` | Nuclear major/minor axis ratio. |

## Present in `cell_nucleus` mode

| Column | Meaning |
|---|---|
| `area_cell_um2` | Area of the associated cell mask in square micrometers. |
| `cell_major_axis_px` | Major axis length of the cell mask in pixels. |
| `cell_minor_axis_px` | Minor axis length of the cell mask in pixels. |
| `cell_aspect_ratio` | Cell major/minor axis ratio. |
| `area_cytoplasm_um2` | Area of the cell mask after excluding the associated nuclear mask. |

## Present when `use_voronoi: true`

These columns describe watershed-derived nuclear territories used as pseudo-cell compartments when no cell/cytoplasmic marker is available.

| Column | Meaning |
|---|---|
| `area_voronoi_territory_um2` | Area of the nuclear territory assigned to the object. Depending on `voronoi_seed_mode`, territories are generated from nuclear centroids or by watershed expansion from labelled nuclear masks. |
| `area_pseudo_cytoplasm_um2` | Area of the nuclear territory after excluding the nuclear mask, or the perinuclear ring if `voronoi_ring_width_px > 0`. |

## Present when `quantify_protein: true`

Protein intensities are measured on the extracted protein/intensity channel after the selected z-projection (`protein_z_projection`, for z-stacks).

| Column | Meaning |
|---|---|
| `mean_intensity_protein_nucleus` | Mean protein-channel intensity inside the nuclear mask. |
| `mean_intensity_protein_cell` | Mean protein-channel intensity in the associated cell mask. Present only in `cell_nucleus` mode. |
| `mean_intensity_protein_cytoplasm` | Mean protein-channel intensity in the cell-derived cytoplasmic mask. Present only in `cell_nucleus` mode. |
| `mean_intensity_protein_voronoi_territory` | Mean protein-channel intensity in the nuclear territory. Present only if `use_voronoi: true`. |
| `mean_intensity_protein_pseudo_cytoplasm` | Mean protein-channel intensity in the pseudo-cytoplasmic region. Present only if `use_voronoi: true`. |

## Related QC files

Each image output folder also contains segmentation masks, QC overlays and a filtering report.

```text
masks/
qc_overlays/
filtering_report.txt
```

`filtering_report.txt` reports how many objects were removed during area/eccentricity filtering and during nucleus-to-cell or nucleus-to-territory association.
