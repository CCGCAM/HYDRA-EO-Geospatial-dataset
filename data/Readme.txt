HYDRA-EO – Data Folder

This folder contains all datasets used in the HYDRA-EO project, organized by processing level, sensor type, and data origin.

Structure:

- raw/
  Original, unprocessed datasets as acquired from different sources:
  - Satellite imagery (e.g., Sentinel-2, PRISMA, EnMAP)
  - Airborne and UAV data
  - Field measurements

- interim/
  Intermediate datasets generated during preprocessing:
  - Temporally aggregated data
  - Pre-processed hyperspectral imagery
  - Harmonized multi-sensor datasets

- processed/
  Final, analysis-ready datasets:
  - Atmospherically corrected imagery
  - Derived products (e.g., vegetation indices, plant traits)
  - Model outputs

- Sentinel-2/
  Sentinel-2 specific datasets:
  - Time series imagery
  - Regions of Interest (ROIs)
  - Aggregated and processed products for temporal analysis


- rasters/
  Raster datasets used across workflows:
  - Multispectral and hyperspectral imagery
  - Derived raster layers (e.g., NDVI, SIF, temperature)

- grids/
  Spatial grid definitions used for analysis:
  - Tiling systems
  - Aggregation units (e.g., pixels, ESUs)

- Others/
  Additional datasets from collaborators or external sources:
  - Ancillary data (e.g., weather, soil)
  - External or experimental datasets not yet integrated

Notes:
- All raw data must remain unchanged.
- Processing workflows should clearly separate interim and final products.
- Ensure consistent coordinate reference systems (CRS) and metadata across datasets.
- Data naming conventions and versioning should be applied when updating datasets.