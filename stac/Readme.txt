HYDRA-EO – STAC Metadata

This folder contains the STAC (SpatioTemporal Asset Catalog) metadata for the HYDRA-EO geospatial datasets.

Purpose:
The STAC structure enables standardized data discovery, access, and interoperability across satellite, airborne, and field datasets. It is designed to support scalable analysis workflows and integration with cloud-based platforms.

Structure:

- catalog.json
  Root catalog describing the overall HYDRA-EO dataset.

- collections/
  Collection-level metadata grouping datasets by sensor or data type (e.g., Sentinel-2, PRISMA, airborne, field).

- items/
  Item-level metadata describing individual scenes, acquisitions, or observations.
  Each item includes:
  - Spatial extent (geometry, bounding box)
  - Temporal information (acquisition date/time)
  - Links to assets (data files)
  - Sensor/platform metadata

- assets/
  References to actual data files stored in the repository or external locations (e.g., GeoTIFF, NetCDF).

Standards:
- STAC version: follow the official STAC specification (https://stacspec.org)
- Coordinate Reference System (CRS): use standard EPSG codes
- Time format: ISO 8601 (UTC)
- Metadata fields should follow STAC core and relevant extensions (e.g., EO, Raster, Projection)

Guidelines:
- Each dataset (satellite, airborne, field) should be represented as a STAC Collection.
- Each acquisition or observation should be defined as a STAC Item.
- Assets must include clear roles (e.g., data, thumbnail, metadata).
- Ensure consistency between STAC metadata and actual data files.

Notes:
- The STAC catalog is intended to be machine-readable and compatible with tools such as pystac, stac-browser, and cloud-native geospatial platforms.
- Keep metadata updated when datasets are modified or new data are added.