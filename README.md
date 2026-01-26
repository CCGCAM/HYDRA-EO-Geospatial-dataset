# HYDRA-EO Geo-spatial  Dataset

<p align="center">

![HYDRA-EO logo](assets/Hydra-eo-logo.png)

</p>

**Geospatial analysis and methodological exploration for multi-stressor crop monitoring using hyperspectral, thermal, and multispectral Earth Observation data across UAV, airborne, and satellite platforms.**

This repository serves as an **exploratory geospatial sandbox** within the ESA **HYDRA-EO** initiative. It hosts open resources for testing, intercomparison, and methodological development, including geospatial data schemas, STAC metadata, analysis scripts, notebooks, and scientific documentation.

⚠️ **Important**\

This repository is **not** intended for operational processing chains. It is designed for **testing concepts, exploring synergies, and developing and validating methods** that may later feed into operational or project-level repositories.

------------------------------------------------------------------------

## Scope and Focus

This repository focuses on **geospatial data analysis and integration** for agricultural monitoring, with emphasis on:

-   Multi-sensor Earth Observation datasets (hyperspectral, thermal, multispectral)
-   Spatial–temporal harmonization across UAV, airborne, and satellite scales
-   STAC-based data cataloguing and discovery
-   Reproducible geospatial workflows in R and Python

The repository is designed to support **exploratory analysis, validation, and algorithm development** rather than operational processing pipelines.

------------------------------------------------------------------------

## Key Objectives

The geospatial analysis objectives of HYDRA-EO include:

1.  **Spatial characterization of crop stress patterns** using EO-derived variables (reflectance, vegetation indices, canopy temperature, SIF).
2.  **Multi-scale data integration**, enabling comparison and upscaling from UAV and airborne imagery to satellite observations (Sentinel-2, PRISMA, EnMAP, FLEX, CHIME).
3.  **Sensor intercomparison and harmonization**, including spectral and spatial resampling across platforms.
4.  **Linkage of EO products with field observations**, supporting validation and spatial attribution of stress signals.
5.  **Provision of FAIR geospatial datasets and metadata**, aligned with STAC, ESA EarthCODE, and open-science practices.

------------------------------------------------------------------------

## Repository Contents

This repository provides reusable geospatial resources for the EO and crop science community:

-   **STAC Catalogs**\
    Static STAC catalogs describing EO datasets, collections, and acquisitions across sensors and sites (`stac/`).

-   **Geospatial Datasets (references)**\
    Metadata and pointers to GeoTIFF/COG raster products (hyperspectral, thermal, multispectral) and associated quicklooks. Large datasets are hosted externally and referenced via STAC.

-   **Analysis Scripts**\
    R and Python scripts for raster processing, spatial extraction, sensor harmonization, and exploratory geospatial analysis (`scripts/`).

-   **Notebooks & Tutorials**\
    R Markdown and Jupyter notebooks demonstrating geospatial workflows, including pixel extraction, time-series analysis, and multi-sensor comparisons (`notebooks/`).

-   **Documentation & Roadmap**\
    Methodological notes, analysis specifications, and the scientific roadmap for geospatial exploitation of HYDRA-EO datasets (`docs/`, `routemap/`).

------------------------------------------------------------------------

## Repository Structure

``` text
HYDRA-EO/
├─ assets/                # figures, logos, visual material
├─ stac/                  # STAC catalog, collections, and items
├─ data/                  # metadata, inventories, and data pointers
│  ├─ raw/                # raw acquisitions (not tracked)
│  ├─ interim/            # intermediate geospatial products
│  └─ processed/          # example outputs / derived products
├─ scripts/               # geospatial analysis code (R / Python)
│  ├─ R/
│  └─ python/
├─ notebooks/             # geospatial analysis notebooks & tutorials
```

## Radiative Transfer Modeling in HYDRA-EO

The HYDRA-EO project builds upon two in-house R packages developed and maintained by the team, which form the backbone of the synthetic simulations used in this repository:

### 🔹 [ToolsRTM](https://gitlab.com/caminoccg/toolsrtm)

ToolsRTM is a comprehensive R package for structural radiative transfer modeling at canopy level. It supports PROSAIL and other 1D RTM simulations, enables the generation of look-up tables (LUTs), and allows band resampling to different sensors such as Sentinel-2, PRISMA, EnMAP, and CHIME.

The package also includes utilities for computing vegetation indices, BRF, and sensitivity analyses, along with Shiny applications for interactive trait–reflectance exploration.

### 🔹 [SCOPEinR](https://gitlab.com/caminoccg/scopeinr)

SCOPEinR is an R implementation of the SCOPE model for functional radiative transfer modeling. It simulates photosynthesis, sun-induced fluorescence (SIF), gross primary production (GPP), canopy temperature, and transpiration.

SCOPEinR links physiological traits with reflectance and flux measurements and integrates meteorological drivers from flux towers. It also provides Shiny applications for interactive exploration of functional traits and outputs.

------------------------------------------------------------------------

Together, **ToolsRTM + SCOPEinR** allow HYDRA-EO to generate **synthetic datasets** that couple **structural (reflectance)** and **functional (photosynthesis, SIF)** signals, providing a robust foundation for algorithm validation, stress detection, and multi-sensor data integration in the ESA monitoring framework.

### Manuals

The manuals are accessible through the [Shiny app](https://carlos-camino.shinyapps.io/0-toolsrtm-simulator/) or directly within the [ToolsRTM](https://carlos-camino.shinyapps.io/0-toolsrtm-simulator/_w_ef4421a7/Notebooks/R/ToolsRTM/ToolsRTM.html) and [SCOPEinR](https://carlos-camino.shinyapps.io/0-toolsrtm-simulator/_w_ef4421a7/Notebooks/R/SCOPEinR/SCOPEinR.html) packages. Vignettes are currently under development.

### Citation

If you use **ToolsRTM** packages, please cite the following references:

Camino et al., (2024). RT-Simulator: An Online Platform to Simulate Canopy Reflectance from Biochemical and Structural Plant Properties Using Radiative Transfer Models, *IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium*, Athens, Greece, 2024, pp. 2811-2814, [doi: 10.1109/IGARSS53475.2024.10642442](https://ieeexplore.ieee.org/document/10642442).

Arano et al., (2024). Enhancing Chlorophyll Content Estimation With Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical Models, *IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium*, Athens, Greece, 2024, pp. 4486-4489, [doi: 10.1109/IGARSS53475.2024.10641613](https://ieeexplore.ieee.org/document/10641613).

Camino et al., (in prep). Integrating physiological plant traits with Sentinel-2 imagery for monitoring gross primary production and detecting forest disturbances.

### License

[![](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

The **ToolsRTM** and **SCOPEinR** package is licensed under the MIT License, allowing for free use, modification, and distribution. This package is available on GitLab, and we encourage contributions and collaborations from the community. For more details, please refer to the LICENSE file in the repository.
