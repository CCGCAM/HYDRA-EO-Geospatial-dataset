
#
rm(list= ls())
raster_temp_dir<-paste(tempdir(),'raster/',sep='')
files_to_be_removed<-list.files (raster_temp_dir,full.names = T)
file.remove(files_to_be_removed)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#	0. Libraries   -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ---------------------------
# Step 1: Define Required Packages
# ---------------------------

required_packages <- c(
  "raster",
  "RColorBrewer",
  "pls",
  "signal",
  "MASS",
  "caret",
  "ggplot2",
  "sf",
  "svMisc",
  "terra",
  "dplyr"
)
# ---------------------------
# Step 1.1: Load All Required Packages
# ---------------------------

invisible(lapply(c(required_packages,'tinytex'), library, character.only = TRUE))

# ------------------------------------------------------------------------------
# Step 2: Install Missing Packages (for Latex and pdf generation)
# ------------------------------------------------------------------------------

installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
} else {
  cat(" All required packages are already installed.\n")
}

# ------------------------------------------------------------------------------
# Step 3: Install & Load ToolsRTM from GitLab
# ------------------------------------------------------------------------------

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  # Public GitLab repo (no token needed). Set upgrade="never" for reproducibility.
  remotes::install_gitlab("caminoccg/toolsrtm", upgrade = "never")
}
library(ToolsRTM)

cat("\n ToolsRTM is ready: ", as.character(packageVersion("ToolsRTM")), "\n", sep = "")

# ------------------------------------------------------------------------------
# Step 4: Install & Load ToolsRTM from GitLab
# ------------------------------------------------------------------------------

if (!requireNamespace("SCOPEinR", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  # Public GitLab repo (no token needed). Set upgrade="never" for reproducibility.
  remotes::install_gitlab("caminoccg/scopeinr", upgrade = "never")
}
library(SCOPEinR)
cat("\n ToolsRTM is ready: ", as.character(packageVersion("ToolsRTM")), "\n", sep = "")


source('scripts/functions/getSeries_v2.R')
source('scripts/functions/getIndices_SE2.R')

# ------------------------------------------------------------------------------
# 1. Sentinel-2 band configuration
# ------------------------------------------------------------------------------

SE_bands <- c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12","SCL")
rfl.bands <- c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12")
factorR <- 1 / 10000


# ------------------------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------------------------

paths.se2a='data/Sentinel-2/Series/Chaparrillo/' # Scenario for fungus outbreaks
paths.shape  <- "data/sites/Olive-Pistachio-CIAG_trees/"
paths.outs  <- "Tables/results/"

if (!dir.exists(paths.outs)) dir.create(paths.outs, recursive = TRUE)

# ------------------------------------------------------------------------------
# 3. Read Sentinel-2 files
# ------------------------------------------------------------------------------

#paths.se2a='Sentinel-2/BarkBeetle/2016-2021/' # Scenario for Bark Beetle outbreaks

files.SE2 <- list.files(paths.se2a, full.names = T) ###Stack SE2

SE2.data <-raster::stack(files.SE2[2]) ###Stack SE2
names(SE2.data) <-SE_bands
SE2.data <-SE2.data * factorR
SE2.data[['SCL']] <- SE2.data[['SCL']] / factorR
raster::plot(SE2.data[[3]])


# ------------------------------------------------------------------------------
# 4. Read Shapefile files
# ------------------------------------------------------------------------------
files_shp <- list.files(paths.shape, pattern = "\\.shp$", full.names = TRUE)
print(files_shp)

shape <- sf::read_sf(dsn = files_shp[1])
shape <- st_transform(shape, crs(SE2.data))
names(shape)
unique(shape$ID)   # change "crop" if different


# ------------------------------------------------------------
# Define disease incidence proportion
# Example: 40% of trees affected
# ------------------------------------------------------------
incidence <- 0.40   # change to 0.60 if needed

# Set seed for reproducibility
set.seed(123)

# Total number of trees
n_trees <- nrow(shape)

# Number of diseased trees
n_diseased <- round(n_trees * incidence)

# Create binary disease column: 0 = healthy, 1 = diseased
shape$disease <- 0

# Randomly select diseased trees
idx_diseased <- sample(1:n_trees, size = n_diseased, replace = FALSE)
shape$disease[idx_diseased] <- 1

# ------------------------------------------------------------------------------
# 5. Coordinate system alignment (VERY IMPORTANT)
# ------------------------------------------------------------------------------

# Check CRS of raster and shapefile
crs(SE2.data)
st_crs(shape)

# Reproject shapefile to match Sentinel-2 CRS
shape <- sf::st_transform(shape, crs = sf::st_crs(crs(SE2.data)))

# ------------------------------------------------------------------------------
# 6. Convert to terra format (recommended for raster operations)
# ------------------------------------------------------------------------------

# Convert raster and vector to terra objects
r <- terra::rast(files.SE2[2])
names(r) <- SE_bands

# Apply reflectance scaling
r[[rfl.bands]] <- r[[rfl.bands]] * factorR

# Convert shapefile to terra
v <- terra::vect(shape)

# Ensure both have same CRS
v <- terra::project(v, terra::crs(r))

# ------------------------------------------------------------------------------
# 7. Zoom to shapefile extent
# ------------------------------------------------------------------------------

# Define buffer (in meters if using UTM)
buffer <- 50

# Get extent of shapefile + buffer
e <- terra::ext(v)
e <- terra::ext(
  e[1] - buffer, e[2] + buffer,
  e[3] - buffer, e[4] + buffer
)

# Plot zoomed raster
plot(r[[3]],
     xlim = c(e[1], e[2]),
     ylim = c(e[3], e[4]),
     main = "Zoom Sentinel-2 with shapefile")

# Overlay shapefile
plot(v, add = TRUE, col = ifelse(v$disease == 1, "red", "green"), pch = 16)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#	8. Getting the time series from Stack files      -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data.serie<- getSeries(pathRaster=paths.se2a, shapefile=shape, band_names=SE_bands,factorR=factorR, get.indices = F)

indices <-getIndicesSE2(data.serie[,rfl.bands], sensor = "Sentinel-2a", df.data = NULL, fast.process =T)
head(indices)

data.export <-cbind(data.serie,indices)
head(data.export)

# keep vegetation-related pixels (SCL == 4 (vegetation) / SCL == 5 (Not-vegetated() / SCL == 7 (Unclassified (often usable over orchards / sparse canopy))
data.export.veg <- data.export %>%
  dplyr::filter(SCL %in% c(4, 5, 7))
head(data.export.veg)

file.to.export<-paste(paths.outs,'TimeSerie_SE2_Olive_2022_2025_SCL_filtered_with_deseases.csv',sep='')
write.table(data.export.veg, file = file.to.export, sep=",", row.names = FALSE, col.names = T,append = F)


# ------------------------------------------------------------
# 9. Control whether to use SCL-filtered data or full data
# ------------------------------------------------------------
use_scl_filter <- FALSE   # TRUE = use data.export.veg ; FALSE = use data.export

data.plot <- if (use_scl_filter) data.export.veg else data.export


# Date in correct format
data.plot$Date <- as.Date(data.plot$Date)

# check available dates
unique(data.plot$Date)

# choose one date
date_sel <- max(as.Date(data.plot$Date))
# choose one date
date_sel <- max(as.Date(data.export$Date))

# subset one date
ndvi_date <- data.export %>%
  dplyr::filter(as.Date(Date) == date_sel) %>%
  dplyr::select(ID, NDVI)

head(ndvi_date)

# ------------------------------------------------------------
# 10. Join NDVI to shapefile
# ------------------------------------------------------------
shape_ndvi <- shape %>%
  dplyr::left_join(ndvi_date, by = "ID")

head(shape_ndvi)

# ------------------------------------------------------------
# 11. Plot shapefile coloured by NDVI
# ------------------------------------------------------------
plot(shape_ndvi["NDVI"], main = paste("NDVI on", date_sel))




# ------------------------------------------------------------
# 12. Calculate summary statistics by date
# ------------------------------------------------------------

use_scl_filter <- TRUE   # TRUE = use data.export.veg ; FALSE = use data.export

data.plot <- if (use_scl_filter) data.export.veg else data.export
# Date in correct format
data.plot$Date <- as.Date(data.plot$Date)


ts_summary <- data.plot %>%
  group_by(Date) %>%
  summarise(
    mean_NDVI   = mean(NDVI, na.rm = TRUE),
    median_NDVI = median(NDVI, na.rm = TRUE),
    p05_NDVI    = quantile(NDVI, 0.05, na.rm = TRUE),
    p25_NDVI    = quantile(NDVI, 0.25, na.rm = TRUE),
    p75_NDVI    = quantile(NDVI, 0.75, na.rm = TRUE),
    p95_NDVI    = quantile(NDVI, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

head(ts_summary)


# ------------------------------------------------------------
# 11. Plot all point/polygon time series in light grey
#     and overlay mean/median/percentiles in red
# ------------------------------------------------------------
ggplot() +
  # all individual time series
  geom_line(
    data = data.plot,
    aes(x = Date, y = NDVI, group = ID),
    color = "grey75",
    alpha = 0.5
  ) +
  # percentile band (25% - 75%)
  geom_ribbon(
    data = ts_summary,
    aes(x = Date, ymin = p25_NDVI, ymax = p75_NDVI),
    fill = "red",
    alpha = 0.2
  ) +
  # mean line
  geom_line(
    data = ts_summary,
    aes(x = Date, y = mean_NDVI),
    color = "red",
    linewidth = 1.2
  ) +
  # median line
  geom_line(
    data = ts_summary,
    aes(x = Date, y = median_NDVI),
    color = "darkred",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    title = "NDVI time series for all polygons",
    subtitle = "Grey lines = individual polygons; Red = mean; Shaded area = 25th–75th percentile",
    x = "Date",
    y = "NDVI"
  ) +
  theme_bw()


ts_disease <- data.plot %>%
  group_by(Date, SCL) %>%
  summarise(mean_NDVI = mean(NDVI, na.rm = TRUE),
            .groups = "drop")

ggplot() +
  # all individual curves
  geom_line(data = data.plot,
            aes(x = Date, y = NDVI, group = ID),
            color = "grey85", alpha = 0.3) +

  # mean per disease class
  geom_line(data = ts_disease,
            aes(x = Date, y = mean_NDVI, color = as.factor(SCL)),
            linewidth = 1.3) +

  labs(title = "NDVI dynamics by SCL",
       color = "Disease level") +
  theme_minimal()
