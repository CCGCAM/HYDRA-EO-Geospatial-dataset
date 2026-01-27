###############################################################################
# Multi-resolution synthetic rasters from an AOI polygon (GeoJSON)
# Robust version: avoids focal() errors for very small rasters (e.g., FLEX 300 m)
###############################################################################

library(sf)
library(terra)
library(leaflet)

# ---------------------------------------------------------------------------
# INPUTS
# ---------------------------------------------------------------------------
aoi_path <- "data/sites/Olive_orchard.geojson"
out_dir  <- "data/rasters"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# READ AOI POLYGON
# ---------------------------------------------------------------------------
aoi_sf <- st_read(aoi_path, quiet = TRUE)
aoi_sf <- st_make_valid(aoi_sf)

# IMPORTANT: resolutions are in meters -> AOI must be in a projected CRS (e.g., UTM).
if (st_is_longlat(aoi_sf)) {
  # Change EPSG to your correct UTM zone if needed
  aoi_sf <- st_transform(aoi_sf, 32630)
}

aoi <- vect(aoi_sf)

# ---------------------------------------------------------------------------
# SAFE FOCAL MEAN (prevents errors when raster is smaller than the window)
# ---------------------------------------------------------------------------
safe_focal_mean <- function(r, max_window = 7) {
  nr <- nrow(r)
  nc <- ncol(r)

  # terra constraint: window_rows <= 2*nr and window_cols <= 2*nc
  wmax_r <- min(max_window, 2 * nr - 1)
  wmax_c <- min(max_window, 2 * nc - 1)
  wsize  <- min(wmax_r, wmax_c)

  # Ensure odd and at least 3; otherwise skip smoothing
  if (wsize < 3) return(r)
  if (wsize %% 2 == 0) wsize <- wsize - 1
  if (wsize < 3) return(r)

  focal(r, w = matrix(1, wsize, wsize), fun = mean, na.policy = "omit")
}

# ---------------------------------------------------------------------------
# TEMPLATE CREATION (masked raster at target resolution)
# ---------------------------------------------------------------------------
make_template <- function(aoi_vect, res_m) {
  r <- rast(ext(aoi_vect), resolution = res_m, crs = crs(aoi_vect))
  if (ncell(r) == 0) stop("Template has 0 cells: check CRS/extent/resolution.")

  # give it values so mask() works
  values(r) <- 1

  r <- extend(r, aoi_vect)
  r <- mask(r, aoi_vect)

  if (all(is.na(values(r)))) stop("All cells are NA after masking: polygon/CRS mismatch?")
  r
}

# ---------------------------------------------------------------------------
# SYNTHETIC MULTI-BAND (hyperspectral-like) cube
# ---------------------------------------------------------------------------
make_hyperspec <- function(r_masked, n_bands = 200, seed = 1) {
  set.seed(seed)

  # Base spatial vegetation pattern
  base <- r_masked
  values(base) <- runif(ncell(base), 0.2, 0.9)

  # Safe smoothing (adapts window or skips if raster is too small)
  base <- safe_focal_mean(base, max_window = 7)

  # "Wavelengths" for realism
  wl <- seq(400, 2500, length.out = n_bands)

  # Simple vegetation-like spectral shape (0..~1)
  spectral_shape <- plogis((wl - 700) / 40) * (1 - 0.15 * plogis((wl - 1400) / 80))

  # Stack bands
  cube <- rast(lapply(1:n_bands, function(i) {
    b <- base * spectral_shape[i]
    values(b) <- values(b) + rnorm(ncell(b), sd = 0.01)
    b
  }))
  names(cube) <- sprintf("B%03d", 1:n_bands)
  cube
}

# ---------------------------------------------------------------------------
# SYNTHETIC SINGLE-BAND (thermal / FLEX-like)
# ---------------------------------------------------------------------------
make_singleband <- function(r_masked, seed = 1, lo = 290, hi = 315) {
  set.seed(seed)

  r <- r_masked
  values(r) <- runif(ncell(r), lo, hi)

  # Safe smoothing (adapts window or skips)
  r <- safe_focal_mean(r, max_window = 5)

  r
}

# ---------------------------------------------------------------------------
# SENSOR / RESOLUTION CONFIG
# NOTE: nband > 1 will produce a multi-band raster with that many bands.
# ---------------------------------------------------------------------------
spec <- list(
  UAV_1m     = list(res = 1,   bands = 200),
#  UAV_subm   = list(res = 0.5, bands = 200),
  S2_20m     = list(res = 20,  bands = 12),
  EnMAP_30m  = list(res = 30,  bands = 200),
  PRISMA_30m = list(res = 30,  bands = 200),
  THERM_70m  = list(res = 70,  bands = 1),
  FLEX_300m  = list(res = 300, bands = 1)
)

# ---------------------------------------------------------------------------
# GENERATE AND SAVE RASTERS (robust for small rasters like FLEX 300 m)
# ---------------------------------------------------------------------------
outputs <- list()

for (nm in names(spec)) {

  res_m <- spec[[nm]]$res
  nband <- spec[[nm]]$bands

  message("Processing: ", nm, " | res=", res_m, " m | bands=", nband)

  rmask <- make_template(aoi, res_m)
  message("  template size: ", nrow(rmask), " x ", ncol(rmask),
          " | cells=", ncell(rmask))

  if (nband > 1) {
    ras <- make_hyperspec(rmask, n_bands = nband, seed = 10 + res_m)
    fn  <- file.path(out_dir, paste0(nm, "_", res_m, "m_", nband, "b.tif"))
  } else {
    ras <- make_singleband(rmask, seed = 100 + res_m)
    fn  <- file.path(out_dir, paste0(nm, "_", res_m, "m_1b.tif"))
  }

  writeRaster(ras, fn, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
  outputs[[nm]] <- fn
}

outputs

# ---------------------------------------------------------------------------
# LEAFLET QUICKLOOK (display one representative band per multi-band raster)
# ---------------------------------------------------------------------------

# Leaflet requires lon/lat
aoi_ll <- st_transform(aoi_sf, 4326)

make_quicklook <- function(tif_path, band = 1) {
  r <- rast(tif_path)
  q <- if (nlyr(r) > 1) r[[band]] else r
  project(q, "EPSG:4326")
}

# Choose display band (e.g., band 6 for S2-like 12 bands; band 60 for 200 bands)
q_uav1  <- make_quicklook(outputs$UAV_1m,     band = 60)
q_s2    <- make_quicklook(outputs$S2_20m,     band = 6)
q_enmap <- make_quicklook(outputs$EnMAP_30m,  band = 60)
q_therm <- make_quicklook(outputs$THERM_70m,  band = 1)
q_flex  <- make_quicklook(outputs$FLEX_300m,  band = 1)

pal_for <- function(r) {
  v <- values(r)
  v <- v[is.finite(v)]   # keep only finite numbers (no NA/Inf)

  # If everything is NA, return NULL (we will skip addRasterImage for that layer)
  if (length(v) == 0) return(NULL)

  colorNumeric("viridis", domain = range(v), na.color = "transparent")
}
addRasterSafe <- function(map, r, group_name, opacity = 0.7) {
  pal <- pal_for(r)
  if (is.null(pal)) {
    message("Skipping layer (all NA): ", group_name)
    return(map)
  }
  map |> addRasterImage(r, colors = pal, opacity = opacity, group = group_name)
}
m <- leaflet() |>
  addProviderTiles(providers$Esri.WorldImagery) |>
  addPolygons(data = aoi_ll, fill = FALSE, color = "cyan", weight = 2, group = "AOI")

m <- addRasterSafe(m, q_uav1,  "UAV 1 m")
m <- addRasterSafe(m, q_s2,    "S2 20 m")
m <- addRasterSafe(m, q_enmap, "EnMAP 30 m")
m <- addRasterSafe(m, q_therm, "Thermal 70 m")
m <- addRasterSafe(m, q_flex,  "FLEX 300 m")

m |> addLayersControl(
  overlayGroups = c("AOI", "UAV 1 m", "S2 20 m", "EnMAP 30 m", "Thermal 70 m", "FLEX 300 m"),
  options = layersControlOptions(collapsed = FALSE)
)
