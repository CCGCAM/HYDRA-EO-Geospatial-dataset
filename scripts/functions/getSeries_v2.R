#' Get time series from raster files for points or polygons
#'
#' @param pathRaster Folder with rasters in TIFF format
#' @param shapefile Vector object: sf, SpatVector, or Spatial*
#' @param band_names Names of raster bands
#' @param factorR Scaling factor for reflectance bands
#' @param get.indices Logical. If TRUE, compute Sentinel-2 indices
#' @param fun Polygon extraction function (e.g. mean, median). Ignored for points
#'
#' @return data.frame with extracted values by date
#' @export

getSeries <- function(pathRaster = NULL,
                      shapefile = NULL,
                      band_names = NULL,
                      factorR = NULL,
                      get.indices = TRUE,
                      fun = mean) {

  options(warn = -1)

  # ------------------------------------------------------------
  # 1. Checks
  # ------------------------------------------------------------
  if (is.null(pathRaster)) stop("Please provide 'pathRaster'")
  if (is.null(shapefile)) stop("Please provide 'shapefile'")
  if (is.null(band_names)) stop("Please provide 'band_names'")
  if (is.null(factorR)) stop("Please provide 'factorR'")

  files <- list.files(pathRaster, pattern = "\\.tif$", full.names = TRUE)
  dates <- list.files(pathRaster, pattern = "\\.tif$", full.names = FALSE)
  dates <- extract_dates_from_tiff_files(dates)

  if (length(files) == 0) {
    stop("No TIFF files found in pathRaster")
  }

  # ------------------------------------------------------------
  # 2. Convert shapefile to terra::SpatVector
  # ------------------------------------------------------------
  if (inherits(shapefile, "sf") || inherits(shapefile, "sfc")) {
    v <- terra::vect(shapefile)
  } else if (inherits(shapefile, "SpatVector")) {
    v <- shapefile
  } else if (inherits(shapefile, "Spatial")) {
    v <- terra::vect(shapefile)
  } else {
    stop("shapefile must be an sf, SpatVector, or Spatial object")
  }

  # assign ID if not present
  v_df <- terra::as.data.frame(v)
  if (!("ID" %in% names(v_df))) {
    v$ID <- seq_len(nrow(v_df))
  }

  # geometry type
  geom_type <- unique(tolower(terra::geomtype(v)))
  is_points <- any(grepl("point", geom_type))

  # ------------------------------------------------------------
  # 3. Loop over rasters
  # ------------------------------------------------------------
  message("Processing rasters ...")
  progress_bar <- txtProgressBar(min = 0, max = length(files), style = 3, char = "=")

  list_out <- vector("list", length(files))

  for (k in seq_along(files)) {
    setTxtProgressBar(progress_bar, k)

    # read raster with terra
    r <- terra::rast(files[k])

    if (terra::nlyr(r) != length(band_names)) {
      warning(paste("Skipping file due to unexpected number of bands:", basename(files[k])))
      next
    }

    names(r) <- band_names

    # project vector to raster CRS if needed
    if (!terra::same.crs(r, v)) {
      v_use <- terra::project(v, terra::crs(r))
    } else {
      v_use <- v
    }

    # extract values
    if (is_points) {
      ext <- terra::extract(r, v_use)
    } else {
      ext <- terra::extract(r, v_use, fun = fun, na.rm = TRUE)
    }

    # convert to data.frame
    ext <- as.data.frame(ext)

    # merge vector attributes if needed
    attrs <- terra::as.data.frame(v_use)

    # ext$ID is internal row match from terra::extract
    data_write <- merge(ext, attrs, by = "ID", all.x = TRUE)

    # add date
    data_write$Date <- dates[k]

    # scale reflectance bands only
    band_cols <- intersect(band_names, names(data_write))
    scl_cols  <- intersect(c("SCL", "QA", "mask"), band_cols)
    rfl_cols  <- setdiff(band_cols, scl_cols)

    if (length(rfl_cols) > 0) {
      data_write[, rfl_cols] <- data_write[, rfl_cols, drop = FALSE] * factorR
    }

    # if SCL was incorrectly scaled before, leave untouched here
    # add indices if requested
    if (isTRUE(get.indices)) {
      message(paste0("\nAdding spectral indices for image: ", k, "/", length(files)))

      # compute only if reflectance bands exist
      se2_bands <- intersect(
        c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"),
        names(data_write)
      )

      if (length(se2_bands) >= 4) {
        idx <- ToolsRTM::getIndicesSE2(
          df = data_write[, se2_bands, drop = FALSE],
          sensor = "Sentinel-2a",
          df.data = NULL,
          fast.process = TRUE
        )
        data_write <- cbind(data_write, idx)
      }
    }

    list_out[[k]] <- data_write
  }

  close(progress_bar)

  data_export <- do.call(rbind, list_out[!sapply(list_out, is.null)])
  rownames(data_export) <- NULL

  return(data_export)
}
