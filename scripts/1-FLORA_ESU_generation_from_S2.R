#
rm(list= ls())
raster_temp_dir<-paste(tempdir(),'raster/',sep='')
files_to_be_removed<-list.files (raster_temp_dir,full.names = T)
file.remove(files_to_be_removed)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#	0. Libraries   -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (!require("ToolsRTM")) { install.packages("ToolsRTM"); require("ToolsRTM") }  ###
if (!require("raster")) { install.packages("raster"); require("raster") }

if (!require("RColorBrewer")) { install.packages("RColorBrewer"); require("RColorBrewer") }  ### colors
if (!require("pls")) { install.packages("pls"); require("pls") }  ### PLSR
if (!require("signal")) { install.packages("signal"); require("signal") }  ### interpolations

if (!require("MASS")) { install.packages("MASS"); require("MASS") }  ## fs
if (!require("caret")) { install.packages("caret"); require("caret") }  ##
if (!require("rgdal")) { install.packages("rgdal"); require("rgdal") }  ##
if (!require("ggplot2")) { install.packages("ggplot2"); require("ggplot2") }  ##
if (!require("sf")) { install.packages("sf"); require("sf") }  ###
if (!require("rgeos")) { install.packages("rgeos"); require("rgeos") }  ###

if (!require("svMisc")) { install.packages("svMisc"); require("svMisc") }  ###
### Shape with polygons
if (!require("terra")) { install.packages("terra"); require("terra") }  ###
if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ###

source('scripts/functions/getSeries.R')
source('scripts/functions/getIndices_SE2.R')
# ======================================================================
# ESU-based heterogeneity mapping from Sentinel-2 (NDRE) at FLEX scale
# - Reads one Sentinel-2 NetCDF (or raster file) with bands + SCL
# - Computes NDRE (B8A, B5)
# - Masks vegetation using SCL == 4
# - Aggregates NDRE to 30 m and 300 m (FLEX) scales
# - Computes within-300 m SD and Cochran sample size (ESUs required)
# - Plots RGB composites (veg only) at 10 m / 30 m / 300 m
# ======================================================================

library(raster)

# -----------------------------
# 0) Sentinel-2 band metadata (optional)
# -----------------------------
Bands   <- c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12')
regions <- c('Coast aerosol','Blue','Green','Red','Red-edge 1','Red-edge 2','Red-edge 3',
             'NIR','Red-edge 4','Water vapour','Cirrus','SWIR1','SWIR2')
wave_2A <- c(442.7,492.4,559.8,664.6,704.1,740.5,782.8,832.8,864.7,945.1,1373.5,1613.7,2202.4)
spatial <- c(60,10,10,10,20,20,20,10,20,60,60,20,20)

db.SE <- data.frame(Band = Bands, Region = regions, Wavelength_nm = wave_2A, Pixel_m = spatial)
print(head(db.SE))

# -----------------------------
# 1) Inputs
# -----------------------------
# NetCDF / raster stack bands expected in this order:
SE_bands <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','SCL')

# Scaling factor used in Sentinel-2 L2A reflectance products
factor_se2 <- 1/10000

# Folder containing your time series (one file per date)
paths.se2a <- 'data/Sentinel-2/Series/NL-Loobos/'

# Which file/date to load (adjust as needed)
i_file <- 75

# -----------------------------
# 2) Read Sentinel-2 stack and scale reflectance
# -----------------------------
files.SE2 <- list.files(paths.se2a, full.names = TRUE)
stopifnot(length(files.SE2) >= i_file)

SE2.data <- raster::stack(files.SE2[i_file])
names(SE2.data) <- SE_bands

# Scale reflectance bands to [0,1]
SE2.data <- SE2.data * factor_se2

# Restore SCL integer classes (undo scaling)
SE2.data[['SCL']] <- SE2.data[['SCL']] / factor_se2

# Quick check
plot(SE2.data[['B4']], main = "Sentinel-2 B4 (Red)")

# Optional: inspect SCL distribution
print(table(values(SE2.data[['SCL']])))

# -----------------------------
# 3) NDRE computation (B8A, B5)
# NDRE = (NIRre - RedEdge) / (NIRre + RedEdge)
# Using Sentinel-2: B8A (NIR narrow) and B5 (Red-edge 1)
# -----------------------------
ndre <- (SE2.data[['B8A']] - SE2.data[['B5']]) / (SE2.data[['B8A']] + SE2.data[['B5']])
names(ndre) <- "NDRE"
plot(ndre, main = "NDRE (unmasked)")

# -----------------------------
# 4) Vegetation mask using SCL == 4
# SCL classes (L2A): 4 vegetation, 5 bare soil, 6 water, 7-10 clouds/shadows, 11 snow
# -----------------------------
veg_mask <- SE2.data[['SCL']] == 4
plot(veg_mask, main = "Vegetation mask (SCL == 4)")

# Apply mask: keep vegetation pixels only
ndre_veg <- mask(ndre, veg_mask, maskvalue = FALSE, updatevalue = NA)
plot(ndre_veg, main = "NDRE (vegetation only)")

# Scene-level stats (vegetation only)
mu_veg    <- cellStats(ndre_veg, stat = "mean", na.rm = TRUE)
sigma_veg <- cellStats(ndre_veg, stat = "sd",   na.rm = TRUE)
cv_veg    <- sigma_veg / abs(mu_veg)
cat("Scene stats (veg only): mean =", mu_veg, "sd =", sigma_veg, "cv =", cv_veg, "\n")

# -----------------------------
# 5) RGB composites (veg only) at 10 m / 30 m / 300 m
# (Using B4,B3,B2 = R,G,B)
# -----------------------------
r10 <- mask(SE2.data[['B4']], veg_mask, maskvalue = FALSE, updatevalue = NA)
g10 <- mask(SE2.data[['B3']], veg_mask, maskvalue = FALSE, updatevalue = NA)
b10 <- mask(SE2.data[['B2']], veg_mask, maskvalue = FALSE, updatevalue = NA)

rgb10 <- stack(r10, g10, b10)

# Aggregate to 30 m (factor 3) and 300 m (factor 30) from 10 m base grid
rgb30  <- stack(
  aggregate(r10, fact = 3,  fun = mean, na.rm = TRUE),
  aggregate(g10, fact = 3,  fun = mean, na.rm = TRUE),
  aggregate(b10, fact = 3,  fun = mean, na.rm = TRUE)
)

rgb300 <- stack(
  aggregate(r10, fact = 30, fun = mean, na.rm = TRUE),
  aggregate(g10, fact = 30, fun = mean, na.rm = TRUE),
  aggregate(b10, fact = 30, fun = mean, na.rm = TRUE)
)

plotRGB(rgb10,  r = 1, g = 2, b = 3, stretch = "lin", main = "RGB Sentinel-2 (10 m, veg only)")
plotRGB(rgb30,  r = 1, g = 2, b = 3, stretch = "lin", main = "RGB Sentinel-2 (30 m, veg only)")
plotRGB(rgb300, r = 1, g = 2, b = 3, stretch = "lin", main = "RGB Sentinel-2 (300 m, veg only / FLEX)")

# -----------------------------
# 6) Local heterogeneity (optional): moving-window SD on NDRE (veg only)
# Example: 11x11 pixels at 10 m ~ 110 m window
# -----------------------------
w <- matrix(1, 11, 11)
ndre_sd_local <- focal(ndre_veg, w = w, fun = sd, na.rm = TRUE, pad = TRUE)
names(ndre_sd_local) <- "NDRE_SD_local"
plot(ndre_sd_local, main = "NDRE local SD (~110 m window, veg only)")

# -----------------------------
# 7) FLEX-scale aggregation (10 m -> 300 m)
# Within each 300 m pixel:
# - mean NDRE (context)
# - SD NDRE (within-pixel variability)
# -----------------------------
ndre_300m_mean <- aggregate(ndre_veg, fact = 30, fun = mean, na.rm = TRUE)
ndre_300m_sd   <- aggregate(ndre_veg, fact = 30, fun = sd,   na.rm = TRUE)

plot(ndre_300m_mean, main = "NDRE mean at 300 m (veg only)")
plot(ndre_300m_sd,   main = "NDRE SD within 300 m (veg only)")

# -----------------------------
# 8) Cochran sample size (ESUs required) per 300 m pixel
# n = (Z * sigma / E)^2
# - Z: Z-score (1.96 for 95% confidence)
# - sigma: within-pixel SD (ndre_300m_sd)
# - E: allowed absolute error in NDRE units (choose tolerance)
# -----------------------------
Z <- 1.96
E_abs <- 0.05  # example tolerance in NDRE units; adjust for your requirements

n_per_flex <- ceiling((Z * ndre_300m_sd / E_abs)^2)

# Clean up pixels where SD could not be computed
n_vals <- values(n_per_flex)
n_vals[is.infinite(n_vals)] <- NA
values(n_per_flex) <- n_vals

# Plot continuous ESU requirement map
plot(n_per_flex, main = "Required ESUs per 300 m FLEX pixel (Cochran, veg only)")

# -----------------------------
# 9) Plot ESU requirement map with discrete classes (more interpretable)
# -----------------------------
brks <- c(0, 2, 5, 7, 10, 20, 40)
cols <- colorRampPalette(c("white","peachpuff","orange","gold","yellowgreen","green"))(length(brks) - 1)

plot(n_per_flex,
     breaks = brks,
     col = cols,
     main = "ESUs required per 300 m FLEX pixel (classified)")

# Optional: print summary
cat("ESUs (scene-level summary):\n")
print(summary(values(n_per_flex)))
