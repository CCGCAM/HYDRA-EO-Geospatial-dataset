rm(list=ls())

#::::::::::::::::::::::::::::::::::::::::::::::::::::
#::	0. load main Libraries   -----
#::::::::::::::::::::::::::::::::::::::::::::::::::::

# Load required packages (install if missing)

pkgs <- c("sf", "data.table", "dplyr", "ggplot2", "lubridate", "tidyr", "scales", "zoo")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# Load  ToolsRTM
require(ToolsRTM)

out_dir <- "plots"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Paths / inputs
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Define input data (shapefile + time-series table) and basic parameters

path_shapes <- "data/sites"
shp_name    <- "Orchards"
ts_file     <- "Tables/results/TimeSerie_SE2_inChaparrillo_2022_2025.csv"

# Crops of interest (filter parcels)
crops_keep <- c("Olive", "Pistachio")

# SCL classes to keep (4=vegetation, 5=not-vegetation/soil) - adjust if needed
scl_keep <- c(2,3,4, 5,7)


#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2. Read & prepare shapefile (attributes only)
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read orchard polygons and extract only the attributes needed for joining
shape <- st_read(file.path(path_shapes, paste0(shp_name, ".shp")), quiet = TRUE) |>
  filter(crop %in% crops_keep)

# Keep only the fields we want + rename id -> ID
shape_att <- shape |>
  st_drop_geometry() |>
  transmute(
    ID      = as.character(id),  # rename id -> ID (to match CSV)
    crop    = as.character(crop),
    stress  = as.character(stress),
    country = as.character(country)
  ) |>
  distinct(ID, .keep_all = TRUE)

#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3. Read & prepare time serie table
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read pixel time-series and apply a basic SCL filter

df_raw <- data.table::fread(ts_file, header = TRUE, sep = ",")

df_raw <- df_raw %>%
  mutate(ID = as.integer(ID) + 1L) %>%    # shift 1..8 -> 2..9
  mutate(ID = as.character(ID)) |>
  dplyr::filter(SCL %in% scl_keep)


#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4. Join attributes -> time-series
#::::::::::::::::::::::::::::::::::::::::::::::::::::

# Merge orchard attributes (crop, stress, country) into pixel time-series
# Add temporal variables for later aggregation (year / month / week)

df <- shape_att |>
  left_join(df_raw, by = "ID") %>%
  drop_na() %>%
  mutate(
    Date  = as.Date(Date),
    Year  = year(Date),
    Month = month(Date),
    Week = floor_date(Date, unit = "week", week_start = 1))



#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 5.1 QC: Number of pixels per SCL class and orchard
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check how many valid pixels exist per parcel and month
# Useful to detect clouds, mixed pixels or missing acquisitions


df_scl_count <- df %>%
  mutate(
    Date = Date,
    SCL = as.factor(SCL),
    Year = Year,
    Month = Month,
    Week = Week,
    ID  = as.factor(ID)
  ) %>%
  count(ID, SCL, Month,name = "n_pix")

plot.n <- ggplot(df_scl_count, aes(x = ID, y = n_pix, fill = SCL)) +
  geom_col() +
  facet_wrap(~Month, scales = "free_y") +
  theme_bw() +
  labs(x = "Orchard (ID)", y = "Number of pixels", fill = "SCL") +
  theme(axis.text.x = element_text(face="bold", angle=45, hjust=1),
        axis.text.y = element_text(face="bold"),
        axis.title  = element_text(face="bold"),
        strip.text  = element_text(face="bold"))

print(plot.n)
ggsave(filename = file.path(out_dir, "1-SCL-pixels_per_orchard_month_2022_2025.png"),
  plot = plot.n, width = 12,height = 7,dpi = 300)

#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 5.2 QC: Vegetation vs Soil fraction per plot
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Evaluate the fraction of soil-visible pixels (SCL=5)
# relative to vegetation pixels (SCL=4) per plot and month.
# This helps quantify canopy cover dynamics and soil influence.

df_cover <- df %>%
  group_by(ID, crop, Month) %>%
  summarise(
    veg  = sum(SCL == 4, na.rm = TRUE),
    soil = sum(SCL == 5, na.rm = TRUE),
    total = veg + soil,
    soil_frac = ifelse(total > 0, soil / total, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    ID_crop = paste(ID, crop, sep = "_")
  )

df_cover <- df_cover %>%
  arrange(crop, ID) %>%   # primero ordena por crop
  mutate(
    crop = factor(crop, levels = c("Olive","Pistachio")),  # fuerza orden
    ID_crop = paste(ID, crop, sep = "_")
  ) %>%
  arrange(crop, ID) %>%
  mutate(
    ID_crop = factor(ID_crop, levels = unique(ID_crop))
  )



plot.soil <- ggplot(df_cover, aes(Month, soil_frac, color = ID_crop)) +
  geom_line(size = 1) +
  theme_bw() +
  labs(
    x = "Months",
    y = "Soil fraction",
    color = "Orchard"
  ) +
  scale_color_manual(values = c(
      "5_Olive" = "#1b7837","6_Olive" = "#5aae61","9_Olive" = "#a6dba0",
      "2_Pistachio" = "#2166ac","3_Pistachio" = "#4393c3", "4_Pistachio" = "#000000",
      "7_Pistachio" = "#053061","8_Pistachio" = "#74add1")) +
  theme(
    axis.text = element_text(face="bold"),
    axis.title = element_text(face="bold"),
    legend.title = element_text(face="bold")
  )

print(plot.soil)

ggsave(filename = file.path(out_dir, "1-Soil-fraction_per_orchard_2022_2025.png"),
       plot = plot.soil, width = 12,height = 7,dpi = 300)

#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 5.3 QC: Vegetation vs Soil fraction per orchard
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Evaluate the fraction of soil-visible pixels (SCL=5)
# relative to vegetation pixels (SCL=4) per orchard and month.
# This helps quantify canopy cover dynamics and soil influence.

df_cover_crop <- df %>%
  group_by(crop, Month) %>%
  summarise(
    veg  = sum(SCL == 4),
    soil = sum(SCL == 5),
    soil_frac = soil / (veg + soil),
    .groups = "drop"
  )

plot.soil <- ggplot(df_cover_crop, aes(Month, soil_frac, color = crop, group = crop)) +
  geom_line(size = 1.2) +
  geom_point() +
  theme_bw() +
  labs(x = "Month", y = "Soil fraction", color = "Crop") +
  theme(
    axis.text = element_text(face="bold"),
    axis.title = element_text(face="bold"),
    legend.title = element_text(face="bold")
  )

print(plot.soil)

ggsave(filename = file.path(out_dir, "1-Soil-fraction_per_crop_2022_2025.png"),
       plot = plot.soil, width = 12,height = 7,dpi = 300)

#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 6. Weekly orchard-level time series
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Reduce pixel noise by aggregating values per parcel and week
# Median used for robustness against mixed pixels



df.week <- df %>%
  group_by(SCL, ID, crop, Week) %>%
  summarise(
    CR.red.nir.6 = median(CR.red.nir.6, na.rm = TRUE),
    CR.SWIR = median(CR.SWIR, na.rm = TRUE),
    NDVI = median(NDVI, na.rm = TRUE),
    n_pix = n(),
    .groups = "drop"
  )   %>%
  mutate(SCL_lab = recode_factor(
    as.factor(SCL),
    `4` = "Covered vegetation",
    `5` = "not vegetation"
  )) |>
  dplyr::filter(SCL %in% c(4,5))  %>%
  filter(n_pix > 15)


#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 6.2 Smooth the time serie
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Reduce pixel noise by aggregating values per parcel and week
# Median used for robustness against mixed pixels



df.week.smooth <- df.week %>%
  arrange(ID, Week) %>%
  group_by(SCL, ID, crop) %>%
  mutate(
    CR.red.nir = rollmedian(CR.red.nir.6, k = 5, fill = NA, align = "center"),
    CR.SWIR = rollmedian(CR.SWIR, k = 5, fill = NA, align = "center"),
    NDVI = rollmedian(NDVI, k = 5, fill = NA, align = "center")
  ) %>%
  ungroup()  %>%
  filter(ID != 9)



plot.credge <-ggplot(df.week.smooth,
       aes(Week, CR.red.nir, color = ID)) +
  geom_line(linewidth = 1) +
  facet_wrap(~SCL_lab+crop) +
  theme_bw()
print(plot.credge)

ggsave(filename = file.path(out_dir, "2-Weekly_CR.rednir_byID_SCL_20222_2025.png"),
       plot = plot.credge, width = 12, height = 6, dpi = 300)



plot.swir <- ggplot(df.week.smooth,
       aes(Week, CR.SWIR, color = ID)) +
  geom_line(linewidth = 1) +
  facet_wrap(~SCL_lab+crop) +
  theme_bw()
print(plot.swir)

ggsave(filename = file.path(out_dir, "2-Weekly_CR.rSWIR_byID_SCL_20222_2025.png"),
       plot = plot.swir, width = 12, height = 6, dpi = 300)



#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 7. Monthly distribution per parcel (boxplot diagnostic)
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Visualise monthly distribution of pixel values per parcel (helps identify soil influence)

df.plot <- df |>
  dplyr::filter(SCL %in% c(4,5))


plot.box <- ggplot(df.plot, aes(x = Month, y = CR.red.nir.6, fill = ID)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~SCL, scales = "free_y") +
  theme_bw() +
  labs(x = "Month", y = "CR-red-nir") +
  theme(
    legend.position = "right",
    strip.text = element_text(face="bold"),
    axis.text = element_text(face="bold"),
    axis.title = element_text(face="bold")
  )

print(plot.box)

ggsave(filename = file.path(out_dir, "3-Boxplot_CR.rednir_month_byID_SC_2022_2025.png"),
  plot = plot.box, width = 12, height = 6, dpi = 300)



#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 8. Monthly time series per parcel (median ± IQR)
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Build robust monthly summaries per parcel (median and interquartile range)


df.ts <- df.plot |>
  group_by(SCL, ID, Year,Month) |>
  summarise(
    across(
      where(is.numeric) & !any_of(c("Year","Month","SCL")),   # avoid resummarising time/SCL
      list(
        med = ~ median(.x, na.rm = TRUE),
        q25 = ~ quantile(.x, 0.25, na.rm = TRUE, names = FALSE),
        q75 = ~ quantile(.x, 0.75, na.rm = TRUE, names = FALSE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(Date_m = as.Date(sprintf("%d-%02d-15", Year, Month)))

plot.ts <- ggplot(df.ts, aes(Date_m, CR.red.nir.6_med, color = factor(ID), group = ID)) +
  geom_ribbon(aes(ymin = CR.red.nir.6_q25, ymax = CR.red.nir.6_q75, fill = factor(ID)),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~SCL, scales = "free_y") +
  theme_bw() +
  labs(x = "", y = "CR-red-nir (median ± IQR)") +
  theme(legend.position="right",
        legend.title = element_blank(),
        strip.text = element_text(face="bold"),
        axis.text = element_text(face="bold"),
        axis.title = element_text(face="bold"))

print(plot.ts)

ggsave(
  filename = file.path(out_dir, "4-Monthly_CRrednir_median_IQR_byID_SCL_2022_2025.png"),
  plot = plot.ts, width = 12, height = 6, dpi = 300)



#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 9. Crop-level CR.SWIR signal (orchard mean -> crop mean)
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Optional: focus on vegetation pixels and compare crop behaviour over time


df.clean <- df |>
  dplyr::filter(SCL %in% c(4,5))  # change to <= if you intentionally want mixed/soil

# One time series per orchard (daily mean across pixels)
df.orchard <- df.clean |>
  group_by(ID, crop, Date) |>
  summarise(
    CR.SWIR = mean(CR.SWIR, na.rm = TRUE),
    .groups = "drop"
  )

# Crop-level signal (mean ± SE across orchards)
df.crop <- df.orchard %>%
  group_by(crop, Date) %>%
  summarise(
    CR.SWIR = mean(CR.SWIR, na.rm = TRUE),
    sd = sd(CR.SWIR, na.rm = TRUE),
    n = dplyr::n(),
    se = dplyr::if_else(n > 1, sd / sqrt(n), 0),
    .groups = "drop"
  )

plot.crop <- ggplot(df.crop, aes(Date, CR.SWIR, color = crop)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    data = df.crop %>% filter(n > 1),
    aes(ymin = CR.SWIR - se, ymax = CR.SWIR + se, fill = crop),
    alpha = 0.2, colour = NA
  ) +
  theme_bw() +
  labs(x = "", y = "CR-SWIR") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(face="bold"),
    axis.title = element_text(face="bold")
  )

print(plot.crop)

ggsave(filename = file.path(out_dir, "5-CR.SWIR_CropMean.png"),
  plot = plot.crop, width = 12, height = 6, dpi = 300)





#::::::::::::::::::::::::::::::::::::::::::::::::::::
# 9. Anomalies
#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Compute simple anomalies and standardized anomalies (Z-score)
# relative to a "healthy" multi-year baseline.
#
# Rationale:
# - Perennial orchards show strong seasonal cycles.
# - A multi-year baseline reduces interannual noise.
# - Weekly climatology (week-of-year) avoids mixing seasons.
#
# Recommended default:
# - Use vegetation pixels only (SCL=4)
# - Build baseline using first 3 years (e.g., 2022–2024)
# - Evaluate anomalies in all years (2022–2026)

#-----------------------------
# 9.1 Inputs / baseline setup
#-----------------------------

baseline_years <- c(2022, 2023, 2024)   # 3-year baseline (adjust if needed)
use_scl <- 5                          # focus on vegetation pixels
var_target <-   "CR.red.nir.6"#"CR.SWIR"                 # change to "CR.red.nir.6" if desired

#-----------------------------
# 9.2 Prepare weekly series
#-----------------------------
# Add Year / Month and Week-of-year index (ISO week)
# Important: week-of-year is used to build a seasonal climatology

df.an <- df.week.smooth %>%
  filter(SCL == use_scl) %>%
  mutate(
    Year     = year(Week),
    Month    = month(Week),
    Week_num = isoweek(Week)
  )


df.an <- df.an %>%
  filter((Month %in% c(3,4,5,6,7,8,9)))
#-----------------------------
# 9.3 Baseline climatology (per crop, per week-of-year)
#-----------------------------
# Compute baseline mean and sd for each crop and week_num.
# This creates an expected seasonal curve from baseline years.

df.base.stats <- df.an %>%
  filter(Year %in% baseline_years) %>%
  group_by(crop, Week_num) %>%
  summarise(
    base_mean = mean(.data[[var_target]], na.rm = TRUE),
    base_sd   = sd(.data[[var_target]],  na.rm = TRUE),
    n_base    = sum(!is.na(.data[[var_target]])),
    .groups   = "drop"
  ) %>%
  # Avoid division-by-zero: if sd is 0 or NA, set to NA
  mutate(base_sd = ifelse(is.na(base_sd) | base_sd == 0, NA_real_, base_sd))

#-----------------------------
# 9.4 Compute anomalies for all years
#-----------------------------
# - Simple anomaly: value - baseline mean
# - Z-score anomaly: (value - baseline mean) / baseline sd

df.anom <- df.an %>%
  left_join(df.base.stats, by = c("crop", "Week_num")) %>%
  mutate(
    anom_simple = .data[[var_target]] - base_mean,
    anom_z      = ( .data[[var_target]] - base_mean ) / base_sd
  )

#-----------------------------
# 9.5 Optional QC filters
#-----------------------------
# Remove weeks where baseline is not reliable:
# - Too few baseline samples
# - Missing sd (cannot compute Z-score)

min_base_n <- 3  # require at least N baseline points per (crop, Week_num)

df.anom.min <- df.anom %>%
  filter(n_base >= min_base_n)

#-----------------------------
# 9.6 Plot: Z-score anomalies per crop (all orchards)
#-----------------------------
# Use facet by crop, color by orchard ID.
# Add reference lines at 0 and +/-1 standard deviation.



plot.anom.z <- ggplot(df.anom.min, aes(x = Week, y = anom_z, color = as.character(ID), group = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~crop, scales = "free_x") +
  theme_bw() +
  labs(
    x = "",
    y = paste0("Z-score anomaly (", var_target, ")"),
    color = "Orchard ID"
  ) +
  theme(
    axis.title  = element_text(face="bold"),
    axis.text   = element_text(face="bold"),
    strip.text  = element_text(face="bold"),
    legend.title = element_text(face="bold")
  )

print(plot.anom.z)
ggsave(filename = file.path(out_dir, paste("9-Anomaly_Z-withSCL-",use_scl,"-",var_target,".png",sep='')),
       plot = plot.anom.z, width = 12, height = 6, dpi = 300)


#-----------------------------
# 9.7 Plot: simple anomalies (optional)
#-----------------------------

plot.anom.simple <- ggplot(df.anom, aes(x = Week, y = anom_simple, color = as.character(ID), group = ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~crop, scales = "free_x") +
  theme_bw() +
  labs(
    x = "",
    y = paste0("Simple anomaly (", var_target, " - baseline mean)"),
    color = "Orchard ID"
  ) +
  theme(
    axis.title  = element_text(face="bold"),
    axis.text   = element_text(face="bold"),
    strip.text  = element_text(face="bold"),
    legend.title = element_text(face="bold")
  )

print(plot.anom.simple)

ggsave(filename = file.path(out_dir, paste("9-Anomaly_Simple_withSCL-",use_scl,"-",var_target,".png",sep='')),
       plot = plot.anom.simple, width = 12, height = 6, dpi = 300)










#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3. Get a look-up table (LUT)  -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

nSamples = 500
inputs.prosail = ToolsRTM::inputsPROSAIL

# Step 1: Define the column names and their initial types
column_names <- c("N", "Cab", "Car", "Anth", "Cbrown", "alpha",
                  "EWT", "LMA", "psoil", "LAI", "TypeLidf" ,"LIDFa", "LIDFb",
                  "hspot", "tts", "tto", "psi")
# Step 2: Create an empty Data Frame with the specified columns
df.LUT <- data.frame(matrix(ncol = length(column_names), nrow = nSamples))
colnames(df.LUT) <- column_names
#set.seed('12345')
# We can use the clock of the computer to set the seed to a random seed every time
#set.seed(Sys.time())
## Adjust PROSPECT-PRO parameters
df.LUT$N <- stats::runif(nSamples,min = 1.0,max=3.0)
df.LUT$Cab <- stats::runif(nSamples,min = 10,max=75)#stats::rnorm(nSamples, mean = 20, sd = 5)
df.LUT$Car = stats::runif(nSamples,min = 0,max=5)
df.LUT$Anth = stats::runif(nSamples,min = 0,max=0.1)
df.LUT$Cbrown = stats::runif(nSamples,min = 0.0,max=1)
df.LUT$alpha = 40 ## internal parameter in PROSPECT
df.LUT$EWT = stats::runif(nSamples,min = 0.0005,max=0.3)
# In case of PROSPECT-PRO: LMA = Prot + CBC
df.LUT$LMA = stats::runif(nSamples,min = 0.0005,max=0.3)
df.LUT$Prot = stats::runif(nSamples,min = 0.0009,max=0.01)
df.LUT$CBC = df.LUT$LMA -df.LUT$Prot

## Adjust FourSAIL parameters
df.LUT$TypeLidf <- 2
df.LUT$LIDFa <- stats::runif(nSamples,min = 30,max= 70)
df.LUT$LIDFb <-  0 #stats::runif(nSamples,min = -0.45,max= 0.1)
df.LUT$LAI = stats::runif(nSamples,min = 0.5,max=4.0)
df.LUT$hspot = stats::runif(nSamples,min = 0.0,max=1)
df.LUT$psoil <-  stats::runif(nSamples,min = 0.0,max=0.95)
## Viewing angles
df.LUT$tts <-  stats::runif(nSamples,min =25,max=30)
df.LUT$tto <- stats::runif(nSamples,min = 25,max=30)
df.LUT$psi <- stats::runif(nSamples,min = 150,max=160)
head(df.LUT)



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4. Get simulations using PROSAIL-model based on our LUT  -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## 4.1 Get simulations at 1 nm using paralell  -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

no_cores <- parallel::detectCores() - 2
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
sim.rfl<-list()

sims<-foreach::foreach(i=1:nrow(df.LUT)) %dopar% {
  data <- ToolsRTM::dataSpec_PDB
  Rsoil.dry  <- data[,11]  # rsoil1 = dry soil
  Rsoil.wet <- data[,12]   #
  psoil = df.LUT[i,'psoil']
  rsoil_<- c(psoil*Rsoil.dry+(1-psoil)*Rsoil.wet)


  data.foursail<- ToolsRTM::foursail(inputLUT=df.LUT[i,],rsoil=rsoil_,LeafModel = 'PROSPECT-D')
  rdot<-data.foursail[[1]]
  rsot<-data.foursail[[2]]

  rfl.prosail<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=df.LUT[i,'tts'],data.light=ToolsRTM::dataSpec_PDB)
  sim.rfl[[i]]<-rfl.prosail

} ##end paralle

# Close the cluster after parallel processing is complete
parallel::stopCluster(cl)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4.2. Convolved the simulated reflectance at Sentinel-2a sensor  -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

wavelength <- dataSpec_PDB$wavelength
sim.canopy<-do.call(rbind,sims)
colnames(sim.canopy)=paste0("X",wavelength)

df.convolved=get.spectra.convolved(rfl=sim.canopy, sensor='Sentinel2a',plot.spectra = T)
colnames(df.convolved) <- c('ID','B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12')
head(df.convolved)
SE2.bands = c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12')

indices <-getIndicesSE2(df.convolved[,SE2.bands], sensor = "Sentinel-2a", df.data = NULL, fast.process =T)

indices.bands = c('GM1','NDVI','TCARI_OSAVI','CR.SWIR','CR.red.nir')

df.convolved.indices<- cbind(df.convolved,indices[,indices.bands])
colnames(df.convolved.indices) <- c(names(df.convolved),indices.bands)
head(df.convolved.indices)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 4.3. Find the closest one profile  -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


SE2.bands = c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12')
selected_bands = c('B6','B7','B8','NDVI','TCARI_OSAVI','CR.SWIR','CR.red.nir')

matrix.sim <- as.matrix(df.convolved.indices[,selected_bands])
matrix.SE2 <-as.matrix((df.smooth[,selected_bands]))

# View the first rows of the normalized matrices
head(matrix.sim)
head(matrix.SE2)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## 4.3.1 Make the inversion using a merit-RMSE -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
m.inversion <-ToolsRTM::get.inversionOpt(rfl.sensor=matrix.SE2, rfl.rtm=matrix.sim, LUT=df.LUT,
                                         wave=selected_bands,method='merit-NRMSE',nOpt=50, custom_stat=NULL)

head(m.inversion[[1]])
head(m.inversion[[2]])

plant.traits <- m.inversion[[2]]
head(plant.traits)

df.plot <- cbind(df.smooth,plant.traits)

#::::::::::::::::::::::::::::::::::::::::::::::::::::
#	5. Time Series   -----
#::::::::::::::::::::::::::::::::::::::::::::::::::::

#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Temporal series NDVI
#::::::::::::::::::::::::::::::::::::::::::::::::::::


axis_y<-  expression(bold('NDVI'))

plot.ndvi <-ggplot(df.plot,aes(x=Date,y=NDVI, color=Area)) +
  geom_rect(aes(xmin = as.Date('2020-06-01'), xmax = as.Date('2020-08-31'), ymin = -Inf, ymax = Inf),
            fill = "grey",alpha = 0.2)  +
  geom_point(size = 1) +
  geom_line(size = 1) +
  labs(title ='', x='',y=axis_y, size=10,face="bold") +
  theme_bw() +
  scale_color_manual(values = c('forestgreen',"darkgoldenrod4")) +
  theme(legend.position="bottom",
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill="grey97"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.text.x=element_text(angle=0, size=12,face="bold")) +
  scale_x_date(limits = c(as.Date("2017-05-01"), as.Date("2021-01-31")),date_breaks = "12 months", date_labels = "%b %Y") +
  geom_point(size = 1, alpha = 0.3) +
  geom_line(size = 1, alpha = 0.3)
print(plot.ndvi)


#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Temporal series CR-SWIR
#::::::::::::::::::::::::::::::::::::::::::::::::::::

axis_y<-  expression(bold('CR-SWIR'))

plot.cr.swir <-ggplot(df.plot,aes(x=Date,y=CR.SWIR, color=Area)) +
  geom_rect(aes(xmin = as.Date('2020-06-01'), xmax = as.Date('2020-08-31'), ymin = -Inf, ymax = Inf),
            fill = "grey",alpha = 0.2)  +
  geom_point(size = 1) +
  geom_line(size = 1) +
  labs(title ='', x='',y=axis_y, size=10,face="bold") +
  theme_bw() +
  scale_color_manual(values = c('forestgreen',"darkgoldenrod4")) +
  theme(legend.position="bottom",
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill="grey97"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.text.x=element_text(angle=0, size=12,face="bold")) +
  scale_x_date(limits = c(as.Date("2017-05-01"), as.Date("2021-01-31")),date_breaks = "12 months", date_labels = "%b %Y") +
  geom_point(size = 1, alpha = 0.3) +
  geom_line(size = 1, alpha = 0.3)
print(plot.cr.swir)


#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Temporal series CR-Red-edge
#::::::::::::::::::::::::::::::::::::::::::::::::::::

axis_y<-  expression(bold('CR-redge.NIR'))


plot.cr.red <-ggplot(df.plot,aes(x=Date,y=CR.red.nir, color=Area)) +
  geom_rect(aes(xmin = as.Date('2020-06-01'), xmax = as.Date('2020-08-31'), ymin = -Inf, ymax = Inf),
            fill = "grey",alpha = 0.2)  +
  geom_point(size = 1) +
  geom_line(size = 1) +
  labs(title ='', x='',y=axis_y, size=10,face="bold") +
  theme_bw() +
  scale_color_manual(values = c('forestgreen',"darkgoldenrod4")) +
  theme(legend.position="bottom",
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill="grey97"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.text.x=element_text(angle=0, size=12,face="bold")) +
  scale_x_date(limits = c(as.Date("2017-05-01"), as.Date("2021-01-31")),date_breaks = "12 months", date_labels = "%b %Y") +
  geom_point(size = 1, alpha = 0.3) +
  geom_line(size = 1, alpha = 0.3)
print(plot.cr.red)


#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Temporal series predicted Cab
#::::::::::::::::::::::::::::::::::::::::::::::::::::

axis_y <- expression(bold("Pred. Chl-ab (" ~ µg ~ cm^ -2 ~")"))


plot.cab <-ggplot(df.plot,aes(x=Date,y=Cab, color=Area)) +
  geom_rect(aes(xmin = as.Date('2020-06-01'), xmax = as.Date('2020-08-31'), ymin = -Inf, ymax = Inf),
            fill = "grey",alpha = 0.2)  +
  geom_point(size = 1) +
  geom_line(size = 1) +
  labs(title ='', x='',y=axis_y, size=10,face="bold") +
  theme_bw() +
  scale_color_manual(values = c('forestgreen',"darkgoldenrod4")) +
  theme(legend.position="bottom",
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill="grey97"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.text.x=element_text(angle=0, size=12,face="bold")) +
  scale_x_date(limits = c(as.Date("2017-05-01"), as.Date("2021-01-31")),date_breaks = "12 months", date_labels = "%b %Y") +
  geom_point(size = 1, alpha = 0.3) +
  geom_line(size = 1, alpha = 0.3)
print(plot.cab)


#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Temporal series predicted EWT
#::::::::::::::::::::::::::::::::::::::::::::::::::::

axis_y <- expression(bold("Pred. ETW (" ~ g ~ cm^ -2 ~")"))


plot.ewt <-ggplot(df.plot,aes(x=Date,y=EWT, color=Area)) +
  geom_rect(aes(xmin = as.Date('2020-06-01'), xmax = as.Date('2020-08-31'), ymin = -Inf, ymax = Inf),
            fill = "grey",alpha = 0.2)  +
  geom_point(size = 1) +
  geom_line(size = 1) +
  labs(title ='', x='',y=axis_y, size=10,face="bold") +
  theme_bw() +
  scale_color_manual(values = c('forestgreen',"darkgoldenrod4")) +
  theme(legend.position="bottom",
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill="grey97"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.text.x=element_text(angle=0, size=12,face="bold")) +
  scale_x_date(limits = c(as.Date("2017-05-01"), as.Date("2021-01-31")),date_breaks = "12 months", date_labels = "%b %Y") +
  geom_point(size = 1, alpha = 0.3) +
  geom_line(size = 1, alpha = 0.3)
print(plot.ewt)


#::::::::::::::::::::::::::::::::::::::::::::::::::::
# Temporal series predicted LAI
#::::::::::::::::::::::::::::::::::::::::::::::::::::

axis_y <- expression(bold("Pred. LAI (" ~ m^ 2 ~ m^ -2 ~")"))


plot.lai <-ggplot(df.plot,aes(x=Date,y=LAI, color=Area)) +
  geom_rect(aes(xmin = as.Date('2020-06-01'), xmax = as.Date('2020-08-31'), ymin = -Inf, ymax = Inf),
            fill = "grey",alpha = 0.2)  +
  geom_point(size = 1) +
  geom_line(size = 1) +
  labs(title ='', x='',y=axis_y, size=10,face="bold") +
  theme_bw() +
  scale_color_manual(values = c('forestgreen',"darkgoldenrod4")) +
  theme(legend.position="bottom",
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill="grey97"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.text.x=element_text(angle=0, size=12,face="bold")) +
  scale_x_date(limits = c(as.Date("2017-05-01"), as.Date("2021-01-31")),date_breaks = "12 months", date_labels = "%b %Y") +
  geom_point(size = 1, alpha = 0.3) +
  geom_line(size = 1, alpha = 0.3)
print(plot.lai)


# Combine the two plots using cowplot
combined.plot <- cowplot::plot_grid(plot.ndvi + theme(legend.position="none"),
                                    plot.lai +  theme(legend.position="none"),
                                    plot.cr.red + theme(legend.position="none"),
                                    plot.cab +  theme(legend.position="none"),
                                    plot.cr.swir +  theme(legend.position="none"),

                                    plot.ewt +  theme(legend.position="none"),


                                    labels = c('A', 'B','C','D','E','F'), ncol = 2, align = "l", axis = "tb")

print(combined.plot)
## for plots (comparison with Field data)
paths.plots = 'Plots/Series/'
ifelse(!dir.exists(paths.plots), dir.create(paths.plots), FALSE)

ggsave(paste(paths.plots,'1-TimeSerie-SR2.png',sep=''),
       plot=combined.plot,
       width = 40, height = 20,  dpi = 300,units = "cm")




#::::::::::::::::::::::::::::::::::::::::::::::::::::
#	6. Reflectance by Seasons and Years  2019 / 2020  -----
#::::::::::::::::::::::::::::::::::::::::::::::::::::

rfl.bands<-c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12') ##all bands
df.sb <- df.plot

Year_ <- 2020
df.previous <- subset(df.sb, Year <= 2019 & (Month == 6 | Month == 7 | Month ==8 ))
df.2020 <- subset(df.sb, Year == 2020 & (Month == 6 | Month == 7 | Month ==8 ))

df.scatter <-rbind(df.previous,df.2020)
head(df.scatter)

#::::::::::::::::::::::::::::::::::::::::::::::::::::
#   ScatterPLot between NDVI vs LAI
#::::::::::::::::::::::::::::::::::::::::::::::::::::


axis_y <- expression(bold("Pred. LAI (" ~ m^ 2 ~ m^ -2 ~")"))
axis_x <- expression(bold("NDVI"))

scatter.lai <-ggplot(df.scatter,aes(x=NDVI,y=LAI, color=as.factor(Period))) +
  geom_point(size = 3, alpha = 0.3)  +
  ggpubr::stat_cor(method = "pearson", aes(),size=4)  +
  geom_smooth(method = 'lm', formula = y ~ x, linetype = 2,se=F) +
  labs(title = '', y=axis_y,x=axis_x, size=10,face="bold") +
  scale_color_manual(values = c('forestgreen','darkgoldenrod4'),
                     labels = c("Healthy", "Infected")) + theme_bw() +
  theme(legend.position="none",
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill="grey97"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.text.x=element_text(angle=0, size=12,face="bold"))

scatter.lai

# Change marginal plot type
scatter.lai.bbox<- ggExtra::ggMarginal(scatter.lai, type = "boxplot",groupColour = T)
print(scatter.lai.bbox)



#::::::::::::::::::::::::::::::::::::::::::::::::::::
#   ScatterPLot between CR.red.nir vs Cab
#::::::::::::::::::::::::::::::::::::::::::::::::::::


axis_y <- expression(bold("Pred. Chl-ab (" ~ µg ~ cm^ -2 ~")"))
axis_x <- expression(bold("CR.red.nir"))

scatter.cab <-ggplot(df.scatter,aes(x=CR.red.nir,y=Cab, color=as.factor(Period))) +
  geom_point(size = 3, alpha = 0.3)  +
  ggpubr::stat_cor(method = "pearson", aes(),size=4)  +
  geom_smooth(method = 'lm', formula = y ~ x, linetype = 2,se=F) +
  labs(title = '', y=axis_y,x=axis_x, size=10,face="bold") +
  scale_color_manual(values = c('forestgreen','darkgoldenrod4'),
                     labels = c("Healthy", "Infected")) + theme_bw() +
  theme(legend.position="none",
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill="grey97"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.text.x=element_text(angle=0, size=12,face="bold"))

print(scatter.cab)

# Change marginal plot type
scatter.cab.bbox<- ggExtra::ggMarginal(scatter.cab, type = "boxplot",groupColour = T)
print(scatter.cab.bbox)



#::::::::::::::::::::::::::::::::::::::::::::::::::::
#   ScatterPLot between EWT vs CR.SWIR
#::::::::::::::::::::::::::::::::::::::::::::::::::::


axis_y <- expression(bold("Pred. EWT (" ~ g ~ cm^ -2 ~")"))
axis_x <- expression(bold("CR-SWIR"))

scatter.ewt <-ggplot(df.scatter,aes(x=CR.SWIR,y=EWT, color=as.factor(Period))) +
  geom_point(size = 3, alpha = 0.3)  +
  ggpubr::stat_cor(method = "pearson", aes(),size=4)  +
  geom_smooth(method = 'lm', formula = y ~ x, linetype = 2,se=F) +
  labs(title = '', y=axis_y,x=axis_x, size=10,face="bold") +
  scale_color_manual(values = c('forestgreen','darkgoldenrod4'),
                     labels = c("Healthy", "Infected")) + theme_bw() +
  theme(legend.position="none",
        plot.background = element_rect(fill = 'white', colour = 'white'),
        panel.background = element_rect(fill="grey97"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        panel.grid.minor = element_line(colour = "grey97"),
        panel.grid.major = element_line(colour = "grey97"),
        axis.text.x=element_text(angle=0, size=12,face="bold"))

print(scatter.ewt)

# Change marginal plot type
scatter.ewt.bbox<- ggExtra::ggMarginal(scatter.ewt, type = "boxplot",groupColour = T)
print(scatter.ewt.bbox)


# Combine the two plots using cowplot
### Here remove axis-X labels to have more space
combined.bb <- cowplot::plot_grid(scatter.lai.bbox,
                                  scatter.cab.bbox,
                                  scatter.ewt.bbox,

                                  labels = c('A', 'B','C'), ncol = 3, align = "v", axis = "tb",
                                  #   rel_heights = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),  # Set relative heights for each row
                                  rel_widths = rep(2, 2))  # Set relative widths for each column

ggsave(paste(paths.plots,'2-Scattersplots-SR2.png',sep=''),
       plot=combined.bb,
       width = 36, height = 12,  dpi = 300,units = "cm")



#::::::::::::::::::::::::::::::::::::::::::::::::::::
#	7. ANOVA  2019 / 2020  -----
#::::::::::::::::::::::::::::::::::::::::::::::::::::

# Find the indices of columns containing "NDVI"

variables_ <- c('NDVI',"CR.red.nir","CR.SWIR",'Cab','EWT','LAI','Cbrown' )
print(variables_)



for (var in c(variables_)){
  print(var)
  fmla <- as.formula(paste(var," ~ ", paste('Area','ID', sep= " + ")))

  anova_ <- summary(aov(fmla, data = df.sb))

  fmla <- as.formula(paste(var," ~ ", paste('ID', collapse= "+")))
  anova.2020 <- summary(aov(fmla, data = df.2020))
  anova.previous <- summary(aov(fmla, data = df.previous))
  print('Anova for 2020 ')
  print(anova.2020)
  print('Anova for previous years ')
  print(anova.previous)
}


