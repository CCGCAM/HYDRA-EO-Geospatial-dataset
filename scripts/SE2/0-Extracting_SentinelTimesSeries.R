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

## this is the band from SE-2
Bands<-c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12')
regions<-c('Coast aerosol','Blue','Green','Red', 'REd-edge1', 'REd-edge2', 'REd-edge3', 'NIR', 'REd-edge4', 'Water-Vapour', 'SWIR-Cirrus', 'SWIR1','SWIR2')
wave_2A<-c(442.7,492.4,559.8,664.6,704.1,740.5,782.8,832.8,864.7,945.1,1373.5,1613.7,2202.4)
spatial <-c(60,10,10,10,20,20,20,10,20,60,60,20,20)
db.SE<-cbind(Bands,regions,wave_2A)
head(db.SE)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#	1.Get the Path for Sentinel/2 images   -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Bands for each NetCDF
SE_bands<-c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','SCL') ##all bands
rfl.bands<-c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12') ##all bands
factor_se2 = 1/10000

#inputs

#paths.se2a='Sentinel-2/BarkBeetle/2016-2021/' # Scenario for Bark Beetle outbreaks
paths.se2a='data/Sentinel-2/Series/Chaparrillo/' # Scenario for fungus outbreaks

files.SE2 <- list.files(paths.se2a, full.names = T) ###Stack SE2

SE2.data <-raster::stack(files.SE2[2]) ###Stack SE2
names(SE2.data) <-SE_bands
SE2.data <-SE2.data * factor_se2
SE2.data[['SCL']] <- SE2.data[['SCL']] / factor_se2
plot(SE2.data[[3]])
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#	2. Extract series from Stack files      -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


path_shapes <- paste('data/sites',sep='')  # Scenario for fungus outbreaks
shapes <- c('Orchards')


shape <- sf::read_sf(dsn = paste(path_shapes,'/',shapes[1],'.shp',sep = ''))
shape <- st_transform(shape, crs(SE2.data))
names(shape)
unique(shape$crop)   # change "crop" if different

# filter crops
shape <- shape %>%
  filter(crop %in% c("Olive", "Pistachio"))

plot(shape)


## convert to Spatial Dataframe
shape.prj = as(shape, "Spatial")

factorR = 1/10000

# output folder
paths.outs='Tables/results/'
ifelse(!dir.exists(paths.outs), dir.create(paths.outs), FALSE)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#	3. Getting the time series from Stack files      -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data.serie<- getSeries(pathRaster=paths.se2a, shapefile=shape.prj, band_names=SE_bands,factorR=factorR, get.indices = F)

indices <-getIndicesSE2(data.serie[,rfl.bands], sensor = "Sentinel-2a", df.data = NULL, fast.process =T)
head(indices)
data.export <-cbind(data.serie,indices)
head(data.export)
data.export$SCL <- data.export$SCL / factorR


file.to.export<-paste(paths.outs,'TimeSerie_SE2_inChaparrillo_2022_2025.csv',sep='')
write.table(data.export, file = file.to.export, sep=",", row.names = FALSE, col.names = T,append = F)


# keep vegetation-related pixels (SCL == 4 (vegetation) / SCL == 5 (Not-vegetated() / SCL == 7 (Unclassified (often usable over orchards / sparse canopy))
data.export.veg <- data.export %>%
  dplyr::filter(SCL %in% c(4, 5, 7))
head(data.export)


file.to.export<-paste(paths.outs,'TimeSerie_SE2_inChaparrillo_2022_2025_SCL_filtered.csv',sep='')
write.table(data.export.veg, file = file.to.export, sep=",", row.names = FALSE, col.names = T,append = F)



