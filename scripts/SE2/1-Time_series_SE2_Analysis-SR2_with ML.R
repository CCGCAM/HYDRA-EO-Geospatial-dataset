rm(list=ls())

#::::::::::::::::::::::::::::::::::::::::::::::::::::
#::	0. load main Libraries   -----
#::::::::::::::::::::::::::::::::::::::::::::::::::::

if (!require("ggplot2")) { install.packages("ggplot2"); require("ggplot2") }  ###
if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ###
if (!require("lubridate")) { install.packages("lubridate"); require("lubridate") }  ###
if (!require("tidyverse")) { install.packages("tidyverse"); require("tidyverse") }  ###
if (!require("zoo")) { install.packages("zoo"); require("zoo") }  ###
if (!require("cowplot")) { install.packages("cowplot"); require("cowplot") }  ###

required_packages <- c("dplyr", "tidyr", "ggplot2","parallel", "doParallel",'reshape2','ggpubr','ggExtra')

missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)

require(ToolsRTM)
#::::::::::::::::::::::::::::::::::::::::::::::::::::
#	1. Import data   -----
#::::::::::::::::::::::::::::::::::::::::::::::::::::

# Data frames
df <- data.table::fread("Tables/results/TimeSerie_SE2_inSR2.csv", header = T, sep = ",")

#::::::::::::::::::::::::::::::::::::::::::::::::::::
#::	2. Getting the mean values in the Time Series
#::::::::::::::::::::::::::::::::::::::::::::::::::::


df.mean <- df %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(
    Year = as.numeric(year(Date)),                         # Extract year from Date
    Month = month(Date),                       # Extract month from Date
    Period = ifelse(Year >= 2020, 2,1),         # Adding a column for seprating both periods
    Area = ifelse(ID == 1, "Infected", ifelse(ID == 2, "Healthy", NA)) # Assign Area
  ) %>%
  group_by(Date, Area) %>%
  summarise_all(mean) %>%
  # Calculate the mean for all variables
  ungroup() %>%                                # Ungroup after summarization
  dplyr::select(ID,Area,Date,Year, Month, Period,everything())    # Reorder columns


dim(df.mean)
head(df.mean)
names(df.mean)


df.mean["SBI"] <- 0.3037 * df.mean["B2"] + 0.2793 * df.mean["B3"] +
  0.4743 * df.mean["B4"] + 0.5585 * df.mean["B8"] + 0.5082 * df.mean["B11"] +
  0.1863 * df.mean["B2"]
df.mean["GVI"] <- -0.2848 * df.mean["B2"] - 0.2435 * df.mean["B3"] -
  0.5436 * df.mean["B4"] + 0.7243 * df.mean["B8"] + 0.084 * df.mean["B11"] -
  0.18 * df.mean["B12"]
df.mean["WET"] <- 0.1509 * df.mean["B2"] + 0.1973 * df.mean["B3"] +
  0.3279 * df.mean["B4"] + 0.3406 * df.mean["B8A"] - 0.7112 * df.mean["B11"] -
  0.4572 * df.mean["B12"]
df.mean["BF.Anth"] <- (df.mean["B3"] - df.mean["B2"])/(df.mean["B3"] +
                                             df.mean["B2"])
df.mean["CR.red.nir"] <- df.mean["B7"]/(df.mean["B8"] + ((df.mean["B6"] -
                                                            df.mean["B8"])/(740.5 - 833)) * (782.8 - 832.8))

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


inputs <- c(SE2.bands,colnames(indices))


LUT.ml <- cbind(df.LUT,df.convolved,indices)



target_trait ='Cab'

# Build the modelling dataset: 1 column for Cab + predictors
dataset <- cbind(
  Cab = LUT.ml[[target_trait]],
  LUT.ml[, inputs]
)

# Remove any duplicated columns (can happen if names overlap)
dataset <- dataset[, !duplicated(names(dataset))]
# Remove columns that are entirely NA
dataset <- dataset[, colSums(is.na(dataset)) == 0]

head(dataset)

# ----------------------------------------------------------------------------
# 8. Optional: remove highly collinear predictors using VIF
# ----------------------------------------------------------------------------

set.seed(42)
rows.vif <- sample(nrow(dataset), min(200, nrow(dataset)))  # subset for speed

vif_keep <- ToolsRTM::getVIF(
  dataset[rows.vif, -1],   # all predictors (exclude first column = Cab)
  thresh = 10              # VIF threshold
)
print(vif_keep)


# Keep only predictors selected by VIF
dataset_vif <- cbind(Cab = dataset$Cab, dataset[, vif_keep])
summary(dataset_vif)

# ----------------------------------------------------------------------------
# 9. Inversion using GB &RF
# ----------------------------------------------------------------------------

model.rf=get.inversion(data = dataset_vif, depVar = "Cab", inputs = vif_keep,
                       algorithm  = "RF", n.samples = 500,seed = 123)

rf_fit <- model.rf$model   # caret::train object

plot(model.rf$importance)

# Predict Cab for each Plot/Date
roll_cols <-  df.mean %>%
  select(B2:last_col()) %>%
  names()
df.smooth <- df.mean %>%
  mutate(Date = as.Date(Date)) %>%
  arrange(ID, Area, Date) %>%
  group_by(ID) %>%
  mutate(across(all_of(roll_cols),
                ~ zoo::rollmean(.x, 3, align = "right", fill = NA))) %>%
  ungroup() %>%
  select(ID, Area, Date, Year, Month, Period, all_of(roll_cols)) %>%
  tidyr::drop_na(all_of(roll_cols))




Cab_pred.rf <- predict(rf_fit, df.smooth[, vif_keep])




df.plot <- cbind(df.smooth,Cab_pred.rf)

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


plot.cab <-ggplot(df.plot,aes(x=Date,y=Cab_pred.rf, color=Area)) +
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

scatter.cab <-ggplot(df.scatter,aes(x=CR.red.nir,y=Cab_pred.rf, color=as.factor(Period))) +
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


