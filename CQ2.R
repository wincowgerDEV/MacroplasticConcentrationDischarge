library(dataRetrieval)
library(dplyr)
library(ggplot2)
library(mgcv)
library(readxl)
library(data.table)
library(stringr)
library(viridis)
library(tidyr)
library(MASS)
library(fitdistrplus)
library(EcoHydRology)
library(hydrostats)
library(plotly)

options(scipen = 999)


BootMean <- function(data) {
  B <- 10000
  mean <- numeric(B)
  n = length(data)
  
  set.seed(34345)
  for (i in 1:B) {
    boot <- sample(1:n, size=n, replace = TRUE)
    mean[i] <- mean(data[boot])
  }
  return(quantile(mean, c(0.025, 0.5, 0.975), na.rm = T))
}

BootSlope <- function(x, y) {
  B <- 10000
  mean <- numeric(B)
  n = length(x)
  
  set.seed(34345)
  for (i in 1:B) {
    boot <- sample(1:n, size=n, replace = TRUE)
    mean[i] <- lm(y[boot] ~ x[boot])$coefficients[2]
  }
  return(quantile(mean, c(0.025, 0.5, 0.975), na.rm = T))
}

BootIntercept <- function(x, y) {
  B <- 10000
  mean <- numeric(B)
  n = length(x)
  
  set.seed(34345)
  for (i in 1:B) {
    boot <- sample(1:n, size=n, replace = TRUE)
    mean[i] <- lm(y[boot] ~ x[boot])$coefficients[1]
  }
  return(quantile(mean, c(0.025, 0.5, 0.975), na.rm = T))
}

#Precip data https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelationships/Data/Raw Data")
precip <- read.csv("Data/72286903171Precip.csv", stringsAsFactors = F)

#GageFieldMeasurements <- read.csv("G:/My Drive/GrayLab/Projects/Plastics/Articles Publish/Active/Stream Modeling/Riverside/GageData.csv")
sampledata <- read.csv("Data/SampleMetaData.csv")
param_cd <- read.csv("Data/param_cd.csv")
#Particles that do not float.
NoFloat <- read.csv("Data/MasterNoFloat_Clean.csv")
Q <- read.csv("Data/SiteQ11066460.csv")
measurements <- read.csv("Data/measurements.csv")
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelationships/Data/Processed Data/Riverside")

RatingCurve <- read_excel("Data/Riverside/GageData.xlsx", sheet = "Rating Curve")
ggplot() + geom_line(data = RatingCurve, aes(x = INDEP, y = DEP)) + geom_point(data = measurements, aes(x = gage_height_va, y = chan_discharge, color = gage_va_change))

#Particles with sinking removed.
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelationships/Data/Raw Data/Santa Ana River SamplesNoFloatRemoved")
csvfiles_1 <- list.files(path = "Data/Santa Ana River SamplesNoFloatRemoved", pattern = ".csv", recursive = T, full.names = T)

filelistcsv_1 <- csvfiles_1 %>%
  lapply(read.csv)

csvfiles <- csvfiles_1
filelistcsv <- filelistcsv_1
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelationships/Data/Processed Data/Riverside")
PlasticMeasurements <- read_excel("Data/Riverside/SantaAna.xlsx", sheet = "PlasticMeasurements")
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/RouseProfilePaper/Data/Processed Data/DataForPublication")
#WSData <- read_excel("Data/Data_WaldschlaegerEdited.xlsx")
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelationships/Data/Processed Data/ParticleSizeConversionData/Plastic Masses")
Master <- read.csv("Data/ParticleSizeConversionData/Plastic Masses/MasterSheet - Sheet1.csv")
Mast <- Master[complete.cases(Master$Mass),]
massmodel <- gam(log10(Mass) ~ log10(Area), data = Mast)

ggplot(Mast, aes(y = log10(Mass), x =  log10(Area))) + geom_point() + geom_smooth(method = "lm")
#Do a test to make sure that 17 particles were removed.
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelam3 stionships/Data/Processed Data/Lab/Santa Ana River Samples")
csvfiles_2 <- list.files("Data/Lab/Macroplastic/Santa Ana River Samples",  pattern = ".csv", recursive = T, full.names = T)
filelistcsv_2 <- csvfiles_2 %>%
  lapply(fread)

#strange, this is looking like more particles are in the updated data. Why? what did I do back then lol.
#Took a look at 3:04 for settling particles, not sure where they are in the image that I have for them. 
#Went through and cleaned the data on 7/12 to make sure it was all in line. Should have all settling particles removed now. not sure what happened before. 
for(n in 1:length(csvfiles_1)){
  print(csvfiles_1[n])
  print(csvfiles_2[n])
  print(nrow(filelistcsv_1[[n]]) - nrow(filelistcsv_2[[n]]))
}

nrow(rbindlist(filelistcsv_1, fill = T)) - nrow(rbindlist(filelistcsv_2, fill = T))

#Don't run this code again, just useful for first grab.
#sites <- c("11066460")
#parameterCd <- "00065"
#startDate <- "1988-01-01"
#endDate <- "2020-01-01"

#availabledata <- whatNWISdata(siteNumber = sites)
#availabledata <- availabledata %>%
#  mutate(parm_cd = as.numeric(parm_cd)) %>%
#  left_join(param_cd)
#ggplot(availabledata, aes(x = begin_date)) + geom_histogram()

#Q <- readNWISuv(sites, parameterCd, startDate, endDate, tz = "America/Los_Angeles")

#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelationships/Data/Raw Data")
#write.csv(Q, "SiteQ11066460.csv")

#ratingdata <- readNWISrating(sites, "base")
#measurements <- readNWISmeas(sites, expanded = T, tz = "America/Los_Angeles")
#write.csv(measurements, "measurements.csv")

velocitycurve <- gam(log10(chan_velocity) ~ log10(chan_discharge), data = measurements)
widthcurve <- gam(log10(chan_width) ~ log10(chan_discharge), data = measurements)
areacurve <- gam(log10(chan_area) ~ log10(chan_discharge), data = measurements)

#Discharge figure ----
Discharge <- Q %>%
  mutate(dateTime = as.POSIXct(strptime(dateTime, format = "%Y-%m-%d %H:%M:%S"), tz = "America/Los_Angeles")) %>%
  rename(INDEP = X_00065_00000) %>%
  left_join(RatingCurve)

ggplot(Discharge)+ 
  geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1) + 
  scale_x_datetime(limits = c(as.POSIXct("2018-10-01 00:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-10-30 23:59:00", tz="America/Los_Angeles"))) + labs(y = "Discharge (cms)", x = "Date Time") + 
  geom_vline(xintercept = c(as.POSIXct("2019-02-02 12:00:00", tz="America/Los_Angeles"), as.POSIXct("2019-01-19 12:00:00", tz="America/Los_Angeles"), as.POSIXct("2019-08-09 12:00:00", tz="America/Los_Angeles"), as.POSIXct("2019-10-21 12:00:00", tz="America/Los_Angeles"), as.POSIXct("2019-2-13 12:00:00", tz="America/Los_Angeles"))) + 
  scale_y_log10(breaks = c(0.1,1,10,100,1000), limits = c(0.1, 1000)) + 
  theme_gray()

#precip is in mm, not sure what date time format is. https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf also probably want to look into the timestep, they are not all 1 hr. It seems like they all relate to the previous hour though so that may be helpful for extracting more granular data from this. May need to average everything within the hour period. 
#Not sure I still need this. not discussing lag time or anything like that. 
#precip2 <- precip %>%
#  dplyr::select(DATE, AA1) %>%
#  filter(grepl("01,", AA1)) %>%
#  mutate(DATE = as.POSIXct(strptime(DATE, format = "%Y-%m-%dT%H:%M:%S"), tz = "America/Los_Angeles")) %>%
#  filter(AA1 != "01,0000,9,5")
  
#for(row in 1:nrow(precip2)) {
#  precip2[row, "precip_mm"] <- strsplit(precip2[row, "AA1"], ",")[[1]][2]
#  precip2[row, "quality_code"] <- strsplit(precip2[row, "AA1"], ",")[[1]][4]
#  precip2[row, "condition_code"] <- strsplit(precip2[row, "AA1"], ",")[[1]][3]
#}

#precip2$precip_mm <- as.numeric(precip2$precip_mm)
#precip2 <- filter(precip2, precip_mm >  0)
#ggplot(precip2, aes(x = DATE, y = precip_mm)) + geom_point() + theme_classic()
#summary(precip2)

#Get particle size info for samples

#Compare size classes float and nofloat ----
datamerge <- data.frame(ParticleID = numeric(), Area = numeric(), SampleName = character(), stringsAsFactors = F)

n = 0
for(csv in csvfiles) {
  n = n + 1
  samplename <- str_match(csv, ".*/(.*)/")[2]
  if(!"Area" %in% colnames(filelistcsv[[n]])) next
  df <- data.frame(ParticleID = filelistcsv[[n]]$X, Area = filelistcsv[[n]]$Area, SampleName = samplename, stringsAsFactors = F)
  datamerge <- rbind(df, datamerge)
}

RemoveThese <- c("Santa Ana 2 2 5 35pm 1 min bomb 24.67",
                 "Santa Ana 2-13-19 6 30pm 54 min integrate 8.40g",
                 "Santa Ana Right Bank 2 2 1 45 10 min", 
                 "win cowger Sample Ballona 11 56am 3 20 19 30 min halfnet", 
                 "Santa Ana 2-13-19 6:30pm 54 min integrate", 
                 "Santa Ana 2/2 5:35pm 1 min bomb", 
                 "Santa Ana 1-17-19 3:04pm 10 min 2/2", 
                 "Santa Ana 1-17-19 3:04pm 10 min 2/3", 
                 "Santa Ana 1-17-19 3:04pm 10 min 2/4", 
                 "Santa Ana Right Bank 2/2 1:45 10 min")

MasterListCompare <- datamerge %>%
  dplyr::filter(!SampleName %in% RemoveThese) %>%
  dplyr::filter(Area > 1) %>%
  mutate(Volume = Area * sqrt(Area)) %>% #Volume is in mm3
  mutate(Mass = Volume * 0.03 * 0.001) #Come back to density to see if it should be 0.03 g/ml

#Look at Particle Sizes, Double check that this has no settling in the buoyant class.
sizecompare <- MasterListCompare %>%
  dplyr::select(Area) %>%
  mutate(class = "Float") %>%
  rbind(dplyr::select(NoFloat, Area, SampleName) %>%
          dplyr::filter(!SampleName %in% RemoveThese) %>%
          dplyr::filter(Area > 1) %>%
          mutate(class = "nofloat") %>%
          dplyr::select(Area, class))

ggplot(sizecompare, aes(x = sqrt(Area), y = class)) + geom_violin() + geom_boxplot(width=.2, notch = T) + scale_x_log10(limits = c(1, 1000), breaks = c(1, 10,100,1000)) + theme_gray()

sizecompare %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(count = n())

#Clean concentration data ----
MasterList <- datamerge %>%
  dplyr::filter(!SampleName %in% RemoveThese) %>%
  dplyr::filter(Area > 25) %>%
  mutate(Volume = Area * sqrt(Area)) %>% #Volume is in mm3
  mutate(Mass = Volume * 0.03 * 0.001) #Come back to density to see if it should be 0.03 g/ml

sampledataclean_pre <- sampledata %>%
  mutate_if(is.factor, as.character) %>%
  mutate(dateTime = as.POSIXct(strptime(paste(Date, Start.Time, sep = " "), format = "%m/%d/%Y %H:%M"), tz = "America/Los_Angeles")) %>%
  rename(SampleName = Sample.ID) %>%
  mutate(chan_discharge = approx(Discharge$dateTime, Discharge$DEP, xout = dateTime, rule = 2, method = "linear", ties=mean)[[2]]) %>%
  right_join(MasterList) 

sampledataclean_pre %>%
  dplyr::select(chan_discharge, SampleName, dateTime) %>%
  distinct() %>%
  ggplot() + geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1)+ geom_point(aes(y = chan_discharge * 0.0283168, x = dateTime), color = "red", size = 3, shape = 23, fill = "red") + scale_x_datetime(limits = c(as.POSIXct("2019-02-02 01:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-02-03 12:00:00", tz="America/Los_Angeles"))) + ylim(0,300) + labs(y = "Discharge (cms)", x = "Date Time") + theme_gray()

sampledataclean_pre %>%
  dplyr::select(chan_discharge, SampleName, dateTime) %>%
  distinct() %>%
  ggplot() + geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1)+ geom_point(aes(y = chan_discharge * 0.0283168, x = dateTime), color = "red", size = 3, shape = 23, fill = "red") + scale_x_datetime(limits = c(as.POSIXct("2019-01-17 01:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-01-18 12:00:00", tz="America/Los_Angeles"))) + ylim(0,300) + labs(y = "Discharge (cms)", x = "Date Time") + theme_gray()

sampledataclean_pre %>%
  dplyr::select(chan_discharge, SampleName, dateTime) %>%
  distinct() %>%
  ggplot() + geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1)+ geom_point(aes(y = chan_discharge * 0.0283168, x = dateTime), color = "red", size = 3, shape = 23, fill = "red") + scale_x_datetime(limits = c(as.POSIXct("2019-02-13 12:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-02-14 01:00:00", tz="America/Los_Angeles"))) + ylim(0,10) + labs(y = "Discharge (cms)", x = "Date Time") + theme_gray()

sampledataclean_pre$chan_velocity <- 10^predict.gam(velocitycurve, sampledataclean_pre, type = "response")
sampledataclean_pre$chan_area <- 10^predict.gam(areacurve, sampledataclean_pre, type = "response")
sampledataclean_pre$chan_width <- 10^predict.gam(widthcurve, sampledataclean_pre, type = "response")

RunoffSamples <- c("Santa Ana 2-2-19 3 32pm 10 mins 2 2 181.63g HH", "Santa Ana 2 2 3 51pm 10min", "Santa Ana 2-2-19 4 29pm 2 mins 277.2g HH", "Santa Ana 2 2 54 seconds 5 03 720.25 grams", "Santa Ana 1-17-19 7 26 3 min HH", "Santa Ana 1-17-19 6 30pm 1 min 27 sec ST HH", "Santa Ana 1-17-19 5 00pm 30 sec ST", "Santa Ana 1-17-19 4 29pm 2 min 2 2" )

sampledataclean <- sampledataclean_pre %>%
  mutate(chan_depth = chan_area/chan_width* 0.3048) %>% #in meters
  mutate(SampleSize = ifelse(chan_depth < 0.4, chan_depth * 0.4 * chan_velocity * 0.3048 * Duration..min.*60, 0.5*0.4* chan_velocity * 0.3048 * Duration..min.*60)) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  mutate(particlespersecond = 1 /(Duration..min. *60)) %>%
  mutate(concentration = 1/SampleSize) %>%
  mutate(chan_discharge_m = chan_discharge * 0.0283168) %>%
  mutate(shear_velocity = sqrt(9.8*0.003954717*chan_depth)) %>% #in meters
  mutate(proportion_sampled = ifelse(chan_depth < 0.4,1 , 0.4/chan_depth)) %>%
  mutate(Runoff = ifelse(SampleName %in% RunoffSamples, "Runoff", "Nonrunoff")) %>%
  mutate(particlemass = 10^predict.gam(massmodel, ., type = "response"))

#Just a qaqc proceedure for making sure we got rid of all the settling particles again. 
sampledatacleannofloat <- NoFloat %>%
  inner_join(sampledataclean, by = "SampleName") %>%
  mutate(AreaDiff = Area.x - Area.y) %>%
  dplyr::filter(Area.x > 25) %>%
  group_by(Area.x) %>%
  dplyr::filter(abs(AreaDiff) == min(abs(AreaDiff))) %>%
  ungroup()#Shows12 that all other particles are not the same size. 

#sampledataclean$particlemass <- 10^predict.gam(massmodel, sampledataclean, type = "response")

#Fluxbyparticlesize showing discharge variation

#Measured Rising Velocities ----
#in meters
velocities_macorplastic <- c(13.59, 8.60, 16.19, 15.87, 4.21, 5.56, 9.55, 2.99, 6.56, 2.21, 4.08, 4.93, 6.11) / 10
mean(velocities_macorplastic)
median(velocities_macorplastic)
min(velocities_macorplastic)
max(velocities_macorplastic)

RouseNumbers <- expand.grid(rising_vel_m_s = velocities_macorplastic, shear_vel_m_s = unique(sampledataclean$shear_velocity)) %>%
  mutate(RouseNum = rising_vel_m_s/shear_vel_m_s)

hist(RouseNumbers$RouseNum)
quantile(RouseNumbers$RouseNum, seq(0,1, by = 0.01))
ecdf(RouseNumbers$RouseNum)

#Particle size comparison runoff-baseflow ----
ggplot(sampledataclean, aes(sqrt(Area), color = Runoff)) + stat_ecdf(size = 3) + scale_color_viridis_d() + scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c(1, 10, 100, 1000), limits = c(1,1000)) + theme_gray() 

runofffit <- filter(sampledataclean, Runoff == "Runoff")
nonrunofffit <- filter(sampledataclean, Runoff == "Nonrunoff")

#null hypothesis is that they are drawn from the same distribution. 
ks.test(runofffit$Area, nonrunofffit$Area, alternative = "two.sided")
#This demonstrates that pdfs are the same

#Not really enough particles in many of these sample days to say much, what does 10 particles really tell us about the particle size distribution?
#ggplot(sampledataclean, aes(sqrt(Area), color = Runoff)) + geom_density(size = 3) + scale_color_viridis_d() + scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c(1, 10, 100, 1000), limits = c(1,1000)) + theme_gray() + facet_wrap(.~Date, nrow = 5)

sampledataclean %>%
  group_by(Runoff, SampleName) %>%
  summarise(count = n())

sampledataclean %>%
  group_by(SampleName) %>%
  summarise(count = n())

sampledataclean %>%
  group_by(Runoff) %>%
  summarise(count = n())

#Corrected concentrations, depth integrated. 
totalConcentrationDischarge <- sampledataclean %>%
  group_by(SampleName, chan_discharge_m, Mass..g., proportion_sampled, SampleSize, dateTime, Duration..min., Date, chan_depth, shear_velocity, Runoff) %>%
  dplyr::summarise(count = n(), sumarea = sum(Area), summass = sum(particlemass, na.rm = T), fiftyquantile = quantile(Area, 0.5), eightyfivequantile = quantile(Area, 0.85), fifteenquantile = quantile(Area, 0.15), minsize = min(Area), maxsize = max(Area)) %>%
  mutate(countconcentration = count/SampleSize*proportion_sampled) %>%
  mutate(areaconcentration = sumarea/SampleSize*proportion_sampled) %>%
  mutate(massconcentration = summass/SampleSize*proportion_sampled) %>%
  mutate(measured_mass = Mass..g./SampleSize*proportion_sampled) %>%
  arrange(desc(dateTime))

#ggplot(totalConcentrationDischarge, aes(x = areaconcentration, y = areaconcentrationadj)) + geom_abline() + geom_point() + scale_y_log10() + scale_x_log10()
ggplot(totalConcentrationDischarge, aes(x = summass, y = Mass..g.)) + geom_point() + scale_y_log10() + scale_x_log10() + geom_abline(intercept = 0, slope = 1)

#Linear Relatoinship, concentration - discharge

ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = areaconcentration, color = log(areaconcentration/countconcentration))) +
  geom_smooth(color = "black") + 
  geom_point() + 
  scale_color_viridis_c() +
  scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + 
  scale_y_log10(breaks = c(1,10,100,1000,10000), labels = c(1,10,100,1000,10000), limits = c(1,10000))  +
  theme_gray(base_size = 18) + 
  labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Area Concentration ("~mm^2~m^-3~")"))+ 
  coord_fixed()


ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = countconcentration, color = log(areaconcentration/countconcentration))) +
  geom_smooth(color = "black") + 
  geom_point() + 
  scale_color_viridis_c() +
  scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + 
  scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100), limits = c(0.01,100)) + 
  theme_gray(base_size = 18) + 
  labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ 
  coord_fixed()

ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = areaconcentration/countconcentration)) + geom_point() + geom_smooth(method = "lm", color = "black") #+ scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100), limits = c(0.01,100)) + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ coord_fixed()

#Adjust the concentrations for the area per particle to all be set to the mean. 
hist(log10(totalConcentrationDischarge$areaconcentration/totalConcentrationDischarge$countconcentration))
mean(totalConcentrationDischarge$areaconcentration/totalConcentrationDischarge$countconcentration)

correctedcount <- as.vector(totalConcentrationDischarge$countconcentration/(totalConcentrationDischarge$areaconcentration/totalConcentrationDischarge$countconcentration)/mean(totalConcentrationDischarge$areaconcentration/totalConcentrationDischarge$countconcentration))

ggplot() + geom_point(aes(x = totalConcentrationDischarge$chan_discharge_m, y = correctedcount)) + geom_smooth(aes(x = totalConcentrationDischarge$chan_discharge_m, y = correctedcount), color = "black") #+ scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100), limits = c(0.01,100)) + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ coord_fixed()

ggplot(totalConcentrationDischarge, aes(x = areaconcentration, y = countconcentration)) +
  geom_smooth(method = "lm", color = "black") + 
  geom_point() + 
  scale_color_viridis_c() +
  scale_x_log10() + 
  scale_y_log10()  +
  theme_gray(base_size = 18)

#Hysteresis behavior
ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = areaconcentration))  + geom_path(aes(color = Date), size = 2) + geom_point() + scale_color_viridis_d() + scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(1,10,100,1000,10000), labels = c(1,10,100,1000,10000), limits = c(1,10000)) + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Area Concentration ("~mm^2~m^-3~")")) + coord_fixed()# + geom_text(aes(x = chan_discharge, y = areaconcentration,label = SampleName))
ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = countconcentration)) + geom_path(aes(color = Date), size = 2) + geom_point() + scale_color_viridis_d() + scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100), limits = c(0.01,100)) + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ coord_fixed()#+ geom_text(aes(label = SampleName))
ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = measured_mass)) + geom_path(aes(color = Date), size = 2) + geom_point() + scale_color_viridis_d() + scale_x_log10() + scale_y_log10() + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Mass Concentration ("~num^1~m^-3~")"))+ coord_fixed()#+ geom_text(aes(label = SampleName))
ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = massconcentration)) + geom_path(aes(color = Date), size = 2) + geom_point() + scale_color_viridis_d() + scale_x_log10() + scale_y_log10() + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Mass Concentration ("~num^1~m^-3~")"))+ coord_fixed()#+ geom_text(aes(label = SampleName))

mass_concentration_model <- lm(log10(totalConcentrationDischarge$massconcentration) ~ log10(totalConcentrationDischarge$chan_discharge_m))
mass_concentration_model_slope_boot <- BootSlope(y = log10(totalConcentrationDischarge$massconcentration), x = log10(totalConcentrationDischarge$chan_discharge_m))
mass_concentration_model_intercept_boot <- BootIntercept(y = log10(totalConcentrationDischarge$massconcentration), x = log10(totalConcentrationDischarge$chan_discharge_m))

summary(mass_concentration_model)

mass_concentration_model$coefficients[2]

predict(mass_concentration_model, c(10, 100))

ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = massconcentration)) + scale_color_viridis_d() + geom_point() + scale_x_log10() + scale_y_log10() + theme_gray() + labs(x = "Discharge (cms)", y = "Mass Concentration (g/m^3)")

#Estimate Flux ----

#Constant Mean
Flux <- Discharge %>%
  filter(dateTime > as.POSIXct("2018-10-01 00:00:00", tz="America/Los_Angeles") & dateTime < as.POSIXct("2019-09-30 23:59:00", tz="America/Los_Angeles")) %>%
  mutate(cubic_m_s = DEP * 0.0283168) 

mean_concentration_bootstrap <- BootMean(totalConcentrationDischarge$massconcentration)

constant_mean_metric_tonnes <- sum(Flux$cubic_m_s * mean(totalConcentrationDischarge$massconcentration) * 15 * 60, na.rm = T)/10^6
min_constant_mean_metric_tonnes <- sum(Flux$cubic_m_s * BootMean(totalConcentrationDischarge$massconcentration)[1] * 15 * 60, na.rm = T)/10^6
max_constant_mean_metric_tonnes <- sum(Flux$cubic_m_s * BootMean(totalConcentrationDischarge$massconcentration)[3] * 15 * 60, na.rm = T)/10^6


#10^BootMean(log10(totalConcentrationDischarge$massconcentration))

Flux$concentration_estimate <- 10^((log10(Flux$cubic_m_s)*
                                       coef(mass_concentration_model)[2] + 
                                       coef(mass_concentration_model)[1])) * 
                                  10^(mean(mass_concentration_model$residuals^2)/2)

Flux$concentration_estimate_lwr <- 10^(log10(Flux$cubic_m_s)*mass_concentration_model_slope_boot[1] + 
                                          mass_concentration_model_intercept_boot[1]) * 
                                  10^(mean(mass_concentration_model$residuals^2)/2)

Flux$concentration_estimate_upr <- 10^(log10(Flux$cubic_m_s)*mass_concentration_model_slope_boot[3] + 
                                         mass_concentration_model_intercept_boot[3]) *
                                        10^(mean(mass_concentration_model$residuals^2)/2)

Flux_uncertainty <- Flux %>%
  mutate(uncertain_concentration = ifelse(concentration_estimate > max(totalConcentrationDischarge$massconcentration) | concentration_estimate < min(totalConcentrationDischarge$massconcentration), 0.1, 0)) %>% #Results in no correction because concentrations are within observed domain
  mutate(uncertain_discharge = ifelse(cubic_m_s > max(measurements$discharge_va, na.rm = T) * 0.0283168 | cubic_m_s < min(measurements$discharge_va, na.rm = T) * 0.0283168, 0.1, 0)) %>%
  mutate(lwr_uncertainty_ci_others = concentration_estimate_lwr - uncertain_discharge*concentration_estimate) %>%
  mutate(upr_uncertainty_ci_others = concentration_estimate_upr + uncertain_discharge*concentration_estimate)

#Flux estimate using regression
discharge_regression_metric_tonnes <- sum(Flux_uncertainty$concentration_estimate * Flux$cubic_m_s * 15 * 60, na.rm = T)/10^6
discharge_regression_metric_tonnes_min <-  sum(Flux_uncertainty$lwr_uncertainty_ci_others * Flux$cubic_m_s * 15 * 60, na.rm = T)/10^6
discharge_regression_metric_tonnes_max <-  sum(Flux_uncertainty$upr_uncertainty_ci_others * Flux$cubic_m_s * 15 * 60, na.rm = T)/10^6

runoff_event_ranges = tibble(
  start = c("2018-10-04 02:15:00", 
            "2018-10-12 22:30:00", 
            "2018-11-29 10:30:00", 
            "2018-12-06 13:00:00", 
            "2019-01-14 13:00:00",
            "2019-01-31 15:15:00",
            "2019-02-02 15:15:00",
            "2019-02-14 01:45:00",
            "2019-02-17 06:45:00",
            "2019-02-21 01:00:00",
            "2019-03-02 05:00:00",
            "2019-03-06 08:00:00",
            "2019-03-20 14:45:00",
            "2019-05-16 09:30:00",
            "2019-05-19 03:45:00",
            "2019-05-22 12:00:00",
            "2019-05-27 00:45:00"),
  end = c( "2018-10-04 12:00:00",
            "2018-10-14 14:30:00",
            "2018-11-30 22:30:00",
            "2018-12-08 07:30:00",
            "2019-01-19 15:15:00",
            "2019-02-01 13:30:00",
            "2019-02-06 09:00:00",
            "2019-02-15 17:45:00",
            "2019-02-18 20:30:00",
            "2019-02-22 08:00:00",
            "2019-03-04 01:30:00",
            "2019-03-08 04:00:00",
            "2019-03-21 11:45:00",
            "2019-05-16 20:45:00",
            "2019-05-20 19:45:00",
            "2019-05-24 13:00:00",
            "2019-05-27 13:15:00")
  )

runoff_dates <- Flux[0,]

for(n in 1:nrow(runoff_event_ranges)){
  runoff_dates <- Flux %>%
    dplyr::filter(dateTime > as.POSIXct(unlist(runoff_event_ranges[n,1]), tz="America/Los_Angeles") & dateTime < as.POSIXct(unlist(runoff_event_ranges[n,2]), tz="America/Los_Angeles")) %>%
    bind_rows(runoff_dates)
}

baseflow_dates <- Flux %>%
  anti_join(runoff_dates)

#Flux during runoff vs baseflow
BootMean(totalConcentrationDischarge$massconcentration[totalConcentrationDischarge$Runoff == "Runoff"])
BootMean(totalConcentrationDischarge$massconcentration[totalConcentrationDischarge$Runoff == "Nonrunoff"])

runoff_flux <- sum(runoff_dates$cubic_m_s * (mean(totalConcentrationDischarge$massconcentration[totalConcentrationDischarge$Runoff == "Runoff"])) * 15 * 60)/10^6
baseflow_flux <- sum(baseflow_dates$cubic_m_s * (mean(totalConcentrationDischarge$massconcentration[totalConcentrationDischarge$Runoff == "Nonrunoff"])) * 15 * 60, na.rm = T)/10^6

runoff_baseflow_flux = runoff_flux + baseflow_flux

runoff_flux_min <- sum(runoff_dates$cubic_m_s * (BootMean(totalConcentrationDischarge$massconcentration[totalConcentrationDischarge$Runoff == "Runoff"])[1]) * 15 * 60)/10^6
baseflow_flux_min <- sum(baseflow_dates$cubic_m_s * (BootMean(totalConcentrationDischarge$massconcentration[totalConcentrationDischarge$Runoff == "Nonrunoff"])[1]) * 15 * 60, na.rm = T)/10^6

runoff_baseflow_flux_min = runoff_flux_min + baseflow_flux_min


runoff_flux_max <- sum(runoff_dates$cubic_m_s * (BootMean(totalConcentrationDischarge$massconcentration[totalConcentrationDischarge$Runoff == "Runoff"])[3]) * 15 * 60)/10^6
baseflow_flux_max <- sum(baseflow_dates$cubic_m_s * (BootMean(totalConcentrationDischarge$massconcentration[totalConcentrationDischarge$Runoff == "Nonrunoff"])[3]) * 15 * 60, na.rm = T)/10^6

runoff_baseflow_flux_max = runoff_flux_max+ baseflow_flux_max

#Runoff flow
sum(runoff_dates$cubic_m_s * 15 * 60) 
sum(baseflow_dates$cubic_m_s * 15 * 60, na.rm = T)

figure_table <- tibble(
  annual_flux_tonnes = c(runoff_baseflow_flux, constant_mean_metric_tonnes, discharge_regression_metric_tonnes),
  min_flux_tonnes    = c(runoff_baseflow_flux_min, min_constant_mean_metric_tonnes, discharge_regression_metric_tonnes_min),
  max_flux_tonnes    = c(runoff_baseflow_flux_max, max_constant_mean_metric_tonnes, discharge_regression_metric_tonnes_max),
  name               = c("runoff baseflow separation", "constant mean", "regression")
)

ggplot(figure_table, aes(y = name, x = annual_flux_tonnes)) + 
  geom_point() + 
  geom_errorbar(aes(xmin = min_flux_tonnes, xmax = max_flux_tonnes)) + 
  scale_x_log10(limits = c(0.01, 1000)) + 
  theme_gray() + 
  labs(y = "", x = "Annual Flux (metric tonnes)")

#approximately 20 times more flux during runoff periods than dry periods. Suggests issues with current managment strategy.

#Only 7 % of the time are we in these runoff periods, but they account for 10X of the flux. 
nrow(runoff_dates)/nrow(Flux)
