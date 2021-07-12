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

header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-c(1:2),]
}

nathantheme <- function() {
  theme(legend.title=element_text(size=14, face="bold"))+
    theme(legend.text = element_text(size = 12, face="bold", family = "Times New Roman"))+
    theme(axis.line = element_line(size=1, linetype = "solid"))+
    theme(axis.text = element_text(size=12, face="bold", color = "black", family = "Times New Roman"))+
    theme(axis.ticks.length = unit(5, "pt"))+
    theme(axis.ticks= element_line(size=1, color = "black"))+
    theme(axis.title = element_text(size=14, face="bold", family = "Times New Roman"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size = 18, face = "bold", family = "Times New Roman"))+
    theme(axis.title.x = element_text(margin = margin(10,10,10,10,"pt"), family = "Times New Roman"))+
    theme(axis.title.y = element_text(margin = margin(10,10,10,10,"pt"), family = "Times New Roman"))+
    theme(panel.background =element_rect(fill="white"))+
    theme(panel.grid=element_line(color="grey100", size = 0.5))#+
    #annotation_logticks(base=10, size = 1)
}

bootks <- function(x, y){
  B <- 10000
  ksd <- numeric(B)
  nx = length(x)
  ny = length(y)
  
  set.seed(34345)
  for (i in 1:B) {
    bootx <- sample(x, size=nx, replace = TRUE)
    booty <- sample(y, size=ny, replace = TRUE)
    ksd[i] <- unname(ks.test(bootx, booty)$statistic)
  }
  return(ksd)
}

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

RouseAverageSurfaceSample <- function(Concentration, RouseNumber, SampleDepth, StreamDepth) {
  SampleLocation <- (StreamDepth-SampleDepth)/StreamDepth
  ReferenceLocation = (SampleLocation+1)/2
  Sample <- seq(SampleLocation, 0.99, by = 0.001)
  SampleRegion <- rep(Concentration, length(Sample))
  NormDepth <- seq(0.05,SampleLocation, by= 0.001) #Will break if the sample location is less than 0.05
  Profile <- Concentration * (((1-NormDepth)/NormDepth)*((ReferenceLocation)/(1-ReferenceLocation)))^RouseNumber
  ProfileFull <- c(Profile, SampleRegion)
  return(mean(ProfileFull))
}


logMaintain <- function(x) {
  x <- ifelse(x == 0, 0, ifelse(x > 0, log10(x*1000000000), -log10(abs(x*1000000000))))
  return(x)
}

reverselogMaintain <- function (x) {
  x <- ifelse(x == 0, 0, ifelse(x > 0, 10^x/1000000000, -1*10^abs(x)/1000000000))
  return(x)
}

RisingConcentrationEstimateMinMax <- function(Discharge, RisingConcentration, StreamDepth, SampleDepth, ShearVelocity, KnownRisingVelocity){
  ConcentrationLoopRising <- c()
  for(Vel in reverselogMaintain(KnownRisingVelocity)) {
    RouseNumberModel <- Vel/(0.4*ShearVelocity)
    if(RouseNumberModel > 2.5) { 
      next
    }
    else if(RouseNumberModel < -2.5) {
      ConcentrationEstimate = Discharge*(SampleDepth/StreamDepth)*RisingConcentration/Discharge 
    }
    else {
      ConcentrationEstimate = RouseAverageSurfaceSample(RisingConcentration, RouseNumberModel, SampleDepth, StreamDepth)
    }
    ConcentrationLoopRising <- c(ConcentrationLoopRising, ConcentrationEstimate)
  }
  RouseAverageConcentrationRising <- BootMean(ConcentrationLoopRising)
  return(c(RouseAverageConcentrationRising, max(ConcentrationLoopRising), min(ConcentrationLoopRising)))
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
WSData <- read_excel("Data/Data_WaldschlaegerEdited.xlsx")
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelationships/Data/Processed Data/ParticleSizeConversionData/Plastic Masses")
Master <- read.csv("Data/ParticleSizeConversionData/Plastic Masses/MasterSheet - Sheet1.csv") #Why isn't this being used as the particle size data? 
Mast <- Master[complete.cases(Master$Mass),]
massmodel <- gam(log10(Mass) ~ log10(Area), data = Mast)

#Do a test to make sure that 17 particles were removed.
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelationships/Data/Processed Data/Lab/Santa Ana River Samples")
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

Discharge <- Q %>%
  mutate(dateTime = as.POSIXct(strptime(dateTime, format = "%Y-%m-%d %H:%M:%S"), tz = "America/Los_Angeles")) %>%
  rename(INDEP = X_00065_00000) %>%
  left_join(RatingCurve)

options(scipen = 999)
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

RemoveThese <- c("Santa Ana 2 2 5 35pm 1 min bomb 24.67","Santa Ana 2-13-19 6 30pm 54 min integrate 8.40g","Santa Ana Right Bank 2 2 1 45 10 min", "Santa Ana 2-13-19 5:37 30 min Full Depth Sample 5mm net 1090.48g", "win cowger Sample Ballona 11 56am 3 20 19 30 min halfnet", "Santa Ana 2-13-19 6:30pm 54 min integrate", "Santa Ana 2/2 5:35pm 1 min bomb", "Santa Ana 1-17-19 3:04pm 10 min 2/2", "Santa Ana 1-17-19 3:04pm 10 min 2/3", "Santa Ana 1-17-19 3:04pm 10 min 2/4", "Santa Ana Right Bank 2/2 1:45 10 min")

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
  mutate(SampleSize = ifelse(chan_depth < 0.4, chan_depth * 0.4 * chan_velocity * 0.3048 * Duration..min.*60, 0.3*0.4* chan_velocity * 0.3048 * Duration..min.*60)) %>% #in cubic meters
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
  arrange(desc(dateTime))

#ggplot(totalConcentrationDischarge, aes(x = areaconcentration, y = areaconcentrationadj)) + geom_abline() + geom_point() + scale_y_log10() + scale_x_log10()
ggplot(totalConcentrationDischarge, aes(x = summass, y = Mass..g.)) + geom_point() + scale_y_log10() + scale_x_log10() + geom_abline(intercept = 0, slope = 1)

MassConcentration <- sampledataclean %>%
  dplyr::select(SampleName, Mass..g., SampleSize, chan_discharge_m, dateTime, Date) %>%
  distinct() %>%
  mutate(massconc = Mass..g./SampleSize) %>%
  arrange(desc(dateTime))

#Linear Relatoinship, concentration - discharge
ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = areaconcentration)) + geom_point() + geom_smooth(method = "lm", color = "black") + scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(1,10,100,1000,10000), labels = c(1,10,100,1000,10000), limits = c(1,10000))  + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Area Concentration ("~mm^2~m^-3~")")) + coord_fixed()
ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = countconcentration)) + geom_point() + geom_smooth(method = "lm", color = "black") + scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100), limits = c(0.01,100)) + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ coord_fixed()

#Power relationship
ggplot(totalConcentrationDischarge,aes(x = chan_discharge_m,y = areaconcentration)) +
  geom_point() + 
  stat_smooth(method = 'nls', formula = 'y~a*x^b', start = list(a = 1,b=1),se=FALSE) + scale_x_log10() + scale_y_log10()

shapiro.test(log10(totalConcentrationDischarge$countconcentration))

countmodel <- gam(log10(totalConcentrationDischarge$countconcentrationadj) ~ log10(totalConcentrationDischarge$chan_discharge_m))

summary(countmodel)

shapiro.test(log10(totalConcentrationDischarge$areaconcentrationadj))

areamodel <- gam(log10(totalConcentrationDischarge$areaconcentrationadj) ~ log10(totalConcentrationDischarge$chan_discharge_m))

summary(areamodel)

#Hysteresis behavior
ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = areaconcentration))  + geom_path(aes(color = Date), size = 2) + geom_point() + scale_color_viridis_d() + scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(1,10,100,1000,10000), labels = c(1,10,100,1000,10000), limits = c(1,10000)) + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Area Concentration ("~mm^2~m^-3~")")) + coord_fixed()# + geom_text(aes(x = chan_discharge, y = areaconcentration,label = SampleName))
ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = countconcentration)) + geom_path(aes(color = Date), size = 2) + geom_point() + scale_color_viridis_d() + scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100), limits = c(0.01,100)) + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ coord_fixed()#+ geom_text(aes(label = SampleName))
ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = massconcentration)) + geom_path(aes(color = Date), size = 2) + geom_point() + scale_color_viridis_d() + scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100), limits = c(0.01,100)) + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ coord_fixed()#+ geom_text(aes(label = SampleName))

ggplot(MassConcentration, aes(x = chan_discharge_m, y = massconc)) + scale_color_viridis_d() + geom_path(aes(color = Date), size = 2) + geom_point() + scale_x_log10() + scale_y_log10() + theme_gray() + labs(x = "Discharge (cms)", y = "Mass Concentration (g/m^3)")#+ geom_text(aes(label = SampleName))

#summofmassperareabin
sampledataclean %>%
  dplyr::group_by(equalintervalarea) %>%
  dplyr::summarise(Mass = sum(particlemass)) %>%
  ggplot(aes(x = equalintervalarea, y = Mass)) + scale_fill_viridis_d() + geom_bar(stat = "identity", position = "dodge") + scale_y_log10() + theme_classic() #+ geom_text(aes(label = count), position = position_dodge(width = 1))


