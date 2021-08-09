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


BootLM <- function(x,y) {
  B <- 10000
  models <- data.frame(intercept = numeric(B), slope = numeric(B))
  n = length(x)
  
  set.seed(34345)
  for (i in 1:B) {
    boot <- sample(1:n, size=n, replace = TRUE)
    model <- lm(y[boot] ~ x[boot])
    models[i, "intercept"] <- model$coefficients[1]
    models[i, "slope"] <- model$coefficients[2]
  }
  return(models)
}

PredictQuantileBootLM <- function(x,bootdf,minormax) {
  values = numeric(length = length(x))
  for(i in 1:length(x)){
      if(minormax == "min"){
        values[i] <- quantile(x[i]*bootdf$slope + bootdf$intercept, c(0.025))
      }
      else{
        values[i] <-  quantile(x[i]*bootdf$slope + bootdf$intercept, c(0.975))
      }
  }
  return(values)
}

RouseAverageSurfaceSample <- function(Concentration, RouseNumber, SampleDepth, StreamDepth) {
  SampleLocation <- (StreamDepth-SampleDepth)/StreamDepth
  ReferenceLocation = (SampleLocation+1)/2
  Sample <- seq(SampleLocation, 1, by = 0.00001)
  SampleRegion <- rep(Concentration, length(Sample))
  NormDepth <- seq(0.05,SampleLocation, by= 0.00001)
  Profile <- Concentration * (((1-NormDepth)/NormDepth)*((ReferenceLocation)/(1-ReferenceLocation)))^RouseNumber
  ProfileFull <- c(Profile, SampleRegion)
  return(mean(ProfileFull))
}

find_min_rouse <- function(shear, settling) {
  values <- numeric(length(shear))
  for(n in 1:length(shear)){
     values[n] <-  mean(settling/shear[n])
  }
  values
}

#PredictQuantileBootLM_vector <- Vectorize(PredictQuantileBootLM)

#Data ----
#Sample data
sampledata <- read.csv("Data/SampleMetaData.csv")
#Particles that do not float.
NoFloat <- read.csv("Data/MasterNoFloat_Clean.csv")
#Continuous discharge
Q <- read.csv("Data/SiteQ11066460.csv")
#Stage discharge relationships.
measurements <- read.csv("Data/measurements.csv") %>%
  mutate(Date = as.Date(measurement_dt)) %>%
  filter(Date > as.Date("2018-01-01")) %>%#this is better because more current.
  mutate(chan_depth_m = chan_area/chan_width * 0.3048) %>%
  mutate(chan_velocity_m = chan_velocity * 0.3048) %>%
  mutate(chan_discharge_m = chan_discharge * 0.0283168)
#need to correct the gage_height_va
RatingCurve <- read_excel("Data/Riverside/GageData.xlsx", sheet = "Rating Curve")

#Higher uncertainties when using rating curves for higher flows because few data points. 

ggplot(measurements, aes(x = gage_height_va, y = chan_depth_m)) + geom_point() + geom_smooth(method = "lm") + scale_y_log10() + scale_x_log10()
ggplot(measurements, aes(x = gage_height_va, y = chan_velocity_m)) + geom_point() + geom_smooth(method = "lm") + scale_y_log10() + scale_x_log10()
ggplot(measurements, aes(x = gage_height_va, y = chan_discharge)) + geom_point() + geom_line(data = RatingCurve, aes(x = INDEP+SHIFT, y = DEP) )+ geom_smooth(method = "lm") + scale_y_log10() + scale_x_log10()
#Particles with sinking removed.
csvfiles_1 <- list.files(path = "Data/Santa Ana River SamplesNoFloatRemoved", pattern = ".csv", recursive = T, full.names = T)

filelistcsv_1 <- csvfiles_1 %>%
  lapply(read.csv)

csvfiles <- csvfiles_1
filelistcsv <- filelistcsv_1

#Plastic Masses
Master <- read.csv("Data/ParticleSizeConversionData/Plastic Masses/MasterSheet - Sheet1.csv")

Mast <- Master[complete.cases(Master$Mass),] %>%
  filter(Area > 25)

massmodel <- lm(log10(Mass) ~ log10(Area), data = Mast)
massmodelcurverange <- BootLM(y = log10(Mast$Mass), x = log10(Mast$Area))

#massmodelboot <- BootLM(x = log10(Mast$Area), y = log10(Mast$Mass))
ggplot(Mast, aes(y = log10(Mass), x =  log10(Area))) + geom_point() + geom_smooth(method = "lm")

#Flow models ----
dischargecurve <- lm(log10(chan_discharge) ~ log10(gage_height_va), data = measurements)
dischargecurverange <- BootLM(y = log10(measurements$chan_discharge), x = log10(measurements$gage_height_va))

velocitycurve <- lm(log10(chan_velocity) ~ log10(gage_height_va), data = measurements)
velocitycurverange <- BootLM(y = log10(measurements$chan_velocity), x = log10(measurements$gage_height_va))

depthcurve <- lm(log10(chan_area/chan_width) ~ log10(gage_height_va), data = measurements)
depthcurverange <- BootLM(y = log10(measurements$chan_area/measurements$chan_width), x = log10(measurements$gage_height_va))

#areacurve <- gam(log10() ~ log10(gage_height_va), data = measurements)

#Discharge figure ----
Discharge <- Q %>%
  mutate(dateTime = as.POSIXct(strptime(dateTime, format = "%Y-%m-%d %H:%M:%S"), tz = "America/Los_Angeles")) %>%
  mutate(DEP = 10^(log10(X_00065_00000)*dischargecurve$coefficients[2] + dischargecurve$coefficients[1]) * 10^(mean(dischargecurve$residuals^2)/2))

ggplot(Discharge)+ 
  geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1) + 
  scale_x_datetime(limits = c(as.POSIXct("2018-10-01 00:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-10-30 23:59:00", tz="America/Los_Angeles"))) + labs(y = "Discharge (cms)", x = "Date Time") + 
  geom_vline(xintercept = c(as.POSIXct("2019-02-02 12:00:00", tz="America/Los_Angeles"), as.POSIXct("2019-01-19 12:00:00", tz="America/Los_Angeles"), as.POSIXct("2019-08-09 12:00:00", tz="America/Los_Angeles"), as.POSIXct("2019-10-21 12:00:00", tz="America/Los_Angeles"), as.POSIXct("2019-2-13 12:00:00", tz="America/Los_Angeles"))) + 
  scale_y_log10(breaks = c(0.1,1,10,100,1000), limits = c(0.1, 1000)) + 
  theme_gray()


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
                 "Santa Ana Right Bank 2/2 1:45 10 min") #Samples with notes not to use.

MasterListCompare <- datamerge %>%
  dplyr::filter(!SampleName %in% RemoveThese) %>%
  dplyr::filter(Area > 1) 

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
#Double checked this, these are only for samples which are included in this study for concentration discharge relationships. 
sizecompare %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(count = n())

#Clean concentration data ----
MasterList <- datamerge %>%
  dplyr::filter(!SampleName %in% RemoveThese) %>%
  dplyr::filter(Area > 25) 


sampledataclean_pre <- sampledata %>%
  mutate_if(is.factor, as.character) %>%
  mutate(dateTime = as.POSIXct(strptime(paste(Date, Start.Time, sep = " "), format = "%m/%d/%Y %H:%M"), tz = "America/Los_Angeles")) %>%
  rename(SampleName = Sample.ID) %>%
  mutate(gage_height = approx(Discharge$dateTime, Discharge$X_00065_00000, xout = dateTime, rule = 2, method = "linear", ties=mean)[[2]]) %>%
  mutate(chan_discharge = 10^(log10(gage_height)*dischargecurve$coefficients[2] + dischargecurve$coefficients[1]) * 10^(mean(dischargecurve$residuals^2)/2)) %>%
  mutate(chan_discharge_min = 10^(PredictQuantileBootLM(x = log10(gage_height), bootdf = dischargecurverange, minormax = "min")) * 10^(mean(dischargecurve$residuals^2)/2)) %>%
  mutate(chan_discharge_max = 10^(PredictQuantileBootLM(x = log10(gage_height), bootdf = dischargecurverange, minormax = "max")) * 10^(mean(dischargecurve$residuals^2)/2)) %>%
  mutate(chan_discharge_m = chan_discharge * 0.0283168) %>%
  mutate(chan_discharge_m_min = chan_discharge_min * 0.0283168) %>%
  mutate(chan_discharge_m_max = chan_discharge_max * 0.0283168) %>%
  right_join(MasterList) 

#PredictQuantileBootLM(x = log10(sampledataclean_pre$gage_height), bootdf = dischargecurverange, minormax = "min")

sampledataclean_pre %>%
  dplyr::select(chan_discharge_m, SampleName, dateTime) %>%
  distinct() %>%
  ggplot() + geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1)+ geom_point(aes(y = chan_discharge_m, x = dateTime), color = "red", size = 3, shape = 23, fill = "red") + scale_x_datetime(limits = c(as.POSIXct("2019-02-02 01:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-02-03 12:00:00", tz="America/Los_Angeles"))) + ylim(0,300) + labs(y = "Discharge (cms)", x = "Date Time") + theme_gray()

sampledataclean_pre %>%
  dplyr::select(chan_discharge_m, SampleName, dateTime) %>%
  distinct() %>%
  ggplot() + geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1)+ geom_point(aes(y = chan_discharge_m, x = dateTime), color = "red", size = 3, shape = 23, fill = "red") + scale_x_datetime(limits = c(as.POSIXct("2019-01-17 01:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-01-18 12:00:00", tz="America/Los_Angeles"))) + ylim(0,1000) + labs(y = "Discharge (cms)", x = "Date Time") + theme_gray()

sampledataclean_pre %>%
  dplyr::select(chan_discharge_m, SampleName, dateTime) %>%
  distinct() %>%
  ggplot() + geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1)+ geom_point(aes(y = chan_discharge_m, x = dateTime), color = "red", size = 3, shape = 23, fill = "red") + scale_x_datetime(limits = c(as.POSIXct("2019-02-13 12:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-02-14 01:00:00", tz="America/Los_Angeles"))) + ylim(0,10) + labs(y = "Discharge (cms)", x = "Date Time") + theme_gray()

#But don't these require the log10 of chan_discharge? Not sure these are working right. 
sampledataclean_pre$chan_velocity <- 10^(log10(sampledataclean_pre$gage_height)*velocitycurve$coefficients[2] + velocitycurve$coefficients[1]) * 10^(mean(velocitycurve$residuals^2)/2)
sampledataclean_pre$chan_velocity_min <- 10^(PredictQuantileBootLM(x = log10(sampledataclean_pre$gage_height), bootdf = velocitycurverange, minormax = "min")) * 10^(mean(velocitycurve$residuals^2)/2)
sampledataclean_pre$chan_velocity_max <- 10^(PredictQuantileBootLM(x = log10(sampledataclean_pre$gage_height), bootdf = velocitycurverange, minormax = "max")) * 10^(mean(velocitycurve$residuals^2)/2)

sampledataclean_pre$chan_depth <- 10^(log10(sampledataclean_pre$gage_height)*depthcurve$coefficients[2] + depthcurve$coefficients[1]) * 10^(mean(depthcurve$residuals^2)/2)
sampledataclean_pre$chan_depth_min <- 10^(PredictQuantileBootLM(x = log10(sampledataclean_pre$gage_height), bootdf = depthcurverange, minormax = "min")) * 10^(mean(depthcurve$residuals^2)/2)
sampledataclean_pre$chan_depth_max <- 10^(PredictQuantileBootLM(x = log10(sampledataclean_pre$gage_height), bootdf = depthcurverange, minormax = "max")) * 10^(mean(depthcurve$residuals^2)/2)

sampledataclean_pre$particlemass <- 10^(log10(sampledataclean_pre$Area)*massmodel$coefficients[2] + massmodel$coefficients[1]) * 10^(mean(massmodel$residuals^2)/2)
sampledataclean_pre$particlemass_min <- 10^(PredictQuantileBootLM(x = log10(sampledataclean_pre$Area), bootdf = massmodelcurverange, minormax = "min")) * 10^(mean(massmodel$residuals^2)/2)
sampledataclean_pre$particlemass_max <- 10^(PredictQuantileBootLM(x = log10(sampledataclean_pre$Area), bootdf = massmodelcurverange, minormax = "max")) * 10^(mean(massmodel$residuals^2)/2)


RunoffSamples <- c("Santa Ana 2-2-19 3 32pm 10 mins 2 2 181.63g HH", "Santa Ana 2 2 3 51pm 10min", "Santa Ana 2-2-19 4 29pm 2 mins 277.2g HH", "Santa Ana 2 2 54 seconds 5 03 720.25 grams", "Santa Ana 1-17-19 7 26 3 min HH", "Santa Ana 1-17-19 6 30pm 1 min 27 sec ST HH", "Santa Ana 1-17-19 5 00pm 30 sec ST", "Santa Ana 1-17-19 4 29pm 2 min 2 2" )


sampledataclean <- sampledataclean_pre %>%
  mutate(chan_depth_m = chan_depth* 0.3048) %>% #in meters
  mutate(chan_velocity_m = chan_velocity* 0.3048) %>% #in meters
  mutate(chan_depth_m_min = chan_depth_min* 0.3048) %>% #in meters
  mutate(chan_velocity_m_min = chan_velocity_min* 0.3048) %>% #in meters
  mutate(chan_depth_m_max = chan_depth_max* 0.3048) %>% #in meters
  mutate(chan_velocity_m_max = chan_velocity_max* 0.3048) %>% #in meters
  mutate(Sampledepth = ifelse(chan_depth_m < 0.2, chan_depth_m , 0.2)) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  mutate(Sampledepth_min = ifelse(chan_depth_m_min < 0.1, chan_depth_m_min, 0.1)) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  mutate(Sampledepth_max = ifelse(chan_depth_m_max < 0.3, chan_depth_m_max, 0.3)) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  mutate(SampleSize = Sampledepth * 0.4 * chan_velocity_m * Duration..min. * 60) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  mutate(SampleSize_min = Sampledepth_min * 0.4 * chan_velocity_m_min * Duration..min. * 60) %>% #min sample depth set to 0.1 in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  mutate(SampleSize_max = Sampledepth_max * 0.4 * chan_velocity_m_max * Duration..min. * 60) %>% #max sample depth set to 0.3 in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  #mutate(SampleSize_min = ifelse(chan_depth_m_min < 0.4, chan_depth_m_min * 0.4 * chan_velocity_m_min * Duration..min. * 60, 0.5 * 0.4 * chan_velocity_m_min * Duration..min. * 60)) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  #mutate(SampleSize_max = ifelse(chan_depth_max < 0.4, chan_depth_m_max * 0.4 * chan_velocity_m_max * Duration..min. * 60, 0.5 * 0.4 * chan_velocity_m_max * Duration..min. * 60)) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  mutate(particlespersecond = 1 /(Duration..min. * 60)) %>%
  mutate(concentration = 1/SampleSize) %>%
  mutate(shear_velocity = sqrt(9.8*0.003954717*chan_depth_m)) %>% #in meters
  mutate(proportion_sampled = Sampledepth/chan_depth_m) %>%
  mutate(proportion_sampled_max = Sampledepth_max/chan_depth_m_min) %>%#Half submerged.
  mutate(proportion_sampled_max = ifelse(proportion_sampled_max > 1, 1, proportion_sampled_max)) %>%
  mutate(proportion_sampled_min = Sampledepth_min/chan_depth_m_max) %>% #May consider rearranging these probabilities to have max by min so that it amplifies the effect. I think currently we are standardizing the uncertainty based on discharge ranges already. 
  mutate(proportion_sampled_min = ifelse(proportion_sampled_min > 1, 1, proportion_sampled_min)) %>%
  mutate(Runoff = ifelse(SampleName %in% RunoffSamples, "Runoff", "Nonrunoff"))



#Possible Rouse Numbers Measured Rising Velocities ----
#in meters per second
velocities_macorplastic <- c(13.59, 8.60, 16.19, 15.87, 4.21, 5.56, 9.55, 2.99, 6.56, 2.21, 4.08, 4.93, 6.11) / 10
mean(velocities_macorplastic)
median(velocities_macorplastic)
min(velocities_macorplastic)
max(velocities_macorplastic)

RouseNumbers <- expand.grid(rising_vel_m_s = velocities_macorplastic, shear_vel_m_s = unique(sampledataclean$shear_velocity)) %>%
  mutate(RouseNum = rising_vel_m_s/shear_vel_m_s)

min(RouseNumbers$RouseNum)
hist(RouseNumbers$RouseNum[RouseNumbers$RouseNum<2.5])
quantile(RouseNumbers$RouseNum, seq(0,1, by = 0.01))

#Particle size comparison runoff-baseflow ----
ggplot(sampledataclean, aes(sqrt(Area), color = Runoff)) + stat_ecdf(size = 3) + scale_color_viridis_d() + scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c(1, 10, 100, 1000), limits = c(1,1000)) + theme_gray() 

runofffit <- filter(sampledataclean, Runoff == "Runoff")
nonrunofffit <- filter(sampledataclean, Runoff == "Nonrunoff")

#null hypothesis is that they are drawn from the same distribution. 
ks.test(runofffit$Area, nonrunofffit$Area, alternative = "two.sided")
#This demonstrates that pdfs are the same

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
  mutate(min_rouse = find_min_rouse(settling = velocities_macorplastic, shear = shear_velocity)) %>%
  group_by(SampleName, Sampledepth_min, Sampledepth_max, chan_discharge_m, chan_discharge_m_min, chan_discharge_m_max, SampleSize_min, SampleSize_max, proportion_sampled_min, proportion_sampled_max, Mass..g., proportion_sampled, SampleSize, dateTime, Duration..min., Date, chan_depth_m,chan_depth_m_min, chan_depth_m_max, chan_velocity_m_max,  chan_velocity_m_min, chan_velocity_m, shear_velocity, Runoff, min_rouse) %>%
  dplyr::summarise(count = n(), sumarea = sum(Area), summass = sum(particlemass, na.rm = T), minsummass = sum(particlemass_min, na.rm = T), maxsummass = sum(particlemass_max, na.rm = T)) %>%
  mutate(countconcentration = count/SampleSize*proportion_sampled) %>%
  mutate(countconcentration_min = count/SampleSize_max*proportion_sampled_min) %>%
  mutate(countconcentration_max = count/SampleSize_min*proportion_sampled_max) %>%
  mutate(areaconcentration = sumarea/SampleSize*proportion_sampled) %>%
  mutate(areaconcentration_min = sumarea/SampleSize_max*proportion_sampled_min) %>%
  mutate(areaconcentration_max = sumarea/SampleSize_min*proportion_sampled_max) %>%
  mutate(massconcentration = summass/SampleSize*proportion_sampled) %>%
  mutate(massconcentration_min = minsummass/SampleSize_max*proportion_sampled_min) %>%
  mutate(massconcentration_max = maxsummass/SampleSize_min*proportion_sampled_max) %>%
  mutate(measured_mass = Mass..g./SampleSize*proportion_sampled) %>%
  arrange(desc(dateTime))

trueconcentration <- numeric(nrow(totalConcentrationDischarge))
for(n in 1:nrow(totalConcentrationDischarge)){
  if(totalConcentrationDischarge$Sampledepth_min[n] >= totalConcentrationDischarge$chan_depth_m_max[n] | totalConcentrationDischarge$min_rouse[n] > 2.5){
    trueconcentration[n] <- 1
  }
  else{
      trueconcentration[n] <- RouseAverageSurfaceSample(Concentration = 1, RouseNumber = totalConcentrationDischarge$min_rouse[n], SampleDepth = totalConcentrationDischarge$Sampledepth_min[n], StreamDepth = totalConcentrationDischarge$chan_depth_m_max[n])/1
  }
}

totalConcentrationDischarge_2 <- totalConcentrationDischarge %>%
  bind_cols(concentration_multiple = trueconcentration) %>%
  mutate(countconcentration_max_rouse = countconcentration_max * concentration_multiple)

#compare the measured mass of some samples to the particle prediction procedure. 
ggplot(totalConcentrationDischarge, aes(x = summass, y = Mass..g.)) + geom_point() + scale_y_log10() + scale_x_log10() + geom_abline(intercept = 0, slope = 1)

#Concentration discharge relationships ----
#Linear Relatoinship, concentration - discharge

ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = areaconcentration)) +
  geom_smooth(color = "black") + 
  geom_point() + 
  geom_linerange(aes(xmin = chan_discharge_m_min, xmax =chan_discharge_m_max)) + 
  geom_linerange(aes(ymin = areaconcentration_min, ymax =areaconcentration_max)) + 
  #geom_errorbarh(aes(xmin = chan_discharge_m_min, xmax =chan_discharge_m_max)) + 
  scale_color_viridis_c() +
  scale_x_log10() + 
  scale_y_log10()  +
  theme_gray(base_size = 18) + 
  labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Area Concentration ("~mm^2~m^-3~")"))+ 
  coord_fixed()


ggplot(totalConcentrationDischarge_2, aes(x = chan_discharge_m, y = countconcentration)) +
  geom_smooth(color = "black") + 
  geom_linerange(aes(xmin = chan_discharge_m_min, xmax =chan_discharge_m_max)) + 
  geom_linerange(aes(ymin = countconcentration_min, ymax =countconcentration_max)) +
  geom_point() + 
  scale_color_viridis_c() +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_gray(base_size = 18) + 
  labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ 
  coord_fixed()

ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = massconcentration)) +
  geom_smooth(color = "black") + 
  geom_linerange(aes(xmin = chan_discharge_m_min, xmax =chan_discharge_m_max)) + 
  geom_linerange(aes(ymin = massconcentration_min, ymax =massconcentration_max)) +
  geom_point() + 
  scale_color_viridis_c() +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_gray(base_size = 18) + 
  labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ 
  coord_fixed()

ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = areaconcentration/countconcentration)) + geom_point() + geom_smooth(method = "lm", color = "black") #+ scale_x_log10(breaks = c(1,10,100,1000), labels = c(1,10,100,1000), limits = c(1,1000)) + scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100), limits = c(0.01,100)) + theme_gray(base_size = 18) + labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")"))+ coord_fixed()

#Fit area to count concentrations.
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


#May actually not do this after all. Should calculate uncertainties of concentrations and discharge above. 
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

#USGS Data pull ----
#Don't run this code again, just useful for first grab.
#param_cd <- read.csv("Data/param_cd.csv")
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