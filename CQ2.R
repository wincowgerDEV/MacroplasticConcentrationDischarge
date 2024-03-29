#Libraries ----
library(dataRetrieval)
library(dplyr)
library(ggplot2)
library(mgcv)
library(readxl)
library(data.table)
library(stringr)
library(viridis)
library(tidyr)
#library(MASS)
#library(fitdistrplus)
#library(EcoHydRology)
#library(hydrostats)
#library(plotly)
library(matrixStats)
library(tidymv)

#Global settings ----
options(scipen = 999)

#Functions ----
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

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}


#Data ----

#Precipitation
precipitation <- read.csv("Data/RiversideMuniAirport_Cleaned.csv") %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  mutate(precip_mm = na_if(Total.Precip..in., "M")) %>%
  mutate(precip_mm = as.numeric(precip_mm) * 25.4)

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
  mutate(gage_height_m = gage_height_va * 0.3048) %>%
  mutate(chan_discharge_m = chan_discharge * 0.0283168)

#Higher uncertainties when using rating curves for higher flows because few data points. 
#Look at USGS rating curve data
ggplot(measurements, aes(x = gage_height_va, y = chan_depth_m)) + geom_point() + geom_smooth(method = "lm") + scale_y_log10() + scale_x_log10()
ggplot(measurements, aes(x = gage_height_va, y = chan_velocity_m)) + geom_point() + geom_smooth(method = "lm") + scale_y_log10() + scale_x_log10()
ggplot(measurements, aes(x = gage_height_va, y = chan_discharge)) + geom_point() + geom_smooth(method = "lm") + scale_y_log10() + scale_x_log10()

#Particles with sinking removed.
csvfiles_1 <- list.files(path = "Data/Santa Ana River SamplesNoFloatRemoved/", pattern = ".csv", recursive = T, full.names = T)

filelistcsv_1 <- csvfiles_1 %>%
  lapply(read.csv)

csvfiles <- csvfiles_1
filelistcsv <- filelistcsv_1

#Analysis ----

##Plastic Masses to particle size relations ----
Master <- read.csv("Data/ParticleSizeConversionData/Plastic Masses/MasterSheet - Sheet1.csv")

Mast <- Master[complete.cases(Master$Mass),] %>%
  filter(Area > 25) 

nrow(Mast)

massmodel <- lm(log10(Mass) ~ log10(Area), data = Mast)

ggplot(Mast, aes(y = Mass, x =  Area)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_gray_etal() + 
  coord_fixed() + 
  scale_x_log10() + 
  scale_y_log10() +
  labs(x = bquote("Particle Projected Area ("~mm^2~")"), y = "Particle Mass (g)")
  

#By material
ggplot(Mast %>% filter(Material %in% c("Film", "Foam", "Hard Plastic")), aes(y = Mass, x =  Area, color = Material)) + 
    geom_point() +
    geom_smooth(method = "lm") + 
    theme_gray_etal() + 
    coord_fixed() + 
    scale_x_log10() + 
    scale_y_log10() +
    labs(x = bquote("Particle Projected Area ("~mm^2~")"), y = "Particle Mass (g)")

##Flow models ----
dischargecurve <- lm(log10(chan_discharge) ~ log10(gage_height_va), data = measurements)
summary(dischargecurve)

velocitycurve <- lm(log10(chan_velocity) ~ log10(gage_height_va), data = measurements)
summary(velocitycurve)
#velocitycurverange <- BootLM(y = log10(measurements$chan_velocity), x = log10(measurements$gage_height_va))

depthcurve <- lm(log10(chan_area/chan_width) ~ log10(gage_height_va), data = measurements)
summary(depthcurve)

#flow model summaries metric
summary(lm(log10(chan_discharge_m) ~ log10(gage_height_m), data = measurements))

summary(lm(log10(chan_velocity_m) ~ log10(gage_height_m), data = measurements))

summary(lm(log10(chan_depth_m) ~ log10(gage_height_m), data = measurements))

#Discharge figure ----
Discharge <- Q %>%
  mutate(dateTime = as.POSIXct(strptime(dateTime, format = "%Y-%m-%d %H:%M:%S"), tz = "America/Los_Angeles")) %>%
  mutate(DEP = 10^(log10(X_00065_00000)*dischargecurve$coefficients[2] + dischargecurve$coefficients[1]) * 10^(mean(dischargecurve$residuals^2)/2))

daily_flow <- Discharge %>%
  mutate(Date = as.Date(dateTime)) %>%
  group_by(Date) %>%
  summarise(daily_flow_m = sum(DEP * 0.0283168 * 86400), daily_flow_m_mean = mean(DEP * 0.0283168))

sample_dates <- as.Date(c(as.POSIXct("2019-02-02 12:00:00", tz="America/Los_Angeles"), 
                  as.POSIXct("2019-01-19 12:00:00", tz="America/Los_Angeles"), 
                  as.POSIXct("2019-08-09 12:00:00", tz="America/Los_Angeles"), 
                  as.POSIXct("2019-10-21 12:00:00", tz="America/Los_Angeles"), 
                  as.POSIXct("2019-2-13 12:00:00", tz="America/Los_Angeles")))
  
ggplot()+ 
  geom_line(data = daily_flow, aes(y = daily_flow_m_mean, x = Date), size = 1) + 
  #geom_col(data = precipitation, aes(x= Date, y = precip_mm), fill = "blue") + 
  scale_x_date(limits = as.Date(c(as.POSIXct("2018-10-01 00:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-10-30 23:59:00", tz="America/Los_Angeles")))) + 
  labs(y = "Mean Daily Discharge ("~m^3~s^-1~")", x = "Date") + 
  geom_vline(xintercept = sample_dates) + 
  scale_y_log10(limits = c(0.1, 1000)) + 
  theme_gray_etal()

ggplot()+ 
  geom_line(data = daily_flow, aes(y = daily_flow_m_mean, x = Date), size = 1) + 
  geom_line(data = precipitation, aes(x= Date, y = precip_mm), color = "blue", size = 1, alpha = 0.3) + 
  scale_x_date(limits = as.Date(c(as.POSIXct("2018-10-01 00:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-10-30 23:59:00", tz="America/Los_Angeles")))) + 
  labs(y = "Daily Precipitation (mm)", x = "Date") + 
  #geom_vline(xintercept = sample_dates) + 
  scale_y_log10(limits = c(0.1, 1000)) + 
  theme_gray_etal()

##Exceedence probabilities
ggplot() + stat_ecdf(aes(x = Discharge$DEP * 0.0283168)) + scale_x_log10() + theme_gray_etal()

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
  mutate(class = "Rising") %>%
  rbind(dplyr::select(NoFloat, Area, SampleName) %>%
          dplyr::filter(!SampleName %in% RemoveThese) %>%
          dplyr::filter(Area > 1) %>%
          mutate(class = "Settling") %>%
          dplyr::select(Area, class))

ggplot(sizecompare, aes(x = sqrt(Area), y = class)) + 
  geom_violin(aes(fill = class)) +
  geom_text(aes(x = 400, y = "Rising", label = "N = 927"), size = 10)+
  geom_text(aes(x = 400, y = "Settling", label = "N = 17"), size = 10)+
  scale_fill_viridis_d() + 
  geom_boxplot(aes(fill = class), color = "black", width=.2, notch = T) + 
  scale_x_log10(limits = c(1, 1000), breaks = c(1, 10,100,1000)) +
  theme_gray_etal(base_size = 20, bgcolor = "gray") +
  labs(x = "Nominal Particle Size (mm)", y = "")
#Double checked this, these are only for samples which are included in this study for concentration discharge relationships. 
sizecompare %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(count = n())

#propotion rising
927/(927+17)

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
  mutate(chan_discharge_m = chan_discharge * 0.0283168) %>%
  right_join(MasterList) 

#hydrograph Feb 2
sampledataclean_pre %>%
  dplyr::select(chan_discharge_m, SampleName, dateTime) %>%
  distinct() %>%
  ggplot() + 
  geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1)+ 
  geom_point(aes(y = chan_discharge_m, x = dateTime), color = "red", size = 3, shape = 23, fill = "red") + 
  scale_x_datetime(limits = c(as.POSIXct("2019-02-02 01:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-02-03 12:00:00", tz="America/Los_Angeles"))) + 
  ylim(0,1000) + 
  scale_y_log10() +
  labs(y = "Discharge (cms)", x = "Date Time") + 
  theme_gray_etal()

#Hydrograph Jan 17
sampledataclean_pre %>%
  dplyr::select(chan_discharge_m, SampleName, dateTime) %>%
  distinct() %>%
  ggplot() + 
  geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1)+ 
  geom_point(aes(y = chan_discharge_m, x = dateTime), color = "red", size = 3, shape = 23, fill = "red") + 
  scale_x_datetime(limits = c(as.POSIXct("2019-01-17 01:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-01-18 12:00:00", tz="America/Los_Angeles"))) + 
  ylim(0,1000) + 
  scale_y_log10() +
  labs(y = "Discharge (cms)", x = "Date Time") + 
  theme_gray_etal()

#Hydrograph lowflow
sampledataclean_pre %>%
  dplyr::select(chan_discharge_m, SampleName, dateTime) %>%
  distinct() %>%
  ggplot() + geom_line(data = Discharge, aes(y = DEP * 0.0283168, x = dateTime), size = 1)+ geom_point(aes(y = chan_discharge_m, x = dateTime), color = "red", size = 3, shape = 23, fill = "red") + scale_x_datetime(limits = c(as.POSIXct("2019-02-13 12:00:00", tz="America/Los_Angeles"),as.POSIXct("2019-02-14 01:00:00", tz="America/Los_Angeles"))) + ylim(0,10) + labs(y = "Discharge (cms)", x = "Date Time") + theme_gray_etal()

#But don't these require the log10 of chan_discharge? Not sure these are working right. 
sampledataclean_pre$chan_velocity <- 10^(log10(sampledataclean_pre$gage_height)*velocitycurve$coefficients[2] + velocitycurve$coefficients[1]) * 10^(mean(velocitycurve$residuals^2)/2)

sampledataclean_pre$chan_depth <- 10^(log10(sampledataclean_pre$gage_height)*depthcurve$coefficients[2] + depthcurve$coefficients[1]) * 10^(mean(depthcurve$residuals^2)/2)

sampledataclean_pre$particlemass <- 10^(log10(sampledataclean_pre$Area)*massmodel$coefficients[2] + massmodel$coefficients[1]) * 10^(mean(massmodel$residuals^2)/2)

RunoffSamples <- c("Santa Ana 2-2-19 3 32pm 10 mins 2 2 181.63g HH", "Santa Ana 2 2 3 51pm 10min", "Santa Ana 2-2-19 4 29pm 2 mins 277.2g HH", "Santa Ana 2 2 54 seconds 5 03 720.25 grams", "Santa Ana 1-17-19 7 26 3 min HH", "Santa Ana 1-17-19 6 30pm 1 min 27 sec ST HH", "Santa Ana 1-17-19 5 00pm 30 sec ST", "Santa Ana 1-17-19 4 29pm 2 min 2 2" )

sampledataclean <- sampledataclean_pre %>%
  mutate(chan_depth_m = chan_depth* 0.3048) %>% #in meters
  mutate(chan_velocity_m = chan_velocity* 0.3048) %>% #in meters
  mutate(Sampledepth = ifelse(chan_depth_m < 0.2, chan_depth_m , 0.2)) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  mutate(SampleSize = Sampledepth * 0.4 * chan_velocity_m * Duration..min. * 60) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
  mutate(particlespersecond = 1 /(Duration..min. * 60)) %>%
  mutate(concentration = 1/SampleSize) %>%
  mutate(shear_velocity = sqrt(9.8*0.003954717*chan_depth_m)) %>% #in meters
  mutate(proportion_sampled = Sampledepth/chan_depth_m) %>%
  mutate(Runoff = ifelse(SampleName %in% RunoffSamples, "Stormflow", "Lowflow"))


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
ggplot(sampledataclean, aes(sqrt(Area), color = Runoff)) + 
  stat_ecdf(size = 3, alpha = 0.8) + 
  scale_color_viridis_d(option = "E") + 
  scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c(1, 10, 100, 1000), limits = c(1,1000)) + 
  theme_gray_etal(base_size = 20, bgcolor = "gray") + 
  labs(x = "Nominal Particle Size (mm)", y = "Percent Smaller")

runofffit <- filter(sampledataclean, Runoff == "Stormflow")
filter(sampledataclean, Runoff == "Stormflow") %>%
  nrow()

filter(sampledataclean, Runoff == "Stormflow") %>%
  distinct(SampleName)

nonrunofffit <- filter(sampledataclean, Runoff == "Lowflow")
filter(sampledataclean, Runoff == "Lowflow") %>%
  nrow()

filter(sampledataclean, Runoff == "Lowflow") %>%
  distinct(SampleName)


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
  group_by(SampleName, chan_discharge_m, Mass..g., proportion_sampled, SampleSize, dateTime, Duration..min., Date, chan_depth_m, chan_velocity_m, shear_velocity, Runoff, min_rouse) %>%
  dplyr::summarise(count = n(), sumarea = sum(Area), summass = sum(particlemass, na.rm = T)) %>%
  mutate(countconcentration = count/SampleSize*proportion_sampled) %>%
  mutate(areaconcentration = sumarea/SampleSize*proportion_sampled) %>%
  mutate(massconcentration = summass/SampleSize*proportion_sampled) %>%
  mutate(measured_mass = Mass..g./SampleSize*proportion_sampled) %>%
  arrange(desc(dateTime))

totalConcentrationDischarge %>%
  filter(Runoff == "Lowflow") %>%
  pull(chan_depth_m) %>%
  mean()

mean(totalConcentrationDischarge$summass - totalConcentrationDischarge$Mass..g., na.rm = T)
#mpare the measured mass of some samples to the particle prediction procedure. 
ggplot(totalConcentrationDischarge, aes(y = summass, x = Mass..g.)) + 
  geom_point() + 
  scale_y_log10() + 
  scale_x_log10() + 
  #geom_text(aes(x = 0.01, y = 10, label = "Mean Residual = -4.37 g"))
  geom_abline(intercept = 0, slope = 1, size = 1) + 
  theme_gray_etal(base_size = 15) + 
  coord_fixed() +
  labs(x = "Sample Mass Observed (g)", y = "Sample Mass Estimated (g)")


#Fit area to count concentrations.
ggplot(totalConcentrationDischarge, aes(x = areaconcentration, y = countconcentration)) +
  geom_smooth(method = "lm", color = "black") + 
  geom_point() + 
  scale_color_viridis_c() +
  scale_x_log10() + 
  scale_y_log10()  +
  theme_gray_etal(base_size = 18)

#Estimate Flux ----
#Uncertaintys of confidence intervals, if variables are correlated, just sum. https://stats.stackexchange.com/questions/223924/how-to-add-up-partial-confidence-intervals-to-create-a-total-confidence-interval
#https://stats.stackexchange.com/questions/305382/how-do-i-calculate-the-confidence-interval-for-the-product-of-two-numbers-with
#Constant Mean
Flux <- Discharge %>%
  filter(dateTime > as.POSIXct("2018-10-01 00:00:00", tz="America/Los_Angeles") & dateTime < as.POSIXct("2019-09-30 23:59:00", tz="America/Los_Angeles"))

Flux$chan_discharge_m <- 10^(log10(Flux$X_00065_00000)*dischargecurve$coefficients[2] + dischargecurve$coefficients[1]) * 10^(mean(dischargecurve$residuals^2)/2)* 0.0283168

constant_mean_metric_tonnes <- sum(Flux$chan_discharge_m * mean(totalConcentrationDischarge$massconcentration) * 15 * 60, na.rm = T)/10^6

#Flux estimate using regression
#Flow duration curve
ggplot() + stat_ecdf(aes(x = Flux$chan_discharge_m)) + scale_x_log10() + theme_gray_etal()
median(Flux$chan_discharge_m)
max(Flux$chan_discharge_m)
min(Flux$chan_discharge_m)

#https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
hist(log10(totalConcentrationDischarge$massconcentration))
hist(log10(totalConcentrationDischarge$chan_discharge_m))
10^mean(log10(totalConcentrationDischarge$chan_discharge_m))
mean(totalConcentrationDischarge$chan_discharge_m)
#geometric mean of discharge seems to be a better fit to the central tendancy of the discharge time series. 

shapiro.test(log10(totalConcentrationDischarge$chan_discharge_m))
shapiro.test(log10(totalConcentrationDischarge$massconcentration))
#Mass concentration is normally distributed but discharge is not. 
shapiro.test(log10(Flux$chan_discharge_m))
hist(log10(Flux$chan_discharge_m))
10^mean(log10(Flux$chan_discharge_m))
mean(Flux$chan_discharge_m)

totalConcentrationDischarge$chan_discharge_m_log10 <- log10(totalConcentrationDischarge$chan_discharge_m)
totalConcentrationDischarge$massconcentration_log10 <- log10(totalConcentrationDischarge$massconcentration)

Flux$chan_discharge_m_log10 <- log10(Flux$chan_discharge_m)

#Concentration discharge relationships ----
#count concentration model.
model <- gam(log10(countconcentration)~s(log10(chan_discharge_m), k = 7), data = totalConcentrationDischarge, method = "REML")
10^(mean(model$residuals^2)/2)
summary.gam(model)
plot(model, shade = TRUE,  seWithMean = TRUE, residuals = TRUE, pch = 16, cex = 0.8)

#mass concentration model.
mass.model <- gam(massconcentration_log10~s(chan_discharge_m_log10, k = 7), data = totalConcentrationDischarge, method = "REML")
summary.gam(mass.model)
10^(mean(mass.model$residuals^2)/2)
plot(mass.model, shade = TRUE,  seWithMean = TRUE, residuals = TRUE, pch = 16, cex = 0.8)
model_p <- predict_gam(mass.model)

#add in mean concentrations and uncertainties if high or low
flux_low = Flux$chan_discharge_m < min(totalConcentrationDischarge$chan_discharge_m)
#table(flux_low)
flux_high = Flux$chan_discharge_m > max(totalConcentrationDischarge$chan_discharge_m)

arranged <- totalConcentrationDischarge %>%
  arrange(chan_discharge_m)

mean_high_flux_concentration <- mean(arranged$massconcentration[17:20])
mean_low_flux_concentration <- mean(arranged$massconcentration[1:4])

#center
flux_pred <- predict(mass.model, Flux, se.fit = TRUE)
flux_fit_concentration <- 10^(flux_pred$fit)* 10^(mean(mass.model$residuals^2)/2)
Flux$flux_fit_metric_tonnes_time <- ifelse(flux_low, mean_low_flux_concentration, ifelse(flux_high, mean_high_flux_concentration, flux_fit_concentration)) * Flux$chan_discharge_m * 15 * 60/10^6
flux_fit_metric_tonnes <- sum(ifelse(flux_low, mean_low_flux_concentration, ifelse(flux_high, mean_high_flux_concentration, flux_fit_concentration))* Flux$chan_discharge_m* 15 * 60)/10^6


#Bootstrap all datasets 10k times to get fluxes---- 
#Pre allocation
mean_flux <- numeric(length = 10000)
regression_flux <- numeric(length = 10000)

sample_discharge_boot <- matrix(data=NA,nrow=nrow(totalConcentrationDischarge),ncol=10000)
sample_countconcentration_boot <- matrix(data=NA,nrow=nrow(totalConcentrationDischarge),ncol=10000)
sample_areaconcentration_boot <- matrix(data=NA,nrow=nrow(totalConcentrationDischarge),ncol=10000)
sample_massconcentration_boot <- matrix(data=NA,nrow=nrow(totalConcentrationDischarge),ncol=10000)

Flux_boot <- Flux

sampledata_boot <- sampledata
sampledataclean_pre_boot <- sampledata_boot %>%
  mutate_if(is.factor, as.character) %>%
  mutate(dateTime = as.POSIXct(strptime(paste(Date, Start.Time, sep = " "), format = "%m/%d/%Y %H:%M"), tz = "America/Los_Angeles")) %>%
  rename(SampleName = Sample.ID) %>%
  mutate(gage_height = approx(Discharge$dateTime, Discharge$X_00065_00000, xout = dateTime, rule = 2, method = "linear", ties=mean)[[2]]) %>%
  right_join(MasterList) 

set.seed(777)

for(n in 1:length(mean_flux)){
  #bootstrapped datasets
  Mast_boot <- Mast[sample(1:nrow(Mast), size = nrow(Mast), replace = T),]
  measurements_boot <- measurements[sample(1:nrow(measurements), size = nrow(measurements), replace = T),]
  
  #Flow models 
  dischargecurve_boot <- lm(log10(chan_discharge) ~ log10(gage_height_va), data = measurements_boot)
  velocitycurve_boot <- lm(log10(chan_velocity) ~ log10(gage_height_va), data = measurements_boot)
  depthcurve_boot <- lm(log10(chan_area/chan_width) ~ log10(gage_height_va), data = measurements_boot)
  massmodel_boot <- lm(log10(Mass) ~ log10(Area), data = Mast_boot)

  
#Data cleanup  
  sampledataclean_pre_boot$chan_discharge <- 10^(log10(sampledataclean_pre_boot$gage_height)*dischargecurve_boot$coefficients[2] + dischargecurve_boot$coefficients[1]) * 10^(mean(dischargecurve_boot$residuals^2)/2)
  sampledataclean_pre_boot$chan_velocity <- 10^(log10(sampledataclean_pre_boot$gage_height)*velocitycurve_boot$coefficients[2] + velocitycurve_boot$coefficients[1]) * 10^(mean(velocitycurve_boot$residuals^2)/2)
  sampledataclean_pre_boot$chan_depth <- 10^(log10(sampledataclean_pre_boot$gage_height)*depthcurve_boot$coefficients[2] + depthcurve_boot$coefficients[1]) * 10^(mean(depthcurve_boot$residuals^2)/2)
  sampledataclean_pre_boot$particlemass <- 10^(log10(sampledataclean_pre_boot$Area)*massmodel_boot$coefficients[2] + massmodel_boot$coefficients[1]) * 10^(mean(massmodel_boot$residuals^2)/2)
  Flux_boot$chan_discharge_m <- 10^(log10(Flux_boot$X_00065_00000)*dischargecurve_boot$coefficients[2] + dischargecurve_boot$coefficients[1]) * 10^(mean(dischargecurve_boot$residuals^2)/2)* 0.0283168
  Flux_boot$chan_discharge_m_log10 <- log10(Flux_boot$chan_discharge_m)
  
  sample_depth_boot <- runif(1, min = 0.1, max = 0.3) #sample depth range
  
  sampledataclean_boot <- sampledataclean_pre_boot %>%
    mutate(chan_discharge_m = chan_discharge* 0.0283168) %>% #in meters
    mutate(chan_depth_m = chan_depth* 0.3048) %>% #in meters
    mutate(chan_velocity_m = chan_velocity* 0.3048) %>% #in meters
    mutate(Sampledepth = ifelse(chan_depth_m < sample_depth_boot, chan_depth_m, sample_depth_boot)) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
    mutate(SampleSize = Sampledepth * 0.4 * chan_velocity_m * Duration..min. * 60) %>% #in cubic meters, assumes sample net is 1/3rd submerged when not sitting on bottom. 
    mutate(shear_velocity = sqrt(9.8*0.003954717*chan_depth_m)) %>% #in meters
    mutate(proportion_sampled = Sampledepth/chan_depth_m) 
  
  sample_amount_boot <- runif(1, min = -0.1, max = 0.1) #Multiple on sample
  
  #Corrected concentrations, depth integrated. 
  totalConcentrationDischarge_boot <- sampledataclean_boot %>%
    group_by(SampleName, chan_discharge_m, Mass..g., proportion_sampled, SampleSize, dateTime, Duration..min., Date, chan_depth_m, chan_velocity_m, shear_velocity) %>%
    dplyr::summarise(count = n(), sumarea = sum(Area), summass = sum(particlemass, na.rm = T)) %>%
    mutate(countconcentration = (count + count*sample_amount_boot)/SampleSize*proportion_sampled) %>%
    mutate(areaconcentration = (sumarea+ sumarea*sample_amount_boot)/SampleSize*proportion_sampled) %>%
    mutate(massconcentration = (summass + summass*sample_amount_boot)/SampleSize*proportion_sampled) %>%
    mutate(measured_mass = Mass..g./SampleSize*proportion_sampled) %>%
    arrange(desc(dateTime))
  
  sample_discharge_boot[,n] <- totalConcentrationDischarge_boot$chan_discharge_m
  sample_countconcentration_boot[,n] <- totalConcentrationDischarge_boot$countconcentration
  sample_areaconcentration_boot[,n] <- totalConcentrationDischarge_boot$areaconcentration
  sample_massconcentration_boot[,n] <- totalConcentrationDischarge_boot$massconcentration
  
  #Bootstrap the end concentrations
  totalConcentrationDischarge_boot <- totalConcentrationDischarge_boot[sample(1:nrow(totalConcentrationDischarge_boot), size = nrow(totalConcentrationDischarge_boot), replace = T),]
    
  #Estimate Flux 
  
  mean_flux[n] <- sum(Flux_boot$chan_discharge_m * mean(totalConcentrationDischarge_boot$massconcentration) * 15 * 60, na.rm = T)/10^6
  
  #Flux estimate using regression
  
  totalConcentrationDischarge_boot$chan_discharge_m_log10 <- log10(totalConcentrationDischarge_boot$chan_discharge_m)
  totalConcentrationDischarge_boot$massconcentration_log10 <- log10(totalConcentrationDischarge_boot$massconcentration)

  mass.model <- gam(massconcentration_log10~s(chan_discharge_m_log10, k = 7), data = totalConcentrationDischarge_boot, method = "REML")
  
  #add in mean concentrations and uncertainties if high or low
  flux_low = Flux_boot$chan_discharge_m < min(totalConcentrationDischarge_boot$chan_discharge_m)
  flux_high = Flux_boot$chan_discharge_m > max(totalConcentrationDischarge_boot$chan_discharge_m)
  
  arranged <- totalConcentrationDischarge_boot %>%
    arrange(chan_discharge_m)
    
  mean_high_flux_concentration <- mean(arranged$massconcentration[17:20])
  mean_low_flux_concentration <- mean(arranged$massconcentration[1:4])

  #center
  flux_pred <- predict(mass.model, Flux, se.fit = TRUE)
  flux_fit_concentration <- 10^(flux_pred$fit)* 10^(mean(mass.model$residuals^2)/2)
  Flux_boot$flux_fit_metric_tonnes_time <- ifelse(flux_low, mean_low_flux_concentration, ifelse(flux_high, mean_high_flux_concentration, flux_fit_concentration)) * Flux_boot$chan_discharge_m * 15 * 60/10^6
  
  regression_flux[n] <- sum(Flux_boot$flux_fit_metric_tonnes_time)
}

figure_table <- tibble(
  annual_flux_tonnes = c(flux_fit_metric_tonnes, constant_mean_metric_tonnes),
  min_flux_tonnes    = c(quantile(regression_flux, c(0.025)), quantile(mean_flux, c(0.025))),
  max_flux_tonnes    = c(quantile(regression_flux, c(0.975)), quantile(mean_flux, c(0.975))),
  name               = c("Generalized Additive Model", "Constant Mean")
)

#if I fit the geometric mean to the constant mean model the predictions are almost identical. 
ggplot(figure_table, aes(y = name, x = annual_flux_tonnes)) + 
  geom_errorbar(aes(xmin = min_flux_tonnes, xmax = max_flux_tonnes)) + 
  geom_point() + 
  scale_x_log10(breaks = c(1, 10, 100, 1000), limits = c(1,1000)) + 
  theme_gray_etal() + 
  labs(y = "", x = "Annual Flux (metric tonnes)")

totalConcentrationDischarge$discharge_min <- as.vector(rowQuantiles(sample_discharge_boot, probs = c(0.025)))
totalConcentrationDischarge$discharge_max <- as.vector(rowQuantiles(sample_discharge_boot, probs = c(0.975)))
totalConcentrationDischarge$countconcentration_min <- as.vector(rowQuantiles(sample_countconcentration_boot, probs = c(0.025)))
totalConcentrationDischarge$countconcentration_max <- as.vector(rowQuantiles(sample_countconcentration_boot, probs = c(0.975)))
totalConcentrationDischarge$massconcentration_min <- as.vector(rowQuantiles(sample_massconcentration_boot, probs = c(0.025)))
totalConcentrationDischarge$massconcentration_max <- as.vector(rowQuantiles(sample_massconcentration_boot, probs = c(0.975)))

write.csv(totalConcentrationDischarge, "totalConcentrationDischarge_booted.csv")

ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = countconcentration)) + 
  geom_path(aes(color = Date), size = 2) + 
  geom_linerange(aes(xmin = discharge_min, xmax =discharge_max), color = "black", alpha = 0.5) + 
  geom_linerange(aes(ymin = countconcentration_min, ymax =countconcentration_max), color = "black", alpha = 0.5) + 
  geom_point(color = "black", alpha = 0.5, size = 2) + 
  scale_color_viridis_d() + 
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), limits = c(1, 10000)) + 
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100), limits = c(0.001,100)) + 
  theme_gray_etal(base_size = 18) + 
  labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")")) + 
  coord_fixed()#+ geom_text(aes(label = SampleName))

ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = countconcentration)) + 
  geom_smooth(color = "black", method = "gam", formula = y ~ s(x, bs = "cs", k = 7)) +
  geom_linerange(aes(xmin = discharge_min, xmax =discharge_max)) + 
  geom_linerange(aes(ymin = countconcentration_min, ymax =countconcentration_max)) + 
  geom_point(color = "black", alpha = 0.5, size = 2) + 
  scale_color_viridis_d() + 
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), limits = c(1, 10000)) + 
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000), limits = c(0.001,1000)) + 
  theme_gray_etal(base_size = 16) + 
  labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Count Concentration ("~num^1~m^-3~")")) + 
  coord_fixed()#+ geom_text(aes(label = SampleName))

ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = massconcentration)) + 
  geom_smooth(color = "black", method = "gam", formula = y ~ s(x, bs = "cs", k = 7)) +
  geom_linerange(aes(xmin = discharge_min, xmax =discharge_max)) + 
  geom_linerange(aes(ymin = massconcentration_min, ymax = massconcentration_max)) + 
  geom_point(color = "black", alpha = 0.5, size = 2) + 
  scale_color_viridis_d() + 
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), limits = c(1, 10000)) + 
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), limits = c(0.0001,100)) + 
  theme_gray_etal(base_size = 16) + 
  labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Mass Concentration ("~g^1~m^-3~")")) + 
  coord_fixed()#+ geom_text(aes(label = SampleName))


ggplot(totalConcentrationDischarge, aes(x = chan_discharge_m, y = massconcentration)) + 
  geom_path(aes(color = Date), size = 2) + 
  geom_point(color = "black", alpha = 0.5, size = 2) + 
  geom_linerange(aes(xmin = discharge_min, xmax =discharge_max)) + 
  geom_linerange(aes(ymin = massconcentration_min, ymax =massconcentration_max)) + 
  scale_color_viridis_d() + 
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), limits = c(1, 10000)) + 
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10), limits = c(0.0001,10)) + 
  theme_gray_etal(base_size = 18) + 
  labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = bquote("Mass Concentration ("~g^1~m^-3~")")) + 
  coord_fixed()#+ geom_text(aes(label = SampleName))


#Effective Discharge 2018 ----
max(Flux$chan_discharge_m)
min(Flux$chan_discharge_m)
bins <- seq(from = min(Flux$chan_discharge_m), to = max(Flux$chan_discharge_m), length.out = 10)
Flux$categories <- cut(Flux$chan_discharge_m, breaks = bins)

Flux_summarized <- Flux %>%
  mutate(Plasticflux = chan_discharge_m * mean(totalConcentrationDischarge$massconcentration) * 60 * 15) %>%
  arrange(chan_discharge_m) %>%
  mutate(cummulative_flux_mean = cumsum(Plasticflux)/sum(Plasticflux)) %>%
  mutate(cummulative_flux_regression = cumsum(flux_fit_metric_tonnes_time)/sum(flux_fit_metric_tonnes_time))

ggplot(Flux_summarized) + 
  geom_line(aes(x = chan_discharge_m, y = cummulative_flux_mean), size = 4) + 
  geom_line(aes(x = chan_discharge_m, y = cummulative_flux_regression), color = "red", size = 4)+
  scale_x_log10() + 
  theme_gray_etal(base_size = 15) + 
  labs(x = "Discharge (cms)", y = "Cummulative Flux")

mean(regression_flux)
hist(log10(regression_flux))
mean(mean_flux)
hist(log10(mean_flux))
quantile(regression_flux, c(0.025, 0.5, 0.975))
quantile(mean_flux, c(0.025, 0.5, 0.975))



#USGS Data pull ----
#Don't run this code again, just useful for first grab.
param_cd <- read.csv("Data/param_cd.csv")
sites <- c("11066460")
#parameterCd <- "00065"
parameterCd_dv <- "00060" #Daily discharge

startDate <- "1990-10-01"
endDate <- "2020-09-30"

availabledata <- whatNWISdata(siteNumber = sites)
availabledata <- availabledata %>%
    mutate(parm_cd = as.numeric(parm_cd)) %>%
    left_join(param_cd)
ggplot(availabledata, aes(x = begin_date)) + geom_histogram()

#Q <- readNWISuv(sites, parameterCd, startDate, endDate, tz = "America/Los_Angeles")
Q_dv <- readNWISdv(sites, parameterCd_dv, startDate, endDate, statCd = "00003")

#ggplot(Q_dv) + stat_ecdf(aes(x = X_00060_00003* 0.0283168)) + scale_x_log10()

water_Flux_summarized <- Q_dv %>%
    mutate(discharge_m = X_00060_00003* 0.0283168) %>%
    arrange(discharge_m) %>%
    mutate(cummulative_flux = cumsum(discharge_m * 86400)/sum(discharge_m * 86400)) 

ggplot(water_Flux_summarized) + 
    geom_line(aes(x = discharge_m, y = cummulative_flux), size = 4) + 
    #geom_line(aes(x = chan_discharge_m, y = cummulative_flux_regression), color = "red", size = 4)+
    scale_x_log10(limits = c(0.1, 1000)) + 
    theme_gray_etal() + 
    labs(x = bquote("Discharge ("~m^3~s^-1~")"), y = "Cumulative Water Volume (% Total)")

#Low Flows 
water_Flux_summarized %>%
    filter(discharge_m < 3) %>%
    pull(discharge_m) %>%
    mean()

#High Flows
water_Flux_summarized %>%
    filter(discharge_m > 3) %>%
    pull(discharge_m) %>%
    mean()

#2 year flow
years <- water_Flux_summarized %>%
    mutate(Year = format(Date, format = "%Y")) %>%
    group_by(Year) %>%
    summarize(maxyearlyflow = max(discharge_m)) %>%
    ungroup() %>%
    mutate(recurrance = (length(Year) + 1)/rank(1/maxyearlyflow)) %>%
    arrange(recurrance)

#max(Q_dv$X_00060_00003* 0.0283168)
#min(Q_dv$X_00060_00003* 0.0283168)
#bins <- seq(from = min(Q_dv$X_00060_00003* 0.0283168), to = max(Q_dv$X_00060_00003* 0.0283168), length.out = 100000)
#Q_dv$categories <- cut(Q_dv$X_00060_00003* 0.0283168, breaks = bins)

#Q_dv_summarized <- Q_dv %>%
#  mutate(flow = X_00060_00003* 0.0283168 *86400) %>%
#  group_by(categories) %>%
#  summarize(sum = sum(flow))

#ggplot(Q_dv_summarized) + geom_col(aes(x = categories, y = sum))
#setwd("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/CQRelationships/Data/Raw Data")
#write.csv(Q, "SiteQ11066460.csv")

#ratingdata <- readNWISrating(sites, "base")
#measurements <- readNWISmeas(sites, expanded = T, tz = "America/Los_Angeles")
#write.csv(measurements, "measurements.csv")
