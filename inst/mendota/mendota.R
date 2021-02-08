rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(gridExtra)
library(ggplot2)
library(pracma)
library(tidyverse)
library(RColorBrewer)
library(zoo)
library(deSolve)
library(LakeMetabolizer)

library(dplyr)
library(readr)
library(LakeEnsemblR)
library(gotmtools)
library(lubridate)

library(thermod)

### GETTING CONFIGURATION INPUT FROM LER YEAML FILE
config_file <- 'LakeEnsemblR.yaml'
folder <- '.'
yaml  <-  file.path(folder, config_file)
lev <- get_yaml_value(config_file, "location", "elevation")
max_depth <- get_yaml_value(config_file, "location", "depth")
hyp_file <- get_yaml_value(config_file, "location", "hypsograph")
hyp <- read.csv(hyp_file)
start_date <- get_yaml_value(config_file, "time", "start")
stop_date <- get_yaml_value(config_file, "time", "stop")
timestep <- get_yaml_value(config_file, "time", "time_step")
bsn_len <- get_yaml_value(config_file, "model_parameters", "bsn_len") 
bsn_wid <- get_yaml_value(config_file, "model_parameters", "bsn_wid") 
therm_depth <- 10 ^ (0.336 * log10(max(bsn_len, bsn_wid)) - 0.245)
simple_therm_depth <- round(therm_depth)

ix <- match(simple_therm_depth, hyp$Depth_meter)

Ve <- trapz(hyp$Depth_meter[1:ix], hyp$Area_meterSquared[1:ix]) * 1e6 # epilimnion volume
Vh <- trapz(hyp$Depth_meter[ix:nrow(hyp)], hyp$Area_meterSquared[ix:nrow(hyp)]) * 1e6 # hypolimnion volume
At <- approx(hyp$Depth_meter, hyp$Area_meterSquared, therm_depth)$y * 1e4 # thermocline area
Ht <- 3 * 100 # thermocline thickness
As <- max(hyp$Area_meterSquared) * 1e4 # surface area
Tin <- 0 # inflow water temperature
Q <- 0 # inflow discharge
Rl <- 0.03 # reflection coefficient (generally small, 0.03)
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <- 11.7 * 10^(-8) # cal / (cm2 d K4) or: 4.9 * 10^(-3) # Stefan-Boltzmann constant in [J (m2 d K4)-1]
eps <- 0.97 # emissivity of water
rho <- 0.9982 # density (g per cm3)
cp <- 0.99 # specific heat (cal per gram per deg C)
c1 <- 0.47 # Bowen's coefficient
a <- 7 # constant
c <- 9e4 # empirical constant
g <- 9.81  # gravity (m/s2)
calParam <- 1 # parameter for calibrating the entrainment over the thermocline depth

### CONVERTING METEO CSV FILE TO THERMOD SPECIFIC DATA
# met <- read.csv('LakeEnsemblR_meteo_standard.csv', stringsAsFactors = F)
# met[, 1] <- as.POSIXct(met[, 1])
# start_date = '1995-01-01 00:00:00'
# stop_date = '2015-12-31 00:00:00'
# tstep <- diff(as.numeric(met[, 1]))
# met$Dewpoint_Air_Temperature_Celsius <- met$Air_Temperature_celsius - ((100 - met$Relative_Humidity_percent)/5.)
# bound <- met %>%
#   dplyr::select(datetime, Shortwave_Radiation_Downwelling_wattPerMeterSquared, 
#                 Air_Temperature_celsius, Dewpoint_Air_Temperature_Celsius, Ten_Meter_Elevation_Wind_Speed_meterPerSecond)
# bound$Shortwave_Radiation_Downwelling_wattPerMeterSquared <- bound$Shortwave_Radiation_Downwelling_wattPerMeterSquared * 2.0E-5 *3600*24
# bound <- bound %>%
#   dplyr::filter(datetime >= start_date & datetime <= stop_date)
# bound = bound %>%
#   group_by(datetime) %>%
#   summarise_all(mean)
# bound$datetime <- seq(1, nrow(bound))
# write_delim(bound, path = paste0('meteo.txt'), delim = '\t')

# load in the boundary data
bound <-read_delim(paste0('meteo.txt'), delim = '\t')

parameters <- c(Ve, Vh, At, Ht, As, Tin, Q, Rl, Acoeff, sigma, eps, rho, cp, c1, a, c, g, simple_therm_depth, calParam)

colnames(bound) <- c('Day','Jsw','Tair','Dew','vW')

# function to calculate wind shear stress (and transforming wind speed from km/h to m/s)
bound$Uw <- 19.0 + 0.95 * (bound$vW * 1000/3600)^2 
bound$vW <- bound$vW * 1000/3600

boundary <- bound

# simulation maximum length
times <- seq(from = 1, to = max(boundary$Day), by = 1)
# initial water temperatures
yini <- c(3,3) 

calc_dens <- function(wtemp){
  dens = 999.842594 + (6.793952 * 10^-2 * wtemp) - (9.095290 * 10^-3 *wtemp^2) + (1.001685 * 10^-4 * wtemp^3) - (1.120083 * 10^-6* wtemp^4) + (6.536336 * 10^-9 * wtemp^5)
  return(dens)
}

parameters[19] = 3.2 # calibration parameter

out <- run_model(bc = boundary, params = parameters, ini = yini, times = times)

result <- data.frame('Time' = out[,1],
                     'WT_epi' = out[,2], 'WT_hyp' = out[,3])
g1 <- ggplot(result) +
  geom_line(aes(x=Time, y=WT_epi, col='Surface Mixed Layer')) +
  geom_line(aes(x=(Time), y=WT_hyp, col='Bottom Layer')) +
  labs(x = 'Simulated Time', y = 'WT in deg C')  +
  theme_bw()+
  guides(col=guide_legend(title="Layer")) +
  theme(legend.position="bottom");g1

output <- read.table('output.txt')

output <- data.frame('qin'=output[,1],'qout'=output[,2],'mix_e'=output[,3],'mix_h'=output[,4],
                     'sw'=output[,5],'lw'=output[,6],'water_lw'=output[,7],'conv'=output[,8],
                     'evap'=output[,9], 'Rh' = output[,10],
                     'entrain' = output[,11], 'Ri' = output[,12],'time' = output[,13])

output$balance <- apply(output[,c(1:9)],1, sum)

g2 <- ggplot(output) +
  geom_line(aes(x = time,y = qin, col = 'Inflow')) +
  geom_line(aes(x = time,y = qout, col = 'Outflow')) +
  geom_line(aes(x = time,y = mix_e, col = 'Mixing into Epilimnion')) +
  geom_line(aes(x = time,y = mix_h, col = 'Mixing into Hypolimnion')) +
  geom_line(aes(x = time,y = sw, col = 'Shortwave')) +
  geom_line(aes(x = time,y = lw, col = 'Longwave')) +
  geom_line(aes(x = time,y = water_lw, col = 'Reflection')) +
  geom_line(aes(x = time,y = conv, col = 'Conduction')) +
  geom_line(aes(x = time,y = evap, col = 'Evaporation')) +
  geom_line(aes(x = time,y = balance, col = 'Sum'), linetype = "dashed") +
  scale_colour_brewer("Energy terms", palette="Set3") +
  labs(x = 'Simulated Time', y = 'Fluxes in cal/(cm2 d)')  +
  theme_bw()+
  theme(legend.position="bottom");g2

g3 <- ggplot(output) +
  geom_line(aes(x = time,y = Ri, col = 'Richardson')) +
  geom_line(aes(x = time,y = entrain, col = 'Entrainment')) +
  scale_colour_brewer("Stability terms", palette="Set1") +
  labs(x = 'Simulated Time', y = 'Entrainment in cm2/s and Ri in [-]')  +
  theme_bw()+
  theme(legend.position="bottom");g3

g4 <- ggplot(boundary) +
  geom_line(aes(x = Day,y = scale(Jsw), col = 'Shortwave')) +
  geom_line(aes(x = Day,y = scale(Tair), col = 'Airtemp')) +
  geom_line(aes(x = Day,y = scale(Dew), col = 'Dew temperature')) +
  geom_line(aes(x = Day,y = scale(vW), col = 'Wind velocity')) +
  geom_line(aes(x = Day,y = scale(Uw), col = 'Wind shear stress')) +
  scale_colour_brewer("Boundary Conditions", palette="Set3") +
  labs(x = 'Simulated Time', y = 'Scaled meteo. boundaries')  +
  theme_bw()+
  theme(legend.position="bottom");g4

g5 <- grid.arrange(g1, g2, g3, g4, ncol =1);g5
ggsave(file='2L_visual_mendota.png', g5, dpi = 300,width = 200,height = 220, units = 'mm')


obs <-read_delim(paste0('obs.txt'), delim = ',')

obs_sfc <- obs %>%
  filter(Depth_meter == 1) %>%
  mutate('date' = yday(datetime),
         'sfc' = Water_Temperature_celsius)
obs_sfc$time = match(as.Date(obs_sfc$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))
obs_btm <- obs %>%
  filter(Depth_meter == 20) %>%
  mutate('date' = yday(datetime),
         'btm' = Water_Temperature_celsius) 
obs_btm$time = match(as.Date(obs_btm$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))


# all_obs <- data.frame('date' = yday(obs_sfc$datetime), 'sfc' = obs_sfc$Water_Temperature_celsius,
#                       'btm' = obs_btm$Water_Temperature_celsius)
# all_obs$date[nrow(all_obs)] = 366

g1 <- ggplot(result) +
  geom_line(aes(x=Time, y=WT_epi, col='Surface Mixed Layer (model)'), col = 'red') +
  geom_line(aes(x=(Time), y=WT_hyp, col='Bottom Layer (model)'), col = 'blue') +
  geom_point(data = obs_sfc, aes(x=time, y=sfc, col='Surface Mixed Layer (obs)'), col = 'red',linetype = "dashed") +
  geom_point(data = obs_btm, aes(x=(time), y=btm, col='Bottom Layer (obs)'), col = 'blue',linetype = "dashed") +
  labs(x = 'Simulated Time', y = 'WT in deg C')  +
  theme_bw()+
  guides(col=guide_legend(title="Layer")) +
  theme(legend.position="bottom");g1
ggsave(file='2L_compare_mendota.png', g1, dpi = 300,width = 300,height = 120, units = 'mm')
