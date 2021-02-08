rm(list = ls())

library(gridExtra)
library(ggplot2)
library(pracma)
library(tidyverse)
library(RColorBrewer)
library(zoo)

library(gotmtools)
library(LakeEnsemblR)
library(thermod)

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd('/home/robert/Projects/AEMONJ/LakeEnsemblR/example/feeagh')
if (file.exists('output.txt')) {
  # delete file if it exists
  file.remove('output.txt')
}
config_file <- 'Feeagh_master_config.yaml'
folder <- '.'
yaml  <-  file.path(folder, config_file)

meteo_file <- get_yaml_value(file = yaml, label = "meteo", key = "file")
met <- read.csv(file.path(folder, meteo_file), stringsAsFactors = F)
met[, 1] <- as.POSIXct(met[, 1])
# Check time step
tstep <- diff(as.numeric(met[, 1]))
lev <- get_yaml_value(config_file, "location", "elevation")
# Maximum Depth
max_depth <- get_yaml_value(config_file, "location", "depth")
# Read in hypsograph data
hyp_file <- get_yaml_value(config_file, "location", "hypsograph")
hyp <- read.csv(hyp_file)
# Start date
start_date <- get_yaml_value(config_file, "time", "start")
# Stop date
stop_date <- get_yaml_value(config_file, "time", "stop")
# Time step
timestep <- get_yaml_value(config_file, "time", "time_step")

bsn_len <- get_yaml_value(config_file, "model_parameters", "bsn_len") 

bsn_wid <- get_yaml_value(config_file, "model_parameters", "bsn_wid") 

therm_depth <- 10 ^ (0.336 * log10(max(bsn_len, bsn_wid)) - 0.245)

simple_therm_depth <- round(therm_depth)

ix <- match(simple_therm_depth, hyp$Depth_meter)
# input data: month, shortwave radiation, air temperature, dew point temperature, wind speed

met$Dewpoint_Air_Temperature_Celsius <- met$Air_Temperature_celsius - ((100 - met$Relative_Humidity_percent)/5.)
bound <- met %>%
  dplyr::select(datetime, Shortwave_Radiation_Downwelling_wattPerMeterSquared, 
         Air_Temperature_celsius, Dewpoint_Air_Temperature_Celsius, Ten_Meter_Elevation_Wind_Speed_meterPerSecond)
bound$Shortwave_Radiation_Downwelling_wattPerMeterSquared <- bound$Shortwave_Radiation_Downwelling_wattPerMeterSquared * 2.0E-5 *3600*24
bound <- bound %>%
  dplyr::filter(datetime >= start_date & datetime <= stop_date)
bound$datetime <- seq(1, nrow(bound))


write_delim(bound, path = paste0('/home/robert/Projects/DSI/thermod/feeagh/meteo.txt'), delim = '\t')


bound <-read_delim(paste0('/home/robert/Projects/DSI/thermod/feeagh/meteo.txt'), delim = '\t')


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

parameters <- c(Ve, Vh, At, Ht, As, Tin, Q, Rl, Acoeff, sigma, eps, rho, cp, c1, a, c, g)

colnames(bound) <- c('Month','Jsw','Tair','Dew','vW')
bound$Uw <- 19.0 + 0.95 * (bound$vW * 1000/3600)^2 # function to calculate wind shear stress (and transforming wind speed from km/h to m/s)
bound$vW <- bound$vW * 1000/3600

boundary <- bound
# boundary <- bound[cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30,31)), ]
# cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30,31))
# boundary <- rbind(bound[1,], apply(bound[1:32,], 2, mean),
#                        apply(bound[32:60,], 2, mean),
#                        apply(bound[60:91,], 2, mean),
#                        apply(bound[91:121,], 2, mean),
#                        apply(bound[121:152,], 2, mean),
#                        apply(bound[152:182,], 2, mean),
#                        apply(bound[182:213,], 2, mean),
#                        apply(bound[213:244,], 2, mean),
#                        apply(bound[244:274,], 2, mean),
#                        apply(bound[274:305,], 2, mean),
#                        apply(bound[305:335,], 2, mean),
#                        apply(bound[335:366,], 2, mean))
# boundary <- as.data.frame(boundary)
# boundary$Month <- cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30,31))

# boundary <- add_noise(boundary)

times <- seq(from = 1, to = max(boundary$Month), by = 1)
yini <- c(3,3) # initial water temperatures

out <- run_model(modelfunc = 'TwoLayer', bc = boundary, params = parameters, ini = yini, times = times)
# out <- run_model(modelfunc = model, bc = boundary, params = parameters, ini = c(yini, 8/1000*Ve, 8/1000*Vh), 
                 # times = times) # c(yini, 8/1000*Ve, 8/1000*Vh)

result <- data.frame('Time' = out[,1],
                     'WT_epi' = out[,2], 'WT_hyp' = out[,3])
g1 <- ggplot(result) +
  geom_line(aes(x=Time, y=WT_epi, col='Surface Mixed Layer')) +
  geom_line(aes(x=(Time), y=WT_hyp, col='Bottom Layer')) +
  labs(x = 'Simulated Time', y = 'WT in deg C')  +
  theme_bw()+
  guides(col=guide_legend(title="Layer")) +
  theme(legend.position="bottom")

output <- read.table('output.txt')

output <- data.frame('qin'=output[,1],'qout'=output[,2],'mix_e'=output[,3],'mix_h'=output[,4],
                     'sw'=output[,5],'lw'=output[,6],'water_lw'=output[,7],'conv'=output[,8],
                     'evap'=output[,9], 'Rh' = output[,10],
                     'entrain' = output[,11], 'Ri' = output[,12],'time' = output[,13])
#Fatm, Sed, PP, VOl, Oflux_epi, Oflux_hypo
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
  theme(legend.position="bottom")

g3 <- ggplot(output) +
  geom_line(aes(x = time,y = Ri, col = 'Richardson')) +
  geom_line(aes(x = time,y = entrain, col = 'Entrainment')) +
  scale_colour_brewer("Stability terms", palette="Set1") +
  labs(x = 'Simulated Time', y = 'Entrainment in cm2/s and Ri in [-]')  +
  theme_bw()+
  theme(legend.position="bottom")

g4 <- ggplot(boundary) +
  geom_line(aes(x = Month,y = scale(Jsw), col = 'Shortwave')) +
  geom_line(aes(x = Month,y = scale(Tair), col = 'Airtemp')) +
  geom_line(aes(x = Month,y = scale(Dew), col = 'Dew temperature')) +
  geom_line(aes(x = Month,y = scale(vW), col = 'Wind velocity')) +
  geom_line(aes(x = Month,y = scale(Uw), col = 'Wind shear stress')) +
  scale_colour_brewer("Boundary Conditions", palette="Set3") +
  labs(x = 'Simulated Time', y = 'Scaled meteo. boundaries')  +
  theme_bw()+
  theme(legend.position="bottom")


g5 <- grid.arrange(g1, g2, g3, g4, ncol =1);g5
ggsave(file='2L_visual_feeagh.png', g5, dpi = 300,width = 200,height = 220, units = 'mm')

sim_wtr <- matrix(NA, ncol = nrow(out), nrow= (length(hyp$Depth_meter)))
for (ii in 1:ncol(sim_wtr)){
  sim_wtr[1:ix, ii] <- out[ii, 2]
  sim_wtr[(ix + Ht/100):nrow(sim_wtr), ii] <- out[ii, 3]
  sim_wtr[, ii] <- na.approx(sim_wtr[, ii])
}

sim_wtr <- as.matrix( arrange(as.data.frame(sim_wtr), -row_number()))
heat_data <- list('time' = output$time, 
                        'depth' = hyp$Depth_meter)

library(RColorBrewer)
spectral.colors <-  colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral') )
inMeter <- function(x) {paste0(x, " m")}
inCelsius <- function(x) {paste0(x, paste(" ",unit,sep=''))}
tlevels <- seq(-5,20,1)
filled.contour(x=seq(1, nrow(out)), y=heat_data$depth, z = t(sim_wtr),
               levels=tlevels,
               col=rev(spectral.colors(length(tlevels))), cex = 1.5, cex.main = 3.)
               # plot.axes = {axis(1, labels=format(pretty(datetime,20), "%Y-%b"), at = as.numeric(pretty(datetime,20)), cex.axis=2,las=1.5);
               #   axis(2, labels=inMeter(rev(seq(0,max.plot.dep,iter.plot.dep))*(-1)), at = rev(seq(0,max(max.plot.dep),iter.plot.dep))*(-1), cex.axis=1.5)},
               # key.axes = {axis(4,at=unique(tlevels),labels=(inCelsius(unique(tlevels))))})
