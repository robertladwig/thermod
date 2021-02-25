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

### GETTING CONFIGURATION INPUT FROM LER YAML FILE
config_file <- 'LakeEnsemblR.yaml'
folder = '.'
parameters <- configure_from_ler(config_file <- config_file, folder = folder)

# load in the boundary data
bound <-read_delim(paste0(folder,'/meteo.txt'), delim = '\t')

colnames(bound) <- c('Day','Jsw','Tair','Dew','vW')

# function to calculate wind shear stress (and transforming wind speed from km/h to m/s)
bound$Uw <- 19.0 + 0.95 * (bound$vW * 1000/3600)^2 
bound$vW <- bound$vW * 1000/3600

boundary <- bound

# simulation maximum length
times <- seq(from = 1, to = max(boundary$Day), by = 1)
# initial water temperatures
yini <- c(3,3) 

if (file.exists('output.txt')){
  file.remove('output.txt')
}

parameters[19] = 3.2 # calibration parameter
ice_on = TRUE # ice "simulation" on or off?
out <- run_model(bc = boundary, params = parameters, ini = yini, times = times, ice = ice_on)

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
                     'entrain' = output[,11], 'Ri' = output[,12],'time' = output[,13],
                     'ice' = output[,14])

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

if (ice_on){
  ice_df = data.frame('time' = output$time)
  ice_df$ice = ifelse(output$ice == 1, "off", "on")
  g5 <- ggplot(ice_df) +
    geom_point(aes(x = time,y = ice)) +
    scale_colour_brewer("Boundary Conditions", palette="Set3") +
    labs(x = 'Simulated Time', y = 'Ice on/off')  +
    theme_bw()+
    theme(legend.position="bottom");g5
  
  g6 <- grid.arrange(g1, g2, g3, g4, g5, ncol =1);g6
  ggsave(file='2L_visual_mendota_ice.png', g6, dpi = 300,width = 200,height = 250, units = 'mm')
} else {
  g5 <- grid.arrange(g1, g2, g3, g4, ncol =1);g5
  ggsave(file='2L_visual_mendota.png', g5, dpi = 300,width = 200,height = 220, units = 'mm')
}




obs <-read_delim(paste0('obs.txt'), delim = ',')
simple_therm_depth = parameters[18]
start_date <- get_yaml_value(config_file, "time", "start")
stop_date <- get_yaml_value(config_file, "time", "stop")

### EITHER COMPARE AGAINST SPECIFIC DEPTHS
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

### OR COMPARE AGAINST AVERAGE OBSERVED DATA 
obs_sfc <- obs %>%
  filter(Depth_meter <= simple_therm_depth) %>%
  mutate('date' = yday(datetime),
         'sfc' = Water_Temperature_celsius) %>%
  group_by(datetime) %>%
  summarise('wtr_avg' = mean(sfc, na.rm = TRUE))
obs_sfc$time = match(as.Date(obs_sfc$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))
obs_btm <- obs %>%
  filter(Depth_meter >= simple_therm_depth) %>%
  mutate('date' = yday(datetime),
         'btm' = Water_Temperature_celsius) %>%
  group_by(datetime) %>%
  summarise('wtr_avg' = mean(btm, na.rm = TRUE))
obs_btm$time = match(as.Date(obs_btm$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))

g1 <- ggplot(result) +
  geom_line(aes(x=Time, y=WT_epi, col='Surface Mixed Layer (model)'), col = 'red') +
  geom_line(aes(x=(Time), y=WT_hyp, col='Bottom Layer (model)'), col = 'blue') +
  geom_point(data = obs_sfc, aes(x=time, y=wtr_avg, col='Surface Mixed Layer (obs)'), col = 'red',linetype = "dashed") + # sfc
  geom_point(data = obs_btm, aes(x=(time), y=wtr_avg, col='Bottom Layer (obs)'), col = 'blue',linetype = "dashed") + # btm
  labs(x = 'Simulated Time', y = 'WT in deg C')  +
  theme_bw()+
  guides(col=guide_legend(title="Layer")) +
  theme(legend.position="bottom");g1
ggsave(file='2L_compare_mendota.png', g1, dpi = 300,width = 300,height = 120, units = 'mm')


# Oxygen test simulation
# Fnep, Fsed, Ased, diffred 
wq_parameters <- append(parameters, c(0.001 / 1000, 
                                       100, 15000 * 1e4, 100))
wq_parameters[19] = parameters[19] # calibration parameter
# simulation maximum length
times <- seq(from = 1, to = max(boundary$Day), by = 1)
# initial water temperatures
yini <- c(3,3, 10 * 1000/1e6  * wq_parameters[1], 10 * 1000/1e6  * wq_parameters[2]) 

if (file.exists('output.txt')){
  file.remove('output.txt')
}


ice_on = TRUE # ice "simulation" on or off?
out <- run_oxygen_model(bc = boundary, params = wq_parameters, ini = yini, times = times, ice = ice_on)

result <- data.frame('Time' = out[,1],
                     'WT_epi' = out[,2], 'WT_hyp' = out[,3],
                     'DO_epi' = out[,4], 'DO_hyp' = out[,5])
head(result)

go1 <- ggplot(result) +
  geom_line(aes(x=Time, y=DO_epi / 1000 /  wq_parameters[1] * 1e6, col='Surface Mixed Layer')) +
  geom_line(aes(x=(Time), y=DO_hyp / 1000 /  wq_parameters[2] * 1e6, col='Bottom Layer')) +
  labs(x = 'Simulated Time', y = 'DO in g/m3')  +
  theme_bw()+
  guides(col=guide_legend(title="Layer")) +
  theme(legend.position="bottom");go1

obs <-read_delim(paste0('obs_oxygen.txt'), delim = ',')
simple_therm_depth = parameters[18]
start_date <- get_yaml_value(config_file, "time", "start")
stop_date <- get_yaml_value(config_file, "time", "stop")

### EITHER COMPARE AGAINST SPECIFIC DEPTHS
obs_sfc <- obs %>%
  filter(Depth_meter == 1) %>%
  mutate('date' = yday(datetime),
         'sfc' = Dissolved_Oxygen_gPerCubicMeter )
obs_sfc$time = match(as.Date(obs_sfc$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))
obs_btm <- obs %>%
  filter(Depth_meter == 20) %>%
  mutate('date' = yday(datetime),
         'btm' = Dissolved_Oxygen_gPerCubicMeter ) 
obs_btm$time = match(as.Date(obs_btm$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))

### OR COMPARE AGAINST AVERAGE OBSERVED DATA 
obs_sfc <- obs %>%
  filter(Depth_meter <= simple_therm_depth) %>%
  mutate('date' = yday(datetime),
         'sfc' = Dissolved_Oxygen_gPerCubicMeter ) %>%
  group_by(datetime) %>%
  summarise('do_avg' = mean(sfc, na.rm = TRUE))
obs_sfc$time = match(as.Date(obs_sfc$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))
obs_btm <- obs %>%
  filter(Depth_meter >= simple_therm_depth) %>%
  mutate('date' = yday(datetime),
         'btm' = Dissolved_Oxygen_gPerCubicMeter ) %>%
  group_by(datetime) %>%
  summarise('do_avg' = mean(btm, na.rm = TRUE))
obs_btm$time = match(as.Date(obs_btm$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))

go2 <- ggplot(result) +
  geom_line(aes(x=Time, y=DO_epi / 1000 /  wq_parameters[1] * 1e6, col='Surface Mixed Layer (model)')) +
  geom_line(aes(x=(Time), y=DO_hyp / 1000 /  wq_parameters[2] * 1e6, col='Bottom Layer (model)')) +
  geom_point(data = obs_sfc, aes(x=time, y=do_avg, col='Surface Mixed Layer (obs)'),linetype = "dashed") + # sfc
  geom_point(data = obs_btm, aes(x=(time), y=do_avg, col='Bottom Layer (obs)'),linetype = "dashed") + # btm
  labs(x = 'Simulated Time', y = 'DO in g/m3')  +
  scale_color_manual(values = c('blue','blue','red','red')) +
  theme_bw()+
  guides(col=guide_legend(title="Layer")) +
  theme(legend.position="bottom");go2
# ggsave(file='../../images/oxygen.png', go2, dpi = 300,width = 300,height = 120, units = 'mm')

go7 <- grid.arrange(go1, go2, ncol =1);go7
ggsave(file='2L_visual_mendota_oxygen.png', go7, dpi = 300,width = 200,height = 250, units = 'mm')
