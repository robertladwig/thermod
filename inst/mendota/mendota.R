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

bound <- bound[, -c(1)]
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
  geom_line(aes(x=Time, y=WT_epi, col='Surface Mixed Layer (model)')) +
  geom_line(aes(x=(Time), y=WT_hyp, col='Bottom Layer (model)')) +
  geom_point(data = obs_sfc, aes(x=time, y=wtr_avg, col='Surface Mixed Layer (obs)'),linetype = "dashed") + # sfc
  geom_point(data = obs_btm, aes(x=(time), y=wtr_avg, col='Bottom Layer (obs)'), linetype = "dashed") + # btm
  labs(x = 'Simulated Time', y = 'WT in deg C')  +
  scale_color_manual(values = c('blue','blue','red','red')) +
  theme_bw()+
  guides(col=guide_legend(title="Layer")) +
  theme(legend.position="bottom");g1
ggsave(file='2L_compare_mendota.png', g1, dpi = 300,width = 300,height = 120, units = 'mm')

result_filter = result
result_filter$WT_epi = kalman_filtering(time = result_filter$Time, series = result_filter$WT_epi)
result_filter$WT_hyp = kalman_filtering(time = result_filter$Time, series = result_filter$WT_hyp)
g1 <- ggplot(result_filter) +
  geom_line(aes(x=as.numeric(Time), y=WT_epi, col='Surface Mixed Layer (model)')) +
  geom_line(aes(x=as.numeric(Time), y=WT_hyp, col='Bottom Layer (model)')) +
  geom_point(data = obs_sfc, aes(x=as.numeric(time), y=wtr_avg, col='Surface Mixed Layer (obs)'),linetype = "dashed") + # sfc
  geom_point(data = obs_btm, aes(x=as.numeric(time), y=wtr_avg, col='Bottom Layer (obs)'), linetype = "dashed") + # btm
  labs(x = 'Simulated Time', y = 'WT in deg C')  +
  scale_color_manual(values = c('blue','blue','red','red')) +
  theme_bw()+
  guides(col=guide_legend(title="Layer")) +
  theme(legend.position="bottom");g1
ggsave(file='2L_compare_mendota_filter.png', g1, dpi = 300,width = 300,height = 120, units = 'mm')

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



# NPZ tests
# Fnep, Fsed, Ased, diffred 
#kg, kra, Cgz, ksa, aca, epsilon, krz ksz, apa, apc, Fsedp, alpha1, alpha2
test <- c(0.165, 0.15, 1.5*1e3, 1e-5, 0.04 * 1000, 0.6, 0.1 , 1e-5, 15, 0.15, -10000 / 1e4, 3e1, 3e1,1000) #0.35 -3000
npz_parameters <- append(parameters, c(0.001 / 1000, 
                                       100, 25000 * 1e4, 100,
                                       test))
                
npz_parameters[19] = parameters[19] # calibration parameter
npz_parameters[23] = wq_parameters[23] * 1e-1 * 2.5 #npz_parameters[23] / 5#5# calibration parameter
# simulation maximum length
times <- seq(from = 1, to = max(boundary$Day), by = 1)
# initial water temperatures
yini <- c(3,3, 10 * 1000/1e6  * wq_parameters[1], 10 * 1000/1e6  * wq_parameters[2],
          1/1e6  * wq_parameters[1], 1/1e6  * wq_parameters[2],
          0.05 * 1000/1e6  * wq_parameters[1], 0.05 * 1000/1e6  * wq_parameters[2],
          20 /1000/(0.001 * 10e6) * wq_parameters[1], 20 /1000/(0.001 * 10e6) * wq_parameters[2]) 

if (file.exists('output.txt')){
  file.remove('output.txt')
}


ice_on = TRUE # ice "simulation" on or off?
out <- run_npz_model(bc = boundary, params = npz_parameters, ini = yini, times = times, ice = ice_on)

output <- read.table('output.txt')
#qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap, Rh,E, Ri, t, ice_param,
# ATM, NEP, SED, Oflux_epi, Oflux_hypo,
# phy_growth,phy_grazing, zoo_grazing,
# zoo_decay, phy_to_nutr, zoo_to_nutr,
# nutr_to_phy)
output <- data.frame('qin'=output[,1],'qout'=output[,2],'mix_e'=output[,3],'mix_h'=output[,4],
                     'sw'=output[,5],'lw'=output[,6],'water_lw'=output[,7],'conv'=output[,8],
                     'evap'=output[,9], 'Rh' = output[,10],
                     'entrain' = output[,11], 'Ri' = output[,12],'time' = output[,13],
                     'ice' = output[,14], 'Fatm' = output[,15], 'Fnep' = output[,16],
                     'Fsed' = output[,17], 'entr_epi' = output[,18], 'entr_hypo' = output[,19],
                     'phy_growth' = output[,20], 'phy_grazing' = output[,21], 
                     'zoo_grazing' =output[,22], 'zoo_decay' = output[,23],
                     'phy_to_nutr' = output[,24], 'zoo_to_nutr' = output[,25],
                     'nutr_to_phy' = output[,26])

result <- data.frame('Time' = out[,1],
                     'WT_epi' = out[,2], 'WT_hyp' = out[,3],
                     'DO_epi' = out[,4] / 1000 /  wq_parameters[1] * 1e6, 
                     'DO_hyp' = out[,5] / 1000 /  wq_parameters[1] * 1e6,
                     'PHY_epi' = out[,6] / 1000 /  wq_parameters[1] * 1e6, 
                     'PHY_hyp' = out[,7] / 1000 /  wq_parameters[1] * 1e6,
                     'ZOO_epi' = out[,8] / 1000 /  wq_parameters[1] * 1e6, 
                     'ZOO_hyp' = out[,9] / 1000 /  wq_parameters[1] * 1e6,
                     'P_epi' = out[,10] / 1000 /  wq_parameters[1] * 1e6, 
                     'P_hyp' = out[,11] / 1000 /  wq_parameters[1] * 1e6)
head(result)

m.result <- reshape2::melt(result, id = 'Time')
g <- ggplot(m.result) +
  geom_line(aes(Time, value)) + # / 1000 /  wq_parameters[1] * 1e6
  facet_wrap(~ variable, scales = 'free') +
  theme_minimal(); g
ggsave(file='NPZ.png', g, dpi = 300,width = 250,height = 200, units = 'mm')

time_seq <- seq.Date(from =   as.Date(get_yaml_value(config_file, "time", "start")),
                     to = as.Date(get_yaml_value(config_file, "time", "stop")),
                     by = 'day')
result$Date = time_seq

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
  geom_line(aes(Date, WT_epi), col = 'red') +
  labs(x = '', y = 'WTR in deg C')  +
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  geom_point(data = obs_sfc, aes(x=as.Date(datetime), y=wtr_avg), col = 'orange', alpha=0.3) +
  theme_bw()+
  theme(legend.position="bottom");g1
g2 <- ggplot(result) +
  geom_line(aes(Date, WT_hyp), col = 'red') +
  labs(x = '', y = 'WTR in deg C')  +
  theme_bw()+
  geom_point(data = obs_btm, aes(x=as.Date(datetime), y=wtr_avg), col = 'orange', alpha=0.3) +
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  theme(legend.position="bottom");g2

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

g3 <- ggplot(result) +
  geom_line(aes(Date, DO_epi), col = 'blue') +
  labs(x = '', y = 'DO in g/m3')  +
  geom_point(data = obs_sfc, aes(x=as.Date(datetime), y=do_avg), col = 'cyan', alpha=0.3) + # sfc
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  theme_bw()+
  theme(legend.position="bottom");g3
g4 <- ggplot(result) +
  geom_line(aes(Date, DO_hyp), col = 'blue') +
  geom_point(data = obs_btm, aes(x=as.Date(datetime), y=do_avg), col = 'cyan', alpha=0.3) + # btm
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  labs(x = '', y = 'DO in g/m3')  +
  theme_bw()+
  theme(legend.position="bottom");g4

# microgramper liter
obs_chla <- read.csv('obs_chla.csv')
obs_chla2 = obs_chla %>%
  filter(depth_range_m %in% c('0', '4', '8')) %>%
  group_by(sampledate) %>%
  summarise(mean(correct_chl_fluor)) %>%
  rename(correct_chl_fluor = `mean(correct_chl_fluor)`)
obs_chla1 = obs_chla %>%
  filter(depth_range_m %in% c('0-8')) %>%
  select(sampledate, correct_chl_fluor)
obs_chla_cleaned = merge(obs_chla2, obs_chla1, by = 'sampledate')

obs_sfc <- obs_chla1
obs_sfc$time = match(as.Date(obs_chla1$sampledate), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))


g5 <- ggplot(result) +
  geom_line(aes(Date, PHY_epi), col = 'darkgreen') +
  labs(x = '', y = 'Chla in g/m3')  +
  geom_point(data = obs_sfc, aes(x=as.Date(sampledate), y=correct_chl_fluor/1000), col = 'lightgreen', alpha=0.3) + # sfc
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  theme_bw()+
  theme(legend.position="bottom");g5
g6 <- ggplot(result) +
  geom_line(aes(Date, PHY_hyp), col = 'darkgreen') +
  labs(x = '', y = 'Chla in g/m3')  +
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  theme_bw()+
  theme(legend.position="bottom");g6

g7 <- ggplot(result) +
  geom_line(aes(Date, ZOO_epi), col = 'brown') +
  labs(x = '', y = 'C in g/m3')  +
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  theme_bw()+
  theme(legend.position="bottom");g7
g8 <- ggplot(result) +
  geom_line(aes(Date, ZOO_hyp), col = 'brown') +
  labs(x = '', y = 'C in g/m3')  +
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  theme_bw()+
  theme(legend.position="bottom");g8

obs_nutrient <- read.csv('obs_nutrient.csv')
obs_nutrient$datetime <- as.POSIXct(obs_nutrient$datetime)

obs_sfc <- obs_nutrient %>%
  filter(depth <= simple_therm_depth) %>%
  mutate('date' = yday(datetime),
         'sfc' = PHS_frp ) %>%
  group_by(datetime) %>%
  summarise('p_avg' = mean(sfc, na.rm = TRUE))
obs_sfc$time = match(as.Date(obs_sfc$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))
obs_btm <- obs_nutrient %>%
  filter(depth >= simple_therm_depth) %>%
  mutate('date' = yday(datetime),
         'btm' = PHS_frp ) %>%
  group_by(datetime) %>%
  summarise('p_avg' = mean(btm, na.rm = TRUE))
obs_btm$time = match(as.Date(obs_btm$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))

g9 <- ggplot(result) +
  geom_line(aes(Date, P_epi), col = 'orange') +
  labs(x = '', y = 'P in g/m3')  +
  geom_point(data = obs_sfc, aes(x=as.Date(datetime), y=p_avg), col = 'yellow', alpha=0.3) + # sfc
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  theme_bw()+
  theme(legend.position="bottom");g9
g10 <- ggplot(result) +
  geom_line(aes(Date, P_hyp), col = 'orange') +
  labs(x = '', y = 'P in g/m3')  +
  geom_point(data = obs_btm, aes(x=as.Date(datetime), y=p_avg), col = 'yellow', alpha=0.3) +
  scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
  theme_bw()+
  theme(legend.position="bottom");g10

library(patchwork)
g <- (g1 / g3 / g5 / g7 / g9) | (g2 / g4 / g6 / g8 / g10)  +
  plot_annotation(
    title = 'Thermod, 2-layer NPZ model, v1.01',
    subtitle = 'Epilimnion (left), Hypolimnion (right)',
    caption = 'Under development'
  ); g
ggsave(file='2L_visual_mendota_nutrientfoodweb_random.png', g, dpi = 300,width = 200,height = 170, units = 'mm')

# animation loop
for (dim in seq(from= 2, to =(nrow(result)), by = 7)){

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
  
  g1 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, WT_epi), col = 'red') +
    labs(x = '', y = 'WTR in deg C')  +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    geom_point(data = obs_sfc, aes(x=as.Date(datetime), y=wtr_avg), col = 'orange', alpha=0.3) +
    ylim(min(result$WT_epi, obs_sfc$wtr_avg), max(result$WT_epi, obs_sfc$wtr_avg)) +
    theme_bw()+
    theme(legend.position="bottom");g1
  g2 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, WT_hyp), col = 'red') +
    labs(x = '', y = 'WTR in deg C')  +
    theme_bw()+
    geom_point(data = obs_btm, aes(x=as.Date(datetime), y=wtr_avg), col = 'orange', alpha=0.3) +
    ylim(min(result$WT_hyp, obs_btm$wtr_avg), max(result$WT_hyp, obs_btm$wtr_avg)) +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    theme(legend.position="bottom");g2
  
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
  
  g3 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, DO_epi), col = 'blue') +
    labs(x = '', y = 'DO in g/m3')  +
    geom_point(data = obs_sfc, aes(x=as.Date(datetime), y=do_avg), col = 'cyan', alpha=0.3) + # sfc
    ylim(min(result$DO_epi, obs_sfc$do_avg), max(result$DO_epi, obs_sfc$do_avg)) +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    theme_bw()+
    theme(legend.position="bottom");g3
  g4 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, DO_hyp), col = 'blue') +
    geom_point(data = obs_btm, aes(x=as.Date(datetime), y=do_avg), col = 'cyan', alpha=0.3) + # btm
    ylim(min(result$DO_hyp, obs_btm$do_avg), max(result$DO_hyp, obs_btm$do_avg)) +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    labs(x = '', y = 'DO in g/m3')  +
    theme_bw()+
    theme(legend.position="bottom");g4
  
  # microgramper liter
  obs_chla <- read.csv('obs_chla.csv')
  obs_chla2 = obs_chla %>%
    filter(depth_range_m %in% c('0', '4', '8')) %>%
    group_by(sampledate) %>%
    summarise(mean(correct_chl_fluor)) %>%
    rename(correct_chl_fluor = `mean(correct_chl_fluor)`)
  obs_chla1 = obs_chla %>%
    filter(depth_range_m %in% c('0-8')) %>%
    select(sampledate, correct_chl_fluor)
  obs_chla_cleaned = merge(obs_chla2, obs_chla1, by = 'sampledate')
  
  obs_sfc <- obs_chla1
  obs_sfc$time = match(as.Date(obs_chla1$sampledate), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))
  
  
  g5 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, PHY_epi), col = 'darkgreen') +
    labs(x = '', y = 'Chla in g/m3')  +
    geom_point(data = obs_sfc, aes(x=as.Date(sampledate), y=correct_chl_fluor/1000), col = 'lightgreen', alpha=0.3) + # sfc
    ylim(min(result$PHY_epi, obs_sfc$correct_chl_fluor/1000), max(result$PHY_epi, obs_sfc$correct_chl_fluor/1000)) +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    theme_bw()+
    theme(legend.position="bottom");g5
  g6 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, PHY_hyp), col = 'darkgreen') +
    labs(x = '', y = 'Chla in g/m3')  +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    ylim(min(result$PHY_hyp), max(result$PHY_hyp)) +
    theme_bw()+
    theme(legend.position="bottom");g6
  
  g7 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, ZOO_epi), col = 'brown') +
    labs(x = '', y = 'C in g/m3')  +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    ylim(min(result$ZOO_epi), max(result$ZOO_epi)) +
    theme_bw()+
    theme(legend.position="bottom");g7
  g8 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, ZOO_hyp), col = 'brown') +
    labs(x = '', y = 'C in g/m3')  +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    ylim(min(result$ZOO_hyp), max(result$ZOO_hyp)) +
    theme_bw()+
    theme(legend.position="bottom");g8
  
  obs_nutrient <- read.csv('obs_nutrient.csv')
  obs_nutrient$datetime <- as.POSIXct(obs_nutrient$datetime)
  
  obs_sfc <- obs_nutrient %>%
    filter(depth <= simple_therm_depth) %>%
    mutate('date' = yday(datetime),
           'sfc' = PHS_frp ) %>%
    group_by(datetime) %>%
    summarise('p_avg' = mean(sfc, na.rm = TRUE))
  obs_sfc$time = match(as.Date(obs_sfc$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))
  obs_btm <- obs_nutrient %>%
    filter(depth >= simple_therm_depth) %>%
    mutate('date' = yday(datetime),
           'btm' = PHS_frp ) %>%
    group_by(datetime) %>%
    summarise('p_avg' = mean(btm, na.rm = TRUE))
  obs_btm$time = match(as.Date(obs_btm$datetime), seq(as.Date(start_date), as.Date(stop_date), by = 'day'))
  
  g9 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, P_epi), col = 'orange') +
    labs(x = '', y = 'P in g/m3')  +
    geom_point(data = obs_sfc, aes(x=as.Date(datetime), y=p_avg), col = 'yellow', alpha=0.3) + # sfc
    ylim(min(result$P_epi, obs_sfc$p_avg), max(result$P_epi, obs_sfc$p_avg)) +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    theme_bw()+
    theme(legend.position="bottom");g9
  g10 <- ggplot(result[1:dim,]) +
    geom_line(aes(Date, P_hyp), col = 'orange') +
    labs(x = '', y = 'P in g/m3')  +
    geom_point(data = obs_btm, aes(x=as.Date(datetime), y=p_avg), col = 'yellow', alpha=0.3) +
    ylim(min(result$P_hyp, obs_btm$p_avg), max(result$P_hyp, obs_btm$p_avg)) +
    scale_x_date(limits= as.Date(c(min(result$Date), max(result$Date))))+
    theme_bw()+
    theme(legend.position="bottom");g10
  
  library(patchwork)
  g <- (g1 / g3 / g5 / g7 / g9) | (g2 / g4 / g6 / g8 / g10)  +
    plot_annotation(
      title = 'Thermod, 2-layer NPZ model, v1.01',
      subtitle = 'Epilimnion (left), Hypolimnion (right)',
      caption = 'Under development'
    ); g
  # ggsave(file='2L_visual_mendota_nutrientfoodweb_random.png', g, dpi = 300,width = 350,height = 250, units = 'mm')
  # ggsave(file='2L_visual_mendota_nutrientfoodweb_random.png', g, dpi = 300,width = 200,height = 170, units = 'mm')
  
ggsave(file=paste0('anim/',dim,'.png'), g, dpi = 300,width = 200,height = 170, units = 'mm')

}
# 
# devtools::install_github("https://github.com/inrae/RCaN.git", subdir="RCaN")
# library(RCaN)
# 
# # install.packages('rEDM')
# library(rEDM)
# str(block_3sp)
# Mview = Multiview(dataFrame = result,
#                   lib = '1 100',
#                   pred = '101 190',
#                   E = 3,
#                   columns = 'WT_epi ZOO_epi P_epi',
#                   target = 'PHY_epi')
