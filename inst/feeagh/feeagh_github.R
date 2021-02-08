rm(list = ls())

library(gridExtra)
library(ggplot2)
library(pracma)
library(tidyverse)
library(RColorBrewer)
library(zoo)
library(deSolve)
library(LakeMetabolizer)

# load in the boundary data
bound <-read_delim(paste0('/home/robert/Projects/DSI/thermod/feeagh/meteo.txt'), delim = '\t')

urlfile = "https://raw.githubusercontent.com/robertladwig/robertladwig.github.io/master/files/meteo.txt"
bound <- read_tsv(url(urlfile))

# approximate the thermocline depth
therm_depth <- 10 ^ (0.336 * log10(3678) - 0.245) 
simple_therm_depth <- round(therm_depth)

Ve <- 2.886548e+13 # epilimnion volume (cm3)
Vh <- 3.421416e+13 # hypolimnion volume (cm3)
Ht <- 3 * 100 # thermocline thickness (cm)
At <- 26850835487 # thermocline area (cm2)
As <- 3.931e+10 # surface area (cm2)
Tin <- 0 # inflow water temperature
Q <- 0 # inflow discharge 
Rl <- 0.03 # reflection coefficient (generally small, 0.03)
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <- 11.7 * 10^(-8) # cal / (cm2 d K4) or: 4.9 * 10^(-3) # Stefan-Boltzmann constant in (J (m2 d K4)-1)
eps <- 0.97 # emissivity of water
rho <- 0.9982 # density (g per cm3)
cp <- 0.99 # specific heat (cal per gram per deg C)
c1 <- 0.47 # Bowen's coefficient
a <- 7 # constant
c <- 9e4 # empirical constant
g <- 9.81  # gravity (m/s2)
simple_therm_depth = simple_therm_depth

parameters <- c(Ve, Vh, At, Ht, As, Tin, Q, Rl, Acoeff, sigma, eps, rho, cp, c1, a, c, g, simple_therm_depth)

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


run_model <- function(bc, params, ini, times){
  
  Ve <- params[1]
  Vh <- params[2]
  At <- params[3]
  Ht <- params[4]
  As <- params[5]
  Tin <- params[6]
  Q <- params[7]
  Rl <- params[8]
  Acoeff <- params[9]
  sigma <- params[10]
  eps <- params[11]
  rho <- params[12]
  cp <- params[13]
  c1 <- params[14]
  a <- params[15]
  c <- params[16]
  g <- params[17]
  thermDep <- params[18]
  
  TwoLayer <- function(t, y, parms){
    eair <- (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t)))) # air vapor pressure
    esat <- 4.596 * exp((17.27 * Tair(t)) / (237.3 + Tair(t))) # saturation vapor pressure
    RH <- eair/esat *100 # relative humidity
    es <- 4.596 * exp((17.27 * y[1])/ (273.3+y[1]))
    # diffusion coefficient
    Cd <- 0.00052 * (vW(t))^(0.44)
    shear <- 1.164/1000 * Cd * (vW(t))^2
    rho_e <- calc_dens(y[1])/1000
    rho_h <- calc_dens(y[2])/1000
    w0 <- sqrt(shear/rho_e) 
    E0  <- c * w0
    Ri <- ((g/rho)*(abs(rho_e-rho_h)/10))/(w0/(thermDep)^2)
    if (rho_e > rho_h){
      dV = 100
    } else {
      dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000)
    }
    
    # epilimnion water temperature change per time unit
    dTe <-  Q / Ve * Tin -              # inflow heat
      Q / Ve * y[1] +                   # outflow heat
      ((dV * At) / Ve) * (y[2] - y[1]) + # mixing between epilimnion and hypolimnion
      + As/(Ve * rho * cp) * (
        Jsw(t)  + # shortwave radiation
          (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
          (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
          (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
          (Uw(t) * ((es) - (eair))) )# evaporation
    
    # hypolimnion water temperature change per time unit
    dTh <-  ((dV * At) / Vh) * (y[1] - y[2]) 

    # diagnostic variables for plotting
    mix_e <- ((dV * At) / Ve) * (y[2] - y[1])
    mix_h <-  ((dV * At) / Vh) * (y[1] - y[2]) 
    qin <- Q / Ve * Tin 
    qout <- - Q / Ve * y[1] 
    sw <- Jsw(t) 
    lw <- (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl))
    water_lw <- - (eps * sigma * (y[1]+ 273)^4)
    conv <- - (c1 * Uw(t) * (y[1] - Tair(t))) 
    evap <- - (Uw(t) * ((esat) - (eair)))
    Rh <- RH
    E <- (E0 / (1 + a * Ri)^(3/2))
    
    write.table(matrix(c(qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap, Rh,E, Ri, t), nrow=1), 'output.txt', append = TRUE,
                quote = FALSE, row.names = FALSE, col.names = FALSE)

    return(list(c(dTe, dTh)))
  }
  
  # approximating all boundary conditions  (linear interpolation)
  Jsw <- approxfun(x = bc$Day, y = bc$Jsw, method = "linear", rule = 2)
  Tair <- approxfun(x = bc$Day, y = bc$Tair, method = "linear", rule = 2)
  Dew <- approxfun(x = bc$Day, y = bc$Dew, method = "linear", rule = 2)
  Uw <- approxfun(x = bc$Day, y = bc$Uw, method = "linear", rule = 2)
  vW <- approxfun(x = bc$Day, y = bc$vW, method = "linear", rule = 2)
  
  # runge-kutta 4th order
  out <- ode(times = times, y = ini, func = TwoLayer, parms = params, method = 'rk4')
  
  return(out)
}


out <- run_model(bc = boundary, params = parameters, ini = yini, times = times)
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
  geom_line(aes(x = Day,y = scale(Jsw), col = 'Shortwave')) +
  geom_line(aes(x = Day,y = scale(Tair), col = 'Airtemp')) +
  geom_line(aes(x = Day,y = scale(Dew), col = 'Dew temperature')) +
  geom_line(aes(x = Day,y = scale(vW), col = 'Wind velocity')) +
  geom_line(aes(x = Day,y = scale(Uw), col = 'Wind shear stress')) +
  scale_colour_brewer("Boundary Conditions", palette="Set3") +
  labs(x = 'Simulated Time', y = 'Scaled meteo. boundaries')  +
  theme_bw()+
  theme(legend.position="bottom")

g5 <- grid.arrange(g1, g2, g3, g4, ncol =1);g5
ggsave(file='2L_visual_feeagh.png', g5, dpi = 300,width = 200,height = 220, units = 'mm')


obs <-read_delim(paste0('/home/robert/Projects/DSI/thermod/feeagh/obs.txt'), delim = ',')
obs_sfc <- obs %>%
  filter(Depth_meter == 0.9)
obs_btm <- obs %>%
  filter(Depth_meter == 42)
all_obs <- data.frame('date' = yday(obs_sfc$datetime), 'sfc' = obs_sfc$Water_Temperature_celsius,
                      'btm' = obs_btm$Water_Temperature_celsius)
all_obs$date[nrow(all_obs)] = 366

g1 <- ggplot(result) +
  geom_line(aes(x=Time, y=WT_epi, col='Surface Mixed Layer (model)'), col = 'red') +
  geom_line(aes(x=(Time), y=WT_hyp, col='Bottom Layer (model)'), col = 'blue') +
  geom_line(data = all_obs, aes(x=date, y=sfc, col='Surface Mixed Layer (obs)'), col = 'red',linetype = "dashed") +
  geom_line(data = all_obs, aes(x=(date), y=btm, col='Bottom Layer (obs)'), col = 'blue',linetype = "dashed") +
  labs(x = 'Simulated Time', y = 'WT in deg C')  +
  theme_bw()+
  guides(col=guide_legend(title="Layer")) +
  theme(legend.position="bottom");g1
ggsave(file='2L_compare.png', g1, dpi = 300,width = 200,height = 120, units = 'mm')
