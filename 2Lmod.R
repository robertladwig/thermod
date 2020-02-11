rm(list = ls())
library(deSolve)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(thermod)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

if (file.exists('output.txt')) 
  #Delete file if it exists
  file.remove('output.txt')

heat.fluxes <-c()

calc_dens <- function(wtemp){
  dens = 999.842594 + (6.793952 * 10^-2 * wtemp) - (9.095290 * 10^-3 *wtemp^2) + (1.001685 * 10^-4 * wtemp^3) - (1.120083 * 10^-6* wtemp^4) + (6.536336 * 10^-9 * wtemp^5)
  return(dens)
}

TwoLayer <- function(t, y, parms){
  eair <- (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t)))) # air vapor pressure
  esat <- 4.596 * exp((17.27 * Tair(t)) / (237.3 + Tair(t))) # saturation vapor pressure
  RH <- eair/esat *100 # relative humidity
  es <- 4.596 * exp((17.27 * y[1])/ (273.3+y[1]))
  # diffusion coefficient
  Cd <- 0.00052 * (vW(t))^(0.44)
  shear <- 1.164/1000 * Cd * (vW(t))^2
#  shear <- Uw(t)
  rho_e <- calc_dens(y[1])/1000
  rho_h <- calc_dens(y[2])/1000
  w0 <- sqrt(shear/rho_e) # sqrt(shear/rho)
  E0  <- c * w0
  Ri <- ((g/rho)*(abs(rho_e-rho_h)/10))/(w0/(10)^2)
  if (rho_e > rho_h){
    dV = 100
  } else {
    dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000)
  }
  # dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000)
  
  # epilimnion water temperature change per time unit
  dTe <-  Q / Ve * Tin -              # inflow heat
    Q / Ve * y[1] +                   # outflow heat
    ((dV * At) / Ve) * (y[2] - y[1]) + # mixing between epilimnion and hypolimnion # ((vt(t) * At) / Ve) * (y[2] - y[1])
    + As/(Ve * rho * cp) * (
    Jsw(t)  + # shortwave radiation
    (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
    (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
    (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
    (Uw(t) * ((es) - (eair))) )# evaporation
  
  # hypolimnion water temperature change per time unit
  dTh <-  ((dV * At) / Vh) * (y[1] - y[2]) # ((vt(t) * At) / Vh) * (y[1] - y[2]) 


  
  
  # diagnostic variables for plotting
  mix_e <- (vt(t) * At * rho * cp) * (y[2] - y[1])  /As#((vt(t) * At) / Ve) * (y[2] - y[1])
  mix_h <- (vt(t) * At *rho * cp) * (y[1] - y[2]) /At#(vt(t) * At) / (Vh) * (y[1] - y[2])
  qin <- Q * rho * cp * Tin /As#Q / Ve * Tin
  qout <- - Q * rho * cp * y[1] /As #- Q / Ve * y[1] 
  sw <- Jsw(t) #*As#* As/(Ve * rho * cp)
  lw <- (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) #* As#* As/(Ve * rho * cp)
  water_lw <- - (eps * sigma * (y[1]+ 273)^4) #*As#*As/(Ve * rho * cp)
  conv <- - (c1 * Uw(t) * (y[1] - Tair(t))) #*As#*As/(Ve * rho * cp)
  evap <- - (Uw(t) * ((esat) - (eair))) #*As #*As/(Ve * rho * cp)
  Rh <- RH
  E <- (E0 / (1 + a * Ri)^(3/2))
  
  write.table(matrix(c(qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap, Rh,E, Ri, t), nrow=1), 'output.txt', append = TRUE,
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  return(list(c(dTe, dTh)))
}

# input data: month, shortwave radiation, air temperature, dew point temperature, wind speed
bound <- matrix(c(seq(1,12,1),
                  169, 274, 414, 552, 651, 684, 642, 537, 397, 259, 160, 127,
                  8.3, 9., 13.5,13.9,21.8,24.7,29.4,26.6,24.9,15.,9.7,6.6,
                  2.8,3.3,4.9,4.,5.3,7.8,11.8,11.5,7.7,6.8,6.5,2.4,
                  11.6,11.7,16.4,15.6,16.6,16.7,12.7,11.7,14.,12.9,14.8,11.6), nrow = 12, byrow = FALSE)
Ve <- 175000 *1e6 # epilimnion volume
Vh <- 75000*1e6 # hypolimnion volume
At <- 11000 *1e4 # thermocline area
Ht <- 3 * 100 # thermocline thickness
As <- 25000 * 1e4 # surface area
# H <- As/Ve#Ve/As *100
Tin <- 10 # inflow water temperature
Q <- 7500 *1e6 # inflow discharge
Rl <- 0.03 # reflection coefficient (generally small, 0.03)
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <- 11.7 * 10^(-8) # cal / (cm2 d K4) or: 4.9 * 10^(-3) # Stefan-Boltzmann constant in [J (m2 d K4)-1]
eps <- 0.97 # emissivity of water
rho <- 0.9982 # density (g per cm3)
cp <- 0.99 # specific heat (cal per gram per deg C)
c1 <- 0.47 # Bowen's coefficient
a <- 7#10e-1
c <- 9e4#25#50
g <- 9.81

parameters <- c(Ve, Vh, At, Ht, As, Tin, Q, Rl, Acoeff, sigma, eps, rho, cp, c1, a, c, g)

Et <- 7.07 * 10^(-4)  * ((Ve+Vh)/As/100)^(1.1505) # vertifcal diffusion coefficient (cm2 per d)
vto <- Et/(Ht/100) * (86400/10000) #*100 # heat change coefficient across thermocline during stratified season
vt_mix <- 100 # high heat change coefficient during overturn and mixing events
bound <- as.data.frame(cbind(bound, c(rep(vt_mix,3),rep(vto,6),rep(vt_mix,3))))

colnames(bound) <- c('Month','Jsw','Tair','Dew','vW','vt')
bound$Uw <- 19.0 + 0.95 * (bound$vW * 1000/3600)^2 # function to calculate wind shear stress (and transforming wind speed from km/h to m/s)
bound$vW <- bound$vW * 1000/3600

# boundary <- bound
# boundary$Month <- seq(1, nrow(boundary),1)
boundary <- rbind(bound[1,],bound)

boundary <- rbind(boundary, boundary, boundary, boundary, boundary,
                  boundary, boundary, boundary, boundary, boundary)

boundary$Month <- cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30,31, 31,31,28,31,30,31,30,31,31,30,31,30,31,
                           31,31,28,31,30,31,30,31,31,30,31,30,31,
                           31,31,28,31,30,31,30,31,31,30,31,30,31,
                           31,31,28,31,30,31,30,31,31,30,31,30,31,
                           31,31,28,31,30,31,30,31,31,30,31,30,31,
                           31,31,28,31,30,31,30,31,31,30,31,30,31,
                           31,31,28,31,30,31,30,31,31,30,31,30,31,
                           31,31,28,31,30,31,30,31,31,30,31,30,31,
                           31,31,28,31,30,31,30,31,31,30,31,30,31))



# approximating all boundary conditions 
Jsw <- approxfun(x = boundary$Month, y = boundary$Jsw, method = "linear", rule = 2)
Tair <- approxfun(x = boundary$Month, y = boundary$Tair, method = "linear", rule = 2)
Dew <- approxfun(x = boundary$Month, y = boundary$Dew, method = "linear", rule = 2)
Uw <- approxfun(x = boundary$Month, y = boundary$Uw, method = "linear", rule = 2)
vW <- approxfun(x = boundary$Month, y = boundary$vW, method = "linear", rule = 2)
vt <- approxfun(x = boundary$Month, y = boundary$vt, method = "linear", rule = 2)

#times <- seq(from = 1, to = nrow(boundary), by = 1/30)
times <- seq(from = 1, to = max(boundary$Month), by = 1)
yini <- c(5,5) # initial water temperatures

# runge-kutta 4th order
out <- ode(times = times, y = yini, func = TwoLayer, parms = parameters, method = 'rk4')

plot(out[,1], out[,2], col = 'red', main = paste('vt_mix =',vt_mix))
lines(out[,1], out[,3], col = 'blue')

result <- data.frame('Time' = out[,1],
                     'WT_epi' = out[,2], 'WT_hyp' = out[,3])
g1 <- ggplot(result) +
  geom_line(aes(x=Time, y=WT_epi, col='Surface Mixed Layer')) +
  geom_line(aes(x=(Time), y=WT_hyp, col='Bottom Layer')) +
  labs(x = 'Simulated Time', y = 'WT in deg C')  +
  theme_bw()+
  theme(legend.position="bottom")

output <- read.table('output.txt')
# qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap
output <- data.frame('qin'=output[,1],'qout'=output[,2],'mix_e'=output[,3],'mix_h'=output[,4],
                     'sw'=output[,5],'lw'=output[,6],'water_lw'=output[,7],'conv'=output[,8],
                     'evap'=output[,9], 'Rh' = output[,10],
                     'entrain' = output[,11], 'Ri' = output[,12],'time' = output[,13])
output$balance <- apply(output[,-c(10,11,12,13)],1, sum)

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
  geom_line(aes(x = time,y = Ri, col = 'richardson')) +
  geom_line(aes(x = time,y = entrain, col = 'entrainment')) +
  scale_colour_brewer("Energy terms", palette="Set1") +
  labs(x = 'Simulated Time', y = 'entrainment in cm2/s and Ri in [-]')  +
  theme_bw()+
  theme(legend.position="bottom")


#pdf('Simulation.pdf')
grid.arrange(g1, g2, g3, ncol =1)
#dev.off()

#pdf('Simulation.pdf')
g4 <- grid.arrange(g1, g2, g3, ncol =1)
ggsave(file='2L_visual_result.png', g4, dpi = 300,width = 200,height = 220, units = 'mm')
#dev.off()