rm(list = ls())
library(deSolve)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(thermod)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

if (file.exists('output.txt')) 
  # delete file if it exists
  file.remove('output.txt')

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
Tin <- 10 # inflow water temperature
Q <- 7500 *1e6 # inflow discharge
Rl <- 0.03 # reflection coefficient (generally small, 0.03)
Acoeff <- 0.6 # coefficient between 0.5 - 0.7
sigma <- 11.7 * 10^(-8) # cal / (cm2 d K4) or: 4.9 * 10^(-3) # Stefan-Boltzmann constant in [J (m2 d K4)-1]
eps <- 0.97 # emissivity of water
rho <- 0.9982 # density (g per cm3)
cp <- 0.99 # specific heat (cal per gram per deg C)
c1 <- 0.47 # Bowen's coefficient
a <- 7 # constant
c <- 9e4 # empirical constant
g <- 9.81 # gravity (m/s2)

parameters <- c(Ve, Vh, At, Ht, As, Tin, Q, Rl, Acoeff, sigma, eps, rho, cp, c1, a, c, g)

Et <- 7.07 * 10^(-4)  * ((Ve+Vh)/As/100)^(1.1505) # vertifcal diffusion coefficient (cm2 per d)
vto <- Et/(Ht/100) * (86400/10000) #*100 # heat change coefficient across thermocline during stratified season
vt_mix <- 100 # high heat change coefficient during overturn and mixing events
bound <- as.data.frame(cbind(bound, c(rep(vt_mix,3),rep(vto,6),rep(vt_mix,3))))

colnames(bound) <- c('Month','Jsw','Tair','Dew','vW','vt')
bound$Uw <- 19.0 + 0.95 * (bound$vW * 1000/3600)^2 # function to calculate wind shear stress (and transforming wind speed from km/h to m/s)
bound$vW <- bound$vW * 1000/3600

# boundary <- bound
boundary <- rbind(bound[1,],bound)

boundary <- rbind(boundary, boundary, boundary, boundary, boundary, # running 10 years
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

boundary <- add_noise(boundary)

times <- seq(from = 1, to = max(boundary$Month), by = 1)
yini <- c(5,5) # initial water temperatures

out <- run_model(bc = boundary, params = parameters, ini = yini, times = times)

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
ggsave(file='images/2L_visual_result.png', g5, dpi = 300,width = 200,height = 220, units = 'mm')