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
parameters <- configure_from_ler(config_file = config_file, folder = folder)

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
                     'WT_epi' = out[,2], 'WT_hyp' = out[,3],
                     'AT' = bound$Tair)

ggplot(result) +
  geom_line(aes(Time, WT_epi, col = 'WT_epi')) +
  geom_line(aes(Time, AT, col = 'AT')) +
  theme_minimal()


# https://www.kalmanfilter.net/kalman1d.html
state = rep(NA, length(result$Time))
uncertainty = state
for (i in 1:length(result$Time)){
  if (i == 1) {
    x = result$WT_epi[1]
    sigma0 = 10
    est_p = sigma0^2
    q = 1e-4
  } else {
    x = x
    est_p = est_p + q
    meas_p = 0.1^2
    K = est_p / (est_p + meas_p)
    x = x + K * (result$WT_epi[i] - x)
    est_p = (1 - K) * est_p
  }
  state[i] = x
  uncertainty[i] = est_p
}

kalman.df <- data.frame('Time' = out[,1],
                                 'WT_epi' = out[,2], 'WT_hyp' = out[,3],
                                 'AT' = bound$Tair,
                                 'WT_kalman' = state,
                                 'Uncertainty' = uncertainty)

ggplot(kalman.df) +
  geom_line(aes(Time, WT_epi, col = 'Simulated WT')) +
  geom_line(aes(Time, WT_kalman, col = 'Kalman'), linetype = 'dashed') +
  theme_minimal()

ggplot(kalman.df) +
  geom_line(aes(Time, Uncertainty))


# x - state
# p - error / covariance
# a - state transition
# q - transition error
# z - measurement
# h - state-to-measurement transformation
# y difference
# kg - kalman gain
# r - measurement error
# https://stackoverflow.com/questions/33384112/kalman-filter-one-dimensional-several-approaches
x = rep(NA, length(result$Time))
p = state
for (k in 1:length(result$Time)){
  if (k == 1){
    x[k] = result$WT_epi[1]
    p[k] = 10^2
  } else {
    # a = rnorm(1, 0, 1e-2)
    a = 1
    # q = 1e-1
    q = 0
    z = result$WT_epi[k]
    # h = 1e-2
    h = 1
    r = 1e1
    # prediction
    x[k] = a * x[k - 1]
    p[k] = a * p[k - 1] * a + q
    # update
    y = z - h * x[k]
    kg = p * h / (h * p * h + r)
    x[k] = x[k] + kg * y
    p[k] = (1 - kg * h) * p[k]
  }
}

kalman.df2 = data.frame('Time' = out[,1],
                                 'WT_epi' = out[,2], 'WT_hyp' = out[,3],
                                 'AT' = bound$Tair,
                                 'WT_kalman' = x,
                                 'Uncertainty' = p)

ggplot(kalman.df2) +
  geom_line(aes(Time, WT_epi, col = 'Simulated WT')) +
  geom_line(aes(Time, WT_kalman, col = 'Kalman'), linetype = 'dashed') +
  theme_minimal()

ggplot(kalman.df2) +
  geom_line(aes(Time, Uncertainty))
