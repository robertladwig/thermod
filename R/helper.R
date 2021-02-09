## Helper functions

#' Calculate water density from temperature
#'
#' Calculate water density from water temperature using the formula from (Millero & Poisson, 1981).
#'
#' @param wtemp vector or matrix; Water temperatures
#' @return vector or matrix; Water densities in kg/m3
#' @export
calc_dens <- function(wtemp){
  dens = 999.842594 + (6.793952 * 10^-2 * wtemp) - (9.095290 * 10^-3 *wtemp^2) + (1.001685 * 10^-4 * wtemp^3) - (1.120083 * 10^-6* wtemp^4) + (6.536336 * 10^-9 * wtemp^5)
  return(dens)
}

#' Add noise to boundary conditions
#'
#' Adds random noise to data.
#'
#' @param bc vector or matrix; meteorological boundary conditions
#' @return vector or matrix; meteorological boundary conditions with noise
#' @export
add_noise <- function(bc){
  for (ii in 2:ncol(bc)){
    corrupt <- rbinom(nrow(bc),1,0.1)    # choose an average of 10% to corrupt at random
    corrupt <- as.logical(corrupt)
    noise <- rnorm(sum(corrupt),mean = mean(bc[,ii]),sd=sd(bc[,ii])/10)
    bc[corrupt,ii] <- bc[corrupt,ii] +  noise
  }
  return(bc)
}


#' Configure model with LakeEnsemblR data
#'
#' Use existing configuration data from LakeEnsemblR to prepare configuration files for the model
#'
#' @param config_file string of your LER configuration file
#' @param folder folder path to the configuration file
#' @return vector of parameters to run thermod
#' @export
#' @import LakeEnsemblR 
#' @import gotmtools
configure_from_ler = function(config_file = 'LakeEnsemblR.yaml', folder = '.'){
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
  met_file <- get_yaml_value(file = config_file, label = "meteo", key = "file")
  met <- read.csv(met_file, stringsAsFactors = F)
  met[, 1] <- as.POSIXct(met[, 1])
  tstep <- diff(as.numeric(met[, 1]))
  met$Dewpoint_Air_Temperature_Celsius <- met$Air_Temperature_celsius - ((100 - met$Relative_Humidity_percent)/5.)
  bound <- met %>%
    dplyr::select(datetime, Shortwave_Radiation_Downwelling_wattPerMeterSquared,
                  Air_Temperature_celsius, Dewpoint_Air_Temperature_Celsius, Ten_Meter_Elevation_Wind_Speed_meterPerSecond)
  bound$Shortwave_Radiation_Downwelling_wattPerMeterSquared <- bound$Shortwave_Radiation_Downwelling_wattPerMeterSquared * 2.0E-5 *3600*24
  bound <- bound %>%
    dplyr::filter(datetime >= start_date & datetime <= stop_date)
  bound = bound %>%
    group_by(datetime) %>%
    summarise_all(mean)
  bound$datetime <- seq(1, nrow(bound))
  write_delim(bound, path = paste0(folder,'/meteo.txt'), delim = '\t')
  
  parameters <- c(Ve, Vh, At, Ht, As, Tin, Q, Rl, Acoeff, sigma, eps, rho, cp, c1, a, c, g, simple_therm_depth, calParam)
  
  return(parameters)
}

#' Simple two-layer water temperature model
#'
#' Runs a simple two-layer water temperature model in dependence of meteorological drivers and
#' lake setup/configuration. The model automatically calculates summer stratification onset and
#' offset
#'
#' @param bc meteorological boundary conditions: day, shortwave radiation, air temperature, dew
#' point temperature, wind speed and wind shear stress
#' @param params configuration parameters 
#' @param ini vector of the initial water temperatures of the epilimnion and hypolimnion
#' @param times vector of time information
#' @return matrix of simulated water temperatures in the pilimnion and hypolimnion
#' @export
#' @import deSolve 
#' @import LakeMetabolizer
run_model <- function(bc, params, ini, times){
  Ve <- params[1] # epilimnion volume (cm3)
  Vh <- params[2] # hypolimnion volume (cm3)
  At <- params[3] # thermocline area (cm2)
  Ht <- params[4] # thermocline thickness (cm)
  As <- params[5] # surface area (cm2)
  Tin <- params[6] # inflow water temperature (deg C)
  Q <- params[7] # inflow discharge (deg C)
  Rl <- params[8] # reflection coefficient (generally small, 0.03)
  Acoeff <- params[9] # coefficient between 0.5 - 0.7
  sigma <- params[10] # cal / (cm2 d K4) or: 4.9 * 10^(-3) # Stefan-Boltzmann constant in (J (m2 d K4)-1)
  eps <- params[11] # emissivity of water
  rho <- params[12] # density (g per cm3)
  cp <- params[13] # specific heat (cal per gram per deg C)
  c1 <- params[14] # Bowen's coefficient
  a <- params[15] # constant
  c <- params[16] # empirical constant
  g <- params[17] # gravity (m/s2)
  thermDep <- params[18] # thermocline depth (cm)
  calParam <- params[19] # multiplier coefficient for calibrating the entrainment over the thermocline depth
  
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
      dV = 100 * calParam
    } else {
      dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000) ** calParam
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
# run_model <- function(modelfunc = 'TwoLayer', bc, params, ini, times){
#   
#   Ve <- params[1]
#   Vh <- params[2]
#   At <- params[3]
#   Ht <- params[4]
#   As <- params[5]
#   Tin <- params[6]
#   Q <- params[7]
#   Rl <- params[8]
#   Acoeff <- params[9]
#   sigma <- params[10]
#   eps <- params[11]
#   rho <- params[12]
#   cp <- params[13]
#   c1 <- params[14]
#   a <- params[15]
#   c <- params[16]
#   g <- params[17]
#   NEP <- params[18]
#   Fsed<- params[19]
#   MINERAL <- params[20]
#   Ased <- params[21]
#   diffred <- params[22]
#   
#   TwoLayer <- function(t, y, parms){
#   eair <- (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t)))) # air vapor pressure
#   esat <- 4.596 * exp((17.27 * Tair(t)) / (237.3 + Tair(t))) # saturation vapor pressure
#   RH <- eair/esat *100 # relative humidity
#   es <- 4.596 * exp((17.27 * y[1])/ (273.3+y[1]))
#   # diffusion coefficient
#   Cd <- 0.00052 * (vW(t))^(0.44)
#   shear <- 1.164/1000 * Cd * (vW(t))^2
#   #  shear <- Uw(t)
#   rho_e <- calc_dens(y[1])/1000
#   rho_h <- calc_dens(y[2])/1000
#   w0 <- sqrt(shear/rho_e) # sqrt(shear/rho)
#   E0  <- c * w0
#   Ri <- ((g/rho)*(abs(rho_e-rho_h)/10))/(w0/(10)^2)
#   if (rho_e > rho_h){
#     dV = 100
#   } else {
#     dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000)
#   }
#   
#   # epilimnion water temperature change per time unit
#   dTe <-  Q / Ve * Tin -              # inflow heat
#     Q / Ve * y[1] +                   # outflow heat
#     ((dV * At) / Ve) * (y[2] - y[1]) + # mixing between epilimnion and hypolimnion # ((vt(t) * At) / Ve) * (y[2] - y[1])
#     + As/(Ve * rho * cp) * (
#       Jsw(t)  + # shortwave radiation
#         (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
#         (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
#         (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
#         (Uw(t) * ((es) - (eair))) )# evaporation
#   
#   # hypolimnion water temperature change per time unit
#   dTh <-  ((dV * At) / Vh) * (y[1] - y[2]) # ((vt(t) * At) / Vh) * (y[1] - y[2]) 
#   
#   
#   
#   
#   # diagnostic variables for plotting
#   mix_e <- (vt(t) * At * rho * cp) * (y[2] - y[1])  /As#((vt(t) * At) / Ve) * (y[2] - y[1])
#   mix_h <- (vt(t) * At *rho * cp) * (y[1] - y[2]) /At#(vt(t) * At) / (Vh) * (y[1] - y[2])
#   qin <- Q * rho * cp * Tin /As#Q / Ve * Tin
#   qout <- - Q * rho * cp * y[1] /As #- Q / Ve * y[1] 
#   sw <- Jsw(t) #*As#* As/(Ve * rho * cp)
#   lw <- (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) #* As#* As/(Ve * rho * cp)
#   water_lw <- - (eps * sigma * (y[1]+ 273)^4) #*As#*As/(Ve * rho * cp)
#   conv <- - (c1 * Uw(t) * (y[1] - Tair(t))) #*As#*As/(Ve * rho * cp)
#   evap <- - (Uw(t) * ((esat) - (eair))) #*As #*As/(Ve * rho * cp)
#   Rh <- RH
#   E <- (E0 / (1 + a * Ri)^(3/2))
#   
#   write.table(matrix(c(qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap, Rh,E, Ri, t), nrow=1), 'output.txt', append = TRUE,
#               quote = FALSE, row.names = FALSE, col.names = FALSE)
#   # write.table(matrix(c(y[1], y[2], dTe, dTh, t), nrow=1), 'test.txt', append = TRUE,
#   #             quote = FALSE, row.names = FALSE, col.names = FALSE)
#   
#   return(list(c(dTe, dTh)))
#   }
#   
#   TwoLayerOxy <- function(t, y, parms){
#     eair <- (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t)))) # air vapor pressure
#     esat <- 4.596 * exp((17.27 * Tair(t)) / (237.3 + Tair(t))) # saturation vapor pressure
#     RH <- eair/esat *100 # relative humidity
#     es <- 4.596 * exp((17.27 * y[1])/ (273.3+y[1]))
#     # diffusion coefficient
#     Cd <- 0.00052 * (vW(t))^(0.44)
#     shear <- 1.164/1000 * Cd * (vW(t))^2
#     #  shear <- Uw(t)
#     rho_e <- calc_dens(y[1])/1000
#     rho_h <- calc_dens(y[2])/1000
#     w0 <- sqrt(shear/rho_e) # sqrt(shear/rho)
#     E0  <- c * w0
#     Ri <- ((g/rho)*(abs(rho_e-rho_h)/10))/(w0/(10)^2)
#     if (rho_e > rho_h){
#       dV = 1e3 #100
#       dV_oxy <- dV/diffred#/1e2
#       mult = 1/1000
#     } else {
#       dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000)
#       dV_oxy <- dV/diffred
#       mult = 1.0
#     }
#     
#     # epilimnion water temperature change per time unit
#     dTe <-  Q / Ve * Tin -              # inflow heat
#       Q / Ve * y[1] +                   # outflow heat
#       ((dV * At) / Ve) * (y[2] - y[1]) + # mixing between epilimnion and hypolimnion # ((vt(t) * At) / Ve) * (y[2] - y[1])
#       + As/(Ve * rho * cp) * (
#         Jsw(t)  + # shortwave radiation
#           (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
#           (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
#           (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
#           (Uw(t) * ((es) - (eair))) )# evaporation
#     
#     # hypolimnion water temperature change per time unit
#     dTh <-  ((dV * At) / Vh) * (y[1] - y[2]) # ((vt(t) * At) / Vh) * (y[1] - y[2]) 
#     
#     U10 <- wind.scale.base(vW(t), wnd.z = 10)
#     K600 <- k.cole.base(U10)
#     water.tmp = y[1]
#     kO2 <- k600.2.kGAS.base(k600=K600, 
#                                  temperature=water.tmp, 
#                                  gas='O2') * 100 # m/d * 100 cm/m
#     o2sat<-o2.at.sat.base(temp= water.tmp, 
#                                altitude = 300)/1000 # mg O2/L
#     Fatm <- kO2*(o2sat  - y[3] /Ve) * (As) # mg * cm/d * cm2/cm3= cmg/d
#     #  kO2*(o2sat - y[3] )/ (As/Ve) # mg/m m/d = mg/d
#     
#     Sed <- Fsed * y[4] * (Ased/Vh)  * 1.03^(y[2]-20) * mult # m/d * mg * m2/m3
#     #  Fsed * y[4]/ (Ased/Vh)  * 1.08^(y[2]-20) * mult
#     
#     PP <- 1.03^(y[1]-20) * NEP * Ve * mult # mg/m3/d * m3
#     # 1.08^(y[1]-20) * NEP * Ve * mult 
#     
#     VOL <- 1.03^(y[2]-20) * MINERAL * Vh * mult # mg/m3/d * m3
#     
#     dOe <-( PP +
#       Fatm +
#       ((dV_oxy  * At)) * (y[4]/Vh - y[3]/Ve)) 
#     
#     dOh <- ( ((dV_oxy * At)) * (y[3]/Ve - y[4]/Vh) - #((dV_oxy * At) / Vh) * (y[3] - y[4])
#                - VOL - 
#       Sed) 
#     
#     
#     # diagnostic variables for plotting
#     mix_e <- (vt(t) * At * rho * cp) * (y[2] - y[1])  /As#((vt(t) * At) / Ve) * (y[2] - y[1])
#     mix_h <- (vt(t) * At *rho * cp) * (y[1] - y[2]) /At#(vt(t) * At) / (Vh) * (y[1] - y[2])
#     qin <- Q * rho * cp * Tin /As#Q / Ve * Tin
#     qout <- - Q * rho * cp * y[1] /As #- Q / Ve * y[1] 
#     sw <- Jsw(t) #*As#* As/(Ve * rho * cp)
#     lw <- (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) #* As#* As/(Ve * rho * cp)
#     water_lw <- - (eps * sigma * (y[1]+ 273)^4) #*As#*As/(Ve * rho * cp)
#     conv <- - (c1 * Uw(t) * (y[1] - Tair(t))) #*As#*As/(Ve * rho * cp)
#     evap <- - (Uw(t) * ((esat) - (eair))) #*As #*As/(Ve * rho * cp)
#     Rh <- RH
#     E <- (E0 / (1 + a * Ri)^(3/2))
#     Oflux_epi <- ((dV_oxy  * At)) * (y[4]/Vh - y[3]/Ve)
#     Oflux_hypo <- ((dV_oxy * At)) * (y[3]/Ve - y[4]/Vh)
#     
#     write.table(matrix(c(qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap, Rh,E, Ri, t,
#                          Fatm, Sed, PP, VOL, Oflux_epi, Oflux_hypo), nrow=1), 'output.txt', append = TRUE,
#                 quote = FALSE, row.names = FALSE, col.names = FALSE)
#     
#     return(list(c(dTe, dTh, dOe, dOh)))
#   }
#   
#   # approximating all boundary conditions 
#   Jsw <- approxfun(x = bc$Month, y = bc$Jsw, method = "linear", rule = 2)
#   Tair <- approxfun(x = bc$Month, y = bc$Tair, method = "linear", rule = 2)
#   Dew <- approxfun(x = bc$Month, y = bc$Dew, method = "linear", rule = 2)
#   Uw <- approxfun(x = bc$Month, y = bc$Uw, method = "linear", rule = 2)
#   vW <- approxfun(x = bc$Month, y = bc$vW, method = "linear", rule = 2)
#   vt <- approxfun(x = bc$Month, y = bc$vt, method = "linear", rule = 2)
# 
#   # runge-kutta 4th order
#   if (modelfunc == 'TwoLayer'){
#     out <- ode(times = times, y = ini, func = TwoLayer, parms = params, method = 'rk4')
#   } else if (modelfunc == 'TwoLayerOxy'){
#     out <- ode(times = times, y = ini, func = TwoLayerOxy, parms = params, method = 'rk4')
#   }
#  
#   
#   return(out)
# }
