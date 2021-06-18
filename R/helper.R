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
    noise <- rnorm(sum(corrupt), mean = mean(as.numeric(unlist(bc[,ii])), na.rm = TRUE),
                   sd = sd(as.numeric(unlist(bc[,ii])),na.rm = TRUE)/10)
    bc[corrupt,ii] <- bc[corrupt,ii] +  noise
  }
  return(bc)
}

#' Smoothing of time series
#'
#' Filters the signal using a Kalman filter
#'
#' @param bc vector or matrix; meteorological boundary conditions
#' @return vector or matrix; meteorological boundary conditions with noise
#' @export
kalman_filtering <- function(time, series){
  state = rep(NA, length(time))
  uncertainty = state
  for (i in 1:length(time)){
    if (i == 1) {
      x = series[1]
      sigma0 = 10
      est_p = sigma0^2
      q = 1e-4
    } else {
      x = x
      est_p = est_p + q
      meas_p = 0.1^2
      K = est_p / (est_p + meas_p)
      x = x + K * (series[i] - x)
      est_p = (1 - K) * est_p
    }
    state[i] = x
    uncertainty[i] = est_p
  }
  return(state)
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
configure_from_ler = function(config_file = 'LakeEnsemblR.yaml', folder = '.', mode = NULL){
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
  
  if (length(mode) > 0){
    if (mode == 'forecast'){
      ix = max_depth
    } else {
      ix <- match(simple_therm_depth, hyp$Depth_meter)
    } 
  } else {
    ix <- match(simple_therm_depth, hyp$Depth_meter)
  }
  
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
    mutate(date = as.Date(datetime)) %>%
    group_by(date) %>%
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
#' @param ice boolean, if TRUE the model will approximate ice on/off set, otherwise no ice
#' @return matrix of simulated water temperatures in the epilimnion and hypolimnion
#' @export
#' @import deSolve 
#' @import LakeMetabolizer
run_model <- function(bc, params, ini, times, ice = FALSE){
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
    if (ice == TRUE){
      if (y[1] <= 0 && Tair(t) <= 0){
        ice_param = 1e-5
      } else {
        ice_param = 1
      }
      # epilimnion water temperature change per time unit
      dTe <-  Q / Ve * Tin -              # inflow heat
        Q / Ve * y[1] +                   # outflow heat
        ((dV * At) / Ve) * (y[2] - y[1]) + # mixing between epilimnion and hypolimnion
        + As/(Ve * rho * cp) * ice_param * (
          Jsw(t)  + # shortwave radiation
            (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
            (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
            (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
            (Uw(t) * ((es) - (eair))) )# evaporation
    } else{
    
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
    }
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
    
    write.table(matrix(c(qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap, Rh,E, Ri, t, ice_param), nrow=1), 
                'output.txt', append = TRUE,
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


#' Simple two-layer water temperature and oxygen model
#'
#' Runs a simple two-layer water temperature and oxygen model in dependence of meteorological 
#' drivers and lake setup/configuration. The model automatically calculates summer stratification 
#' onset and offset
#'
#' @param bc meteorological boundary conditions: day, shortwave radiation, air temperature, dew
#' point temperature, wind speed and wind shear stress
#' @param params configuration parameters 
#' @param ini vector of the initial water temperatures of the epilimnion and hypolimnion
#' @param times vector of time information
#' @param ice boolean, if TRUE the model will approximate ice on/off set, otherwise no ice
#' @return matrix of simulated water temperature and oxygen mass in the epilimnion and hypolimnion
#' @export
#' @import deSolve 
#' @import LakeMetabolizer
run_oxygen_model <- function(bc, params, ini, times, ice = FALSE){
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
  Fnep <- params[20]
  Fsed<- params[21]
  Ased <- params[22]
  diffred <- params[23]
  
  TwoLayer_oxy <- function(t, y, parms){
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
      dV_oxy = dV/diffred
      mult = 1.#1/1000
    } else {
      dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000) ** calParam
      dV_oxy = dV/diffred
      mult = 1
    }
    U10 <- wind.scale.base(vW(t), wnd.z = 10)
    K600 <- k.cole.base(U10)
    water.tmp = y[1]
    kO2 <- k600.2.kGAS.base(k600=K600,
                            temperature=water.tmp,
                            gas='O2') * 100 # m/d * 100 cm/m
    o2sat<-o2.at.sat.base(temp= water.tmp,
                          altitude = 300)/1e3 # mg/L ==> g/m3 ==> mg O2 * 1000/cm3 * 1e6
    
    if (ice == TRUE){
      if (y[1] <= 0 && Tair(t) <= 0){
        ice_param = 1e-5
      } else {
        ice_param = 1
      }
      # epilimnion water temperature change per time unit
      dTe <-  Q / Ve * Tin -              # inflow heat
        Q / Ve * y[1] +                   # outflow heat
        ((dV * At) / Ve) * (y[2] - y[1]) + # mixing between epilimnion and hypolimnion
        + As/(Ve * rho * cp) * ice_param * (
          Jsw(t)  + # shortwave radiation
            (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
            (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
            (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
            (Uw(t) * ((es) - (eair))) )# evaporation
      ATM <- kO2*(o2sat  - y[3] /Ve) * (As) * ice_param# mg * cm/d * cm2/cm3= mg/d
    } else{
      
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
      ATM <- kO2*(o2sat  - y[3] /Ve) * (As) # mg/cm3 * cm/d * cm2= mg/d
    }
    # hypolimnion water temperature change per time unit
    dTh <-  ((dV * At) / Vh) * (y[1] - y[2]) 
    
    
        SED <- Fsed * Ased * 1.03^(y[2]-20) * (y[4]/Vh/(0.5/1000 + y[4]/Vh))* mult # mg/m2/d * m2 

        NEP <- 1.03^(y[1]-20) * Fnep * Ve * mult # mg/m3/d * m3

        dOe <- ( NEP +
          ATM +
          ((dV_oxy  * At)) * (y[4]/Vh - y[3]/Ve))

        dOh <- (( ((dV_oxy * At)) * (y[3]/Ve - y[4]/Vh) - SED))
    
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
    Oflux_epi <- ((dV_oxy  * At)) * (y[4]/Vh - y[3]/Ve)
    Oflux_hypo <- ((dV_oxy * At)) * (y[3]/Ve - y[4]/Vh)

    write.table(matrix(c(qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap, Rh,E, Ri, t, ice_param,
                         ATM, NEP, SED, Oflux_epi, Oflux_hypo), nrow=1), 
                'output.txt', append = TRUE,
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    return(list(c(dTe, dTh,  dOe, dOh)))
  }
  
  # approximating all boundary conditions  (linear interpolation)
  Jsw <- approxfun(x = bc$Day, y = bc$Jsw, method = "linear", rule = 2)
  Tair <- approxfun(x = bc$Day, y = bc$Tair, method = "linear", rule = 2)
  Dew <- approxfun(x = bc$Day, y = bc$Dew, method = "linear", rule = 2)
  Uw <- approxfun(x = bc$Day, y = bc$Uw, method = "linear", rule = 2)
  vW <- approxfun(x = bc$Day, y = bc$vW, method = "linear", rule = 2)
  
  # runge-kutta 4th order
  out <- ode(times = times, y = ini, func = TwoLayer_oxy, parms = params, method = 'rk4')
  
  return(out)
}

#' Simple two-layer water temperature, oxygen, phytoplankton, zooplankton and nutrient model
#'
#' Runs a simple two-layer water NPZ model in dependence of meteorological 
#' drivers and lake setup/configuration. The model automatically calculates summer stratification 
#' onset and offset
#'
#' @param bc meteorological boundary conditions: day, shortwave radiation, air temperature, dew
#' point temperature, wind speed and wind shear stress
#' @param params configuration parameters 
#' @param ini vector of the initial water temperatures of the epilimnion and hypolimnion
#' @param times vector of time information
#' @param ice boolean, if TRUE the model will approximate ice on/off set, otherwise no ice
#' @return matrix of simulated water temperature and oxygen mass in the epilimnion and hypolimnion
#' @export
#' @import deSolve 
#' @import LakeMetabolizer
run_npz_model <- function(bc, params, ini, times, ice = FALSE){
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
  Fnep <- params[20]
  Fsed<- params[21]
  Ased <- params[22]
  diffred <- params[23]
  #kg, kra, Cgz, ksa, aca, epsilon, krz ksz, apa, apc, Fsedp, alpha1, alpha2
  kg <- params[24] # day-1
  kra <- params[25] # day-1
  Cgz <- params[26] # cm3 per mg C per day
  ksa <- params[27] # day-1
  aca <- params[28] # mg C per mg Chl-a
  epsilon <- params[29] # dimensionless
  krz <- params[30] # day-1
  ksz <- params[31] # day-1
  apa <- params[32] # mg P per mg Chl-a
  apc <- params[33] # mg C per mg P
  Fsedp <- params[34] # mg P per cm3 per day
  alpha1 <- params[35] # dimensionless
  alpha2 <- params[36] # dimensionless
  
  
  TwoLayer_npz <- function(t, y, parms){
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
      dV_oxy = dV/diffred
      mult = 1.#1/1000
    } else {
      dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000) ** calParam
      dV_oxy = dV/diffred
      mult = 1
    }
    U10 <- wind.scale.base(vW(t), wnd.z = 10)
    K600 <- k.cole.base(U10)
    water.tmp = y[1]
    kO2 <- k600.2.kGAS.base(k600=K600,
                            temperature=water.tmp,
                            gas='O2') * 100 # m/d * 100 cm/m
    o2sat<-o2.at.sat.base(temp= water.tmp,
                          altitude = 300)/1e3 # mg/L ==> g/m3 ==> mg O2 * 1000/cm3 * 1e6
    
    if (ice == TRUE){
      if (y[1] <= 0 && Tair(t) <= 0){
        ice_param = 1e-5
      } else {
        ice_param = 1
      }
      # epilimnion water temperature change per time unit
      dTe <-  Q / Ve * Tin -              # inflow heat
        Q / Ve * y[1] +                   # outflow heat
        ((dV * At) / Ve) * (y[2] - y[1]) + # mixing between epilimnion and hypolimnion
        + As/(Ve * rho * cp) * ice_param * (
          Jsw(t)  + # shortwave radiation
            (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
            (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
            (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
            (Uw(t) * ((es) - (eair))) )# evaporation
      ATM <- kO2*(o2sat  - y[3] /Ve) * (As) * ice_param# mg * cm/d * cm2/cm3= mg/d
      
    } else{
      
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
      ATM <- kO2*(o2sat  - y[3] /Ve) * (As) # mg/cm3 * cm/d * cm2= mg/d
    }
    # hypolimnion water temperature change per time unit
    dTh <-  ((dV * At) / Vh) * (y[1] - y[2]) 
    
    
    SED <- Fsed * Ased * 1.03^(y[2]-20) * (y[4]/Vh/(0.5/1000 + y[4]/Vh))* mult # mg/m2/d * m2 
    
    NEP <- 1.03^(y[1]-20) * Fnep * Ve * mult # mg/m3/d * m3
    NEP <- 1.03^(y[1]-20) * (alpha1 * kg  * (y[9])/(y[9] + 2 * Ve /10e6) - alpha2 * kra) * y[5] * Ve * mult #
    
    dOe <- ( NEP +
               ATM +
               ((dV_oxy  * At)) * (y[4]/Vh - y[3]/Ve))
    
    dOh <- (( ((dV_oxy * At)) * (y[3]/Ve - y[4]/Vh) - SED))
    
    # 5 Pe, 6 Ph, 7 Ze, 8 Zh, 9 Ne, 10 Nh
    # kg, kra, Cgz, ksa, aca, epsilon, krz ksz, apa, apc, Fsedp, alpha1, alpha2
    dPe <- ((kg  * (y[9])/(y[9] + 2 * Ve /10e6) - kra) * y[5] - Cgz / Ve * y[7] * y[5]  ) * 1.03^(y[1]-20) - ksa * y[5] 
    dZe <- ((aca * epsilon * Cgz / Ve) * y[7] * y[5] - krz * y[7])* 1.03^(y[1]-20)  - ksz * y[7]
    dNe <-( apa * (1- epsilon) * Cgz / Ve * y[7] * y[5] + apc * krz * y[7] - apa * (kg * (y[9])/(y[9] + 2 * Ve /10e6) - kra) * y[5]) * 1.03^(y[1]-20)+
      ((dV_oxy  * At)) * (y[10]/Vh - y[9]/Ve)
    
    dPh <- ((kg * (y[10])/(y[10] + 2 * Vh /10e6) - kra) * y[6] - Cgz / Vh * y[8] * y[6]) * 1.03^(y[2]-20) + ksa * y[5] - ksa * y[6] 
    dZh <- ((aca * epsilon * Cgz / Vh) * y[8] * y[6] - krz * y[8]) * 1.03^(y[2]-20) + ksz * y[7] - ksz * y[8]
    dNh <- (apa * (1- epsilon) * Cgz / Vh * y[8] * y[6] + apc * krz * y[8] - apa * (kg * (y[10])/(y[10] + 2 * Vh /10e6)  - kra) * y[6]) * 1.03^(y[2]-20)  +
      ((dV_oxy  * At)) * (y[9]/Ve - y[10]/Vh) -
      Fsedp * Ased * 1.03^(y[2]-20) * (y[4]/Vh/(0.5/1000 + y[4]/Vh))* mult
    
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
    Oflux_epi <- ((dV_oxy  * At)) * (y[4]/Vh - y[3]/Ve)
    Oflux_hypo <- ((dV_oxy * At)) * (y[3]/Ve - y[4]/Vh)
    
    write.table(matrix(c(qin, qout, mix_e, mix_h, sw, lw, water_lw, conv, evap, Rh,E, Ri, t, ice_param,
                         ATM, NEP, SED, Oflux_epi, Oflux_hypo), nrow=1), 
                'output.txt', append = TRUE,
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    return(list(c(dTe, dTh,  dOe, dOh, dPe, dPh, dZe, dZh, dNe, dNh)))
  }
  
  # approximating all boundary conditions  (linear interpolation)
  Jsw <- approxfun(x = bc$Day, y = bc$Jsw, method = "linear", rule = 2)
  Tair <- approxfun(x = bc$Day, y = bc$Tair, method = "linear", rule = 2)
  Dew <- approxfun(x = bc$Day, y = bc$Dew, method = "linear", rule = 2)
  Uw <- approxfun(x = bc$Day, y = bc$Uw, method = "linear", rule = 2)
  vW <- approxfun(x = bc$Day, y = bc$vW, method = "linear", rule = 2)
  
  # runge-kutta 4th order
  out <- ode(times = times, y = ini, func = TwoLayer_npz, parms = params, method = 'rk4')
  
  return(out)
}



###FORECASTING WORK
#' Simple one-layer water temperature and oxygen model in forecast mode
#'
#' Runs a simple one-layer water temperature and oxygen model in dependence of meteorological 
#' drivers and lake setup/configuration. 
#'
#' @param bc meteorological boundary conditions: day, shortwave radiation, air temperature, dew
#' point temperature, wind speed and wind shear stress
#' @param params configuration parameters 
#' @param ini vector of the initial water temperatures of the epilimnion and hypolimnion
#' @param times vector of time information
#' @param ice boolean, if TRUE the model will approximate ice on/off set, otherwise no ice
#' @return matrix of simulated water temperature and oxygen mass in the epilimnion and hypolimnion
#' @export
#' @import deSolve 
#' @import LakeMetabolizer
run_temp_oxygen_forecast <- function(bc, params, ini, times, ice = FALSE, 
                                     observed){
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
  Fnep <- params[20]
  Fsed<- params[21]
  Ased <- params[22]
  diffred <- params[23]
  
  OneLayer_forecast <- function(t, y, parms){
    Q <- parms[7]
    Fnep <- parms[20]
    eair <- (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t)))) # air vapor pressure
    esat <- 4.596 * exp((17.27 * Tair(t)) / (237.3 + Tair(t))) # saturation vapor pressure
    RH <- eair/esat *100 # relative humidity
    es <- 4.596 * exp((17.27 * y[1])/ (273.3+y[1]))
    Cd <- 0.00052 * (vW(t))^(0.44)
    mult = 1
    U10 <- wind.scale.base(vW(t), wnd.z = 10)
    K600 <- k.cole.base(U10)
    water.tmp = y[1]
    kO2 <- k600.2.kGAS.base(k600=K600,
                            temperature=water.tmp,
                            gas='O2') * 100 # m/d * 100 cm/m
    o2sat<-o2.at.sat.base(temp= water.tmp,
                          altitude = 300)/1e3 # mg/L ==> g/m3 ==> mg O2 * 1000/cm3 * 1e6
    
    if (ice == TRUE){
      if (y[1] <= 0 && Tair(t) <= 0){
        ice_param = 1e-5
      } else {
        ice_param = 1
      }
      # epilimnion water temperature change per time unit
      dTe <- Q / Ve * y[1] + As/(Ve * rho * cp) * ice_param * (
          Jsw(t)  + # shortwave radiation
            (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
            (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
            (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
            (Uw(t) * ((es) - (eair))) )# evaporation
      ATM <- kO2*(o2sat  - y[2] /Ve) * (As) * ice_param# mg * cm/d * cm2/cm3= mg/d
    } else{
      
      # epilimnion water temperature change per time unit
      dTe <-  Q / Ve * y[1] + As/(Ve * rho * cp) * (
          Jsw(t)  + # shortwave radiation
            (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl)) - # longwave radiation into the lake
            (eps * sigma * (y[1] + 273)^4)  - # backscattering longwave radiation from the lake
            (c1 * Uw(t) * (y[1] - Tair(t))) - # convection
            (Uw(t) * ((es) - (eair))) )# evaporation
      ATM <- kO2*(o2sat  - y[2] /Ve) * (As) # mg/cm3 * cm/d * cm2= mg/d
    }
    # hypolimnion water temperature change per time unit
    
    

    
    NEP <- 1.03^(y[1]-20) * Fnep * Ve * mult # mg/m3/d * m3
    
    dOe <- ( NEP +
               ATM )
    
    # diagnostic variables for plotting
    sw <- Jsw(t) 
    lw <- (sigma * (Tair(t) + 273)^4 * (Acoeff + 0.031 * sqrt(eair)) * (1 - Rl))
    water_lw <- - (eps * sigma * (y[1]+ 273)^4)
    conv <- - (c1 * Uw(t) * (y[1] - Tair(t))) 
    evap <- - (Uw(t) * ((esat) - (eair)))
    Rh <- RH
    # E <- (E0 / (1 + a * Ri)^(3/2))
    
    write.table(matrix(c(sw, lw, water_lw, conv, evap, Rh, t, ice_param,
                         ATM, NEP), nrow=1), 
                'output.txt', append = TRUE,
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    return(list(c(dTe,  dOe)))
  }
  
  # approximating all boundary conditions  (linear interpolation)
  Jsw <- approxfun(x = bc$Day, y = bc$Jsw, method = "linear", rule = 2)
  Tair <- approxfun(x = bc$Day, y = bc$Tair, method = "linear", rule = 2)
  Dew <- approxfun(x = bc$Day, y = bc$Dew, method = "linear", rule = 2)
  Uw <- approxfun(x = bc$Day, y = bc$Uw, method = "linear", rule = 2)
  vW <- approxfun(x = bc$Day, y = bc$vW, method = "linear", rule = 2)
  
  # runge-kutta 4th order
  # out <- ode(times = times, y = ini, func = OneLayer_forecast, parms = params, method = 'rk4')
  
  out_total <- c()
  idx <- which(!is.na(observed$WT_obs))
  idstart <- 1
  q_par <- params[7]
  nep_par <- params[20]
 
  for (nn in which(!is.na(observed[,2]) | !is.na(observed[,3]))[2:length( which(!is.na(observed[,2]) | !is.na(observed[,3])))]){
    #which(!is.na(observed$WT_obs))[2:length(which(!is.na(observed$WT_obs)))]
    # which(!is.na(observed[,2:3]))[2:length(which(!is.na(observed[,2:3])))]
    idstop = nn
    result <- matrix(NA, nrow = 100, ncol =2)
    param_matrix <- matrix(NA, nrow = 100, ncol =2)
    id_row1 = 1
    id_row2 = 1
    
    if ( !is.na(observed$WT_obs[idstop])){
      for (random_run in 1:100){
        params[7] <- rnorm(n = 1, mean = q_par, sd = 1e7)
        while(is.na(params[7])){
          params[7] <- rnorm(n = 1, mean = q_par, sd =  1e7)
        }
        out <- ode(times = c(idstart:idstop), y = ini, func = OneLayer_forecast, parms = params, method = 'rk4')
        nrmse_temp <- sqrt(sum((out[,2]- observed$WT_obs[idstart:idstop])^2, na.rm = TRUE)/(nrow(out)))/(max(out[,2]) - min(out[,2]))
        result[random_run,1] <- nrmse_temp 
        param_matrix[random_run,1] <- params[7] 
      }
      
      id_row1 <- which.min(result[,1])
      params[7] <- param_matrix[id_row1, 1]
      q_par = params[7]
      
    }

    if (!is.na(observed$DO_obs[idstop])){
      for (random_run2 in 1:100){
        params[20] <-  rnorm(n = 1, mean = nep_par, sd =  1e-5)
        while(is.na(params[20])){
          params[20] <-  rnorm(n = 1, mean = nep_par, sd = 1e-5)
        }
        out <- ode(times = c(idstart:idstop), y = ini, func = OneLayer_forecast, parms = params, method = 'rk4')
        nrmse_do <- sqrt(sum((out[,3]/ 1000 /  params[1] * 1e6 - observed$DO_obs[idstart:idstop])^2, na.rm = TRUE)/nrow(out))/((max(out[,3]) - min(out[,3]))/ 1000 /  params[1] * 1e6 )    
        
        result[random_run2,2] <- nrmse_do #, nrmse_do)
        param_matrix[random_run2,2] <- params[20] # params[20])
      }
      
      id_row2 <- which.min(result[,2])
      params[20] <- param_matrix[id_row2, 2]
      nep_par <- params[20]
    }

    # params[7] <- param_matrix[id_row1, 1]
    # 
    # if (length(id_row2) > 0){
    #   
    # } else {
    #   params[20] <- param_matrix[id_row1, 2]
    #   nep_par <- params[20]
    # }
    
    out <- ode(times = c(idstart:idstop), y = ini, func = OneLayer_forecast, parms = params, method = 'rk4')
    
    if (match(nn, which(!is.na(observed[,2]) | !is.na(observed[,3]))[2:length( which(!is.na(observed[,2]) | !is.na(observed[,3])))]) == 1){
      out_total <- rbind(out_total, out) 
    } else {
      out_total <- rbind(out_total, out[-c(1),]) 
    }
    
    print(paste0(match(nn, which(!is.na(observed[,2]) | !is.na(observed[,3]))[2:length( which(!is.na(observed[,2]) | !is.na(observed[,3])))]),'/',
                 length(which(!is.na(observed[,2]) | !is.na(observed[,3]))[2:length( which(!is.na(observed[,2]) | !is.na(observed[,3])))]),
                 '; current WTR NRMSE: ',round(result[id_row1,1],5),
                 '; current DO NRMSE: ',round(result[id_row2,2],5)))
    # print(q_par)
    # print(nep_par)
    
    idstart = nn
    ini <- out[nrow(out), 2:3]
    
    if (!is.na(observed$WT_obs[idstop])){
    ini[1] <- rnorm(n = 1, mean = observed$WT_obs[idstop], sd = 0.1)
    }
    if (!is.na(observed$DO_obs[idstop])){
    ini[2] <- rnorm(n = 1, mean = observed$DO_obs[idstop] * 1000 *  params[1] / 1e6, sd = 10)
    # 1000/1e6  * wq_parameters[1]
    }
  }
  
  if (max(times) > idstart){
    out_forecast <- list()
    for (random_run3 in 1:100){
      params[7] <- rnorm(n = 1, mean = q_par, sd = 1e7)
      while(is.na(params[7])){
        params[7] <- rnorm(n = 1, mean = q_par, sd =  1e7)
      }
      params[20] <-  rnorm(n = 1, mean = nep_par, sd =  1e-5)
      while(is.na(params[20])){
        params[20] <-  rnorm(n = 1, mean = nep_par, sd = 1e-5)
      }
      out <- ode(times = c(idstart:max(times)), y = ini, func = OneLayer_forecast, parms = params, method = 'rk4')
      out_df <- rbind(out_total,  out[-c(1),])
      
      out_forecast[[match(random_run3, seq(1,100))]] <- out_df
    }
    
  } else { out_forecast = out_total}
  
  return(out_forecast)
}



#' retreive the model time steps based on start and stop dates and time step
#'
#' @param model_start model start date in date class
#' @param model_stop model stop date in date class
#' @param time_step model time step, defaults to daily timestep
get_model_dates = function(model_start, model_stop, time_step = 'days'){
  
  model_dates = seq.Date(from = as.Date(model_start), to = as.Date(model_stop), by = time_step)
  
  return(model_dates)
}

#' vector for holding states and parameters for updating
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_params_est number of parameters we're calibrating
#' @param n_step number of model timesteps
#' @param n_en number of ensembles
get_Y_vector = function(n_states, n_params_est, n_step, n_en){
  
  Y = array(dim = c(n_states + n_params_est, n_step, n_en))
  
  return(Y)
}

#' observation error matrix, should be a square matrix where
#'   col & row = the number of states and params for which you have observations
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_param_obs number of parameters for which we have observations
#' @param n_step number of model timesteps
#' @param state_sd vector of state observation standard deviation; assuming sd is constant through time
#' @param param_sd vector of parmaeter observation standard deviation; assuming sd is constant through time
get_obs_error_matrix = function(n_states, n_params_obs, n_step, state_sd, param_sd){
  
  R = array(0, dim = c(n_states + n_params_obs, n_states + n_params_obs, n_step))
  
  state_var = state_sd^2 #variance of temperature observations
  
  param_var = param_sd^2
  
  if(n_params_obs > 0){
    all_var = c(state_var, param_var)
  }else{
    all_var = state_var
  }
  
  for(i in 1:n_step){
    # variance is the same for each depth and time step; could make dynamic or varying by time step if we have good reason to do so
    R[,,i] = diag(all_var, n_states + n_params_obs, n_states + n_params_obs)
  }
  
  return(R)
}

#' Measurement operator matrix saying 1 if there is observation data available, 0 otherwise
#'
#' @param n_states number of states we're updating in data assimilation routine
#' @param n_param_obs number of parameters for which we have observations
#' @param n_params_est number of parameters we're calibrating
#' @param n_step number of model timesteps
#' @param obs observation matrix created with get_obs_matrix function
get_obs_id_matrix = function(n_states, n_params_obs, n_params_est, n_step, obs){
  
  H = array(0, dim=c(n_states + n_params_obs, n_states + n_params_est, n_step))
  
  # order goes 1) states, 2)params for which we have obs, 3) params for which we're estimating but don't have obs
  
  for(t in 1:n_step){
    H[1:(n_states + n_params_obs), 1:(n_states + n_params_obs), t] = diag(ifelse(is.na(obs[,,t]),0, 1), n_states + n_params_obs, n_states + n_params_obs)
  }
  
  return(H)
}


#' turn observation dataframe into matrix
#'
#' @param obs_df observation data frame
#' @param model_dates dates over which you're modeling
#' @param n_step number of model time steps
#' @param n_states number of states we're updating in data assimilation routine
get_obs_matrix = function(obs_df, model_dates, n_step, n_states){
  
  # need to know location and time of observation
  
  obs_df_filtered = obs_df %>%
    dplyr::filter(as.Date(datetime) %in% model_dates) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, doc_lake) %>%
    mutate(date_step = which(model_dates %in% date))
  
  obs_matrix = array(NA, dim = c(n_states, 1, n_step))
  
  for(j in obs_df_filtered$date_step){
    obs_matrix[1, 1, j] = dplyr::filter(obs_df_filtered,
                                        date_step == j) %>%
      pull(doc_lake)
  }
  
  return(obs_matrix)
}



##' @param Y vector for holding states and parameters you're estimating
##' @param R observation error matrix
##' @param obs observations at current timestep
##' @param H observation identity matrix
##' @param n_en number of ensembles
##' @param cur_step current model timestep
kalman_filter = function(Y, R, obs, H, n_en, cur_step){
  
  cur_obs = obs[ , , cur_step]
  
  cur_obs = ifelse(is.na(cur_obs), 0, cur_obs) # setting NA's to zero so there is no 'error' when compared to estimated states
  
  ###### estimate the spread of your ensembles #####
  Y_mean = matrix(apply(Y[ , cur_step, ], MARGIN = 1, FUN = mean), nrow = length(Y[ , 1, 1])) # calculating the mean of each temp and parameter estimate
  delta_Y = Y[ , cur_step, ] - matrix(rep(Y_mean, n_en), nrow = length(Y[ , 1, 1])) # difference in ensemble state/parameter and mean of all ensemble states/parameters
  
  ###### estimate Kalman gain #########
  K = ((1 / (n_en - 1)) * delta_Y %*% t(delta_Y) %*% matrix(t(H[, , cur_step]))) %*%
    qr.solve(((1 / (n_en - 1)) * H[, , cur_step] %*% delta_Y %*% t(delta_Y) %*% matrix(t(H[, , cur_step])) + R[, , cur_step]))
  
  ###### update Y vector ######
  for(q in 1:n_en){
    Y[, cur_step, q] = Y[, cur_step, q] + K %*% (cur_obs - H[, , cur_step] %*% Y[, cur_step, q]) # adjusting each ensemble using kalman gain and observations
  }
  return(Y)
}



#' initialize Y vector with draws from distribution of obs
#'
#' @param Y Y vector
#' @param obs observation matrix
initialize_Y = function(Y, obs, init_params, n_states_est, n_params_est, n_params_obs, n_step, n_en, state_sd, param_sd){
  
  # initializing states with earliest observations and parameters
  first_obs = coalesce(!!!lapply(seq_len(dim(obs)[3]), function(i){obs[,,i]})) %>% # turning array into list, then using coalesce to find first obs in each position.
    ifelse(is.na(.), mean(., na.rm = T), .) # setting initial temp state to mean of earliest temp obs from other sites if no obs
  
  if(n_params_est > 0){
    ## update this later ***********************
    first_params = init_params
  }else{
    first_params = NULL
  }
  
  Y[ , 1, ] = array(abs(rnorm(n = n_en * (n_states_est + n_params_est),
                              mean = c(first_obs, first_params),
                              sd = c(state_sd, param_sd))),
                    dim = c(c(n_states_est + n_params_est), n_en))
  
  return(Y)
}


#' matrix for holding driver data
#'
#' @param drivers_df dataframe which holds all the driver data 
#' @param model_dates dates for model run 
#' @param n_drivers number of model drivers 
#' @param driver_colnames column names of the drivers in the driver dataframe 
#' @param driver_cv coefficient of variation for each driver data 
#' @param n_step number of model timesteps
#' @param n_en number of ensembles
get_drivers = function(drivers_df, model_dates, n_drivers, driver_colnames, driver_cv, n_step, n_en){
  
  drivers_filtered = drivers_df %>% 
    dplyr::filter(as.Date(datetime) %in% model_dates)
  
  drivers_out = array(NA, dim = c(n_step, n_drivers, n_en))
  
  for(i in 1:n_drivers){
    for(j in 1:n_step){
      drivers_out[j,i,] = rnorm(n = n_en, 
                                mean = as.numeric(drivers_filtered[j, driver_colnames[i]]),
                                sd = as.numeric(driver_cv[i] * drivers_filtered[j, driver_colnames[i]]))
    }
  }
  
  return(drivers_out) 
}


#' wrapper for running EnKF 
#' 
#' @param n_en number of model ensembles 
#' @param start start date of model run 
#' @param stop date of model run
#' @param time_step model time step, defaults to days 
#' @param obs_file observation file 
#' @param driver_file driver data file 
#' @param n_states_est number of states we're estimating 
#' @param n_params_est number of parameters we're estimating 
#' @param n_params_obs number of parameters for which we have observations 
#' @param decay_init initial decay rate of DOC 
#' @param obs_cv coefficient of variation of observations 
#' @param param_cv coefficient of variation of parameters 
#' @param driver_cv coefficient of variation of driver data for DOC Load, Discharge out, and Lake volume, respectively 
#' @param init_cond_cv initial condition CV (what we're )
EnKF = function(n_en = 100, 
                start,
                stop, 
                time_step = 'days', 
                obs_file = 'A_EcoForecast/Data/lake_c_data.rds',
                driver_file = 'A_EcoForecast/Data/lake_c_data.rds',
                n_states_est = 1, 
                n_params_est = 1,
                n_params_obs = 0, 
                decay_init = 0.005, 
                obs_cv = 0.1,
                param_cv = 0.25,
                driver_cv = c(0.2, 0.2, 0.2),
                init_cond_cv = 0.1){
  
  
  n_en = n_en
  start = as.Date(start)
  stop = as.Date(stop)
  time_step = 'days' 
  dates = get_model_dates(model_start = start, model_stop = stop, time_step = time_step)
  n_step = length(dates)
  
  # get observation matrix
  obs_df = readRDS(obs_file) %>% 
    select(datetime, doc_lake) 
  
  drivers_df = readRDS(driver_file) %>% 
    select(datetime, doc_load, water_out, lake_vol) 
  
  n_states_est = n_states_est # number of states we're estimating 
  
  n_params_est = n_params_est # number of parameters we're calibrating
  
  n_params_obs = n_params_obs # number of parameters for which we have observations
  
  decay_init = decay_init # Initial estimate of DOC decay rate day^-1 
  
  doc_init = obs_df$doc_lake[min(which(!is.na(obs_df$doc_lake)))]
  
  state_cv = obs_cv #coefficient of variation of DOC observations 
  state_sd = state_cv * doc_init
  init_cond_sd = init_cond_cv * doc_init
  
  param_cv = param_cv #coefficient of variation of DOC decay 
  param_sd = param_cv * decay_init
  
  # driver data coefficient of variation for DOC Load, Discharge out, and Lake volume, respectively 
  driver_cv = driver_cv 
  
  
  # setting up matrices
  # observations as matrix
  obs = get_obs_matrix(obs_df = obs_df,
                       model_dates = dates,
                       n_step = n_step,
                       n_states = n_states_est)
  
  # Y vector for storing state / param estimates and updates
  Y = get_Y_vector(n_states = n_states_est,
                   n_params_est = n_params_est,
                   n_step = n_step,
                   n_en = n_en)
  
  # observation error matrix
  R = get_obs_error_matrix(n_states = n_states_est,
                           n_params_obs = n_params_obs,
                           n_step = n_step,
                           state_sd = state_sd,
                           param_sd = param_sd)
  
  # observation identity matrix
  H = get_obs_id_matrix(n_states = n_states_est,
                        n_params_obs = n_params_obs,
                        n_params_est = n_params_est,
                        n_step = n_step,
                        obs = obs)
  
  # initialize Y vector
  Y = initialize_Y(Y = Y, obs = obs, init_params = decay_init, n_states_est = n_states_est,
                   n_params_est = n_params_est, n_params_obs = n_params_obs,
                   n_step = n_step, n_en = n_en, state_sd = init_cond_sd, param_sd = param_sd)
  
  # get driver data with uncertainty - dim = c(n_step, driver, n_en) 
  drivers = get_drivers(drivers_df = drivers_df, 
                        model_dates = dates,
                        n_drivers = 3, 
                        driver_colnames = c('doc_load', 'water_out', 'lake_vol'), 
                        driver_cv = driver_cv, 
                        n_step = n_step, 
                        n_en = n_en) 
  
  # start modeling
  for(t in 2:n_step){
    for(n in 1:n_en){
      
      # run model; 
      model_output = predict_lake_doc(doc_load = drivers[t-1, 1, n], 
                                      doc_lake = Y[1, t-1, n], 
                                      lake_vol = drivers[t-1, 3, n], 
                                      water_out = drivers[t-1, 2, n],
                                      decay = Y[2, t-1, n])
      
      Y[1 , t, n] = model_output$doc_predict # store in Y vector
      Y[2 , t, n] = model_output$decay
    }
    # check if there are any observations to assimilate 
    if(any(!is.na(obs[ , , t]))){
      Y = kalman_filter(Y = Y,
                        R = R,
                        obs = obs,
                        H = H,
                        n_en = n_en,
                        cur_step = t) # updating params / states if obs available
    }
  }
  out = list(Y = Y, dates = dates, drivers = drivers, R = R, obs = obs, state_sd = state_sd)
  
  return(out)
}