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


#' Extract time and space information
#'
#' Extracts time (from date column) and space (aka depth) information
#'
#' @param t matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param y matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param parns matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return list of datetimes and depths
#' @export
run_model <- function(modelfunc = 'TwoLayer', bc, params, ini, times){
  
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
  NEP <- params[18]
  Fsed<- params[19]
  MINERAL <- params[20]
  Ased <- params[21]
  
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
  
  TwoLayerOxy <- function(t, y, parms){
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
      dV = 1e2 #100
      mult = 1/100#1/1000
    } else {
      dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000)
      mult = 1.0
    }
    
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
    
    K600 <- k.cole.base(2)
    water.tmp = y[1]
    kO2 <- k600.2.kGAS.base(k600=K600, 
                                 temperature=water.tmp, 
                                 gas='O2') # m/d
    o2sat<-o2.at.sat.base(temp= water.tmp, 
                               altitude = 300)*1000 # mg O2/L
    Fatm <- kO2*(o2sat * Ve - y[3] ) * (As/Ve) # mg * m/d * m2/m3= mg/d
    #  kO2*(o2sat - y[3] )/ (As/Ve) # mg/m m/d = mg/d
    
    Sed <- Fsed * y[4] * (Ased/Vh)  * 1.03^(y[2]-20) #* mult # m/d * mg * m2/m3
    #  Fsed * y[4]/ (Ased/Vh)  * 1.08^(y[2]-20) * mult
    
    PP <- 1.08^(y[1]-20) * NEP * Ve * mult # mg/m3/d * m3
    # 1.08^(y[1]-20) * NEP * Ve * mult 
    
    VOL <- 1.08^(y[2]-20) * MINERAL * Vh * mult # mg/m3/d * m3
    
    dOe <-( PP +
      Fatm +
      ((dV * At) / Ve) * (y[4] - y[3]) ) 
    
    dOh <- ( ((dV * At) / Vh) * (y[3] - y[4]) - 
               - VOL - 
      Sed) 
    
    
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
    
    return(list(c(dTe, dTh, dOe, dOh)))
  }
  
  # approximating all boundary conditions 
  Jsw <- approxfun(x = bc$Month, y = bc$Jsw, method = "linear", rule = 2)
  Tair <- approxfun(x = bc$Month, y = bc$Tair, method = "linear", rule = 2)
  Dew <- approxfun(x = bc$Month, y = bc$Dew, method = "linear", rule = 2)
  Uw <- approxfun(x = bc$Month, y = bc$Uw, method = "linear", rule = 2)
  vW <- approxfun(x = bc$Month, y = bc$vW, method = "linear", rule = 2)
  vt <- approxfun(x = bc$Month, y = bc$vt, method = "linear", rule = 2)

  # runge-kutta 4th order
  if (modelfunc == 'TwoLayer'){
    out <- ode(times = times, y = ini, func = TwoLayer, parms = params, method = 'rk4')
  } else if (modelfunc == 'TwoLayerOxy'){
    out <- ode(times = times, y = ini, func = TwoLayerOxy, parms = params, method = 'rk4')
  }
 
  
  return(out)
}
