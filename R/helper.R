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


#' Extract time and space information
#'
#' Extracts time (from date column) and space (aka depth) information
#'
#' @param t matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param y matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param parns matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return list of datetimes and depths
#' @export
TwoLayer <- function(t, y, parms){
  eair <- (4.596 * exp((17.27 * Dew(t)) / (237.3 + Dew(t)))) # air vapor pressure
  esat <- 4.596 * exp((17.27 * Tair(t)) / (237.3 + Tair(t))) # saturation vapor pressure
  RH <- eair/esat *100 # relative humidity
  es <- 4.596 * exp((17.27 * y[1])/ (273.3+y[1]))
  # diffusion coefficient
  Cd <- 0.00052 * (Uw(t))^(0.44)
  shear <- 1.164/1000 * Cd * (Uw(t))^2
  w0 <- sqrt(shear/rho)
  E0  <- c * w0
  
  rho_e <- calc_dens(y[1])/1000
  rho_h <- calc_dens(y[2])/1000
  Ri <- ((g/rho)*(abs(rho_e-rho_h)/10))/(w0/(10)^2)
  dV <- (E0 / (1 + a * Ri)^(3/2))/(Ht/100) * (86400/10000)
  
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