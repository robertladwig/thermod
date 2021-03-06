% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{run_temp_oxygen_forecast}
\alias{run_temp_oxygen_forecast}
\title{Simple one-layer water temperature and oxygen model in forecast mode}
\usage{
run_temp_oxygen_forecast(bc, params, ini, times, ice = FALSE, observed)
}
\arguments{
\item{bc}{meteorological boundary conditions: day, shortwave radiation, air temperature, dew
point temperature, wind speed and wind shear stress}

\item{params}{configuration parameters}

\item{ini}{vector of the initial water temperatures of the epilimnion and hypolimnion}

\item{times}{vector of time information}

\item{ice}{boolean, if TRUE the model will approximate ice on/off set, otherwise no ice}
}
\value{
matrix of simulated water temperature and oxygen mass in the epilimnion and hypolimnion
}
\description{
Runs a simple one-layer water temperature and oxygen model in dependence of meteorological 
drivers and lake setup/configuration.
}
