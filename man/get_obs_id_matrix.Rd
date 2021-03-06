% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{get_obs_id_matrix}
\alias{get_obs_id_matrix}
\title{Measurement operator matrix saying 1 if there is observation data available, 0 otherwise}
\usage{
get_obs_id_matrix(n_states, n_params_obs, n_params_est, n_step, obs)
}
\arguments{
\item{n_states}{number of states we're updating in data assimilation routine}

\item{n_params_est}{number of parameters we're calibrating}

\item{n_step}{number of model timesteps}

\item{obs}{observation matrix created with get_obs_matrix function}

\item{n_param_obs}{number of parameters for which we have observations}
}
\description{
Measurement operator matrix saying 1 if there is observation data available, 0 otherwise
}
