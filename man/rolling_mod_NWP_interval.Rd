% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rolling_mod_NWP_interval.R
\name{rolling_mod_NWP_interval}
\alias{rolling_mod_NWP_interval}
\title{Rolling model}
\usage{
rolling_mod_NWP_interval(
  forc_start = as.Date("2007-01-01"),
  forc_end = as.Date("2023-01-01"),
  q = 0.9,
  ERA_NWP,
  model = "qreg",
  window = 125,
  reweight = FALSE,
  incl_climatology = FALSE,
  formula = "PC1 ~ NWP1_90",
  coef_to_print = c(),
  interval_k = 4,
  skill_interval = 0,
  cores = 4
)
}
\arguments{
\item{q}{Float. Quantile of interest.}

\item{ERA_NWP}{Data.table. Data of ERA observations and quantiles of NWP forecasts.}

\item{model}{String. Specify type of model, e.g. qreg, qgam, etc.}

\item{window}{Integer. Max number of forecast days from initialization.}

\item{reweight}{Boolean. Use reweighted data.}

\item{incl_climatology}{Boolean. Include climatology model.}

\item{formula}{String. Regression formula.}

\item{coef_to_print}{List of coefficients to include in output.}

\item{interval_k}{Integer Number of days in training interval}

\item{skill_interval}{}

\item{cores}{Integer. Number of cores to be used.}
}
\value{
Model output with predictions, loss and model coefficients.
}
\description{
Rolling model
}
\examples{
mod = rolling_mod_NWP_interval('2007-01-01', '2022-05-01', 0.9, ERA_NWP, 1, model = 'qreg', window = 60, reweight=FALSE)
}
