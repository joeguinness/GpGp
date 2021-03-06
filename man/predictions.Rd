% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictions.R
\name{predictions}
\alias{predictions}
\title{Compute Gaussian process predictions using Vecchia's approximations}
\usage{
predictions(
  fit = NULL,
  locs_pred,
  X_pred,
  y_obs = fit$y,
  locs_obs = fit$locs,
  X_obs = fit$X,
  beta = fit$betahat,
  covparms = fit$covparms,
  covfun_name = fit$covfun_name,
  m = 60,
  reorder = TRUE,
  st_scale = NULL
)
}
\arguments{
\item{fit}{GpGp_fit object, the result of \code{\link{fit_model}}}

\item{locs_pred}{prediction locations}

\item{X_pred}{Design matrix for predictions}

\item{y_obs}{Observations associated with locs_obs}

\item{locs_obs}{observation locations}

\item{X_obs}{Design matrix for observations}

\item{beta}{Linear mean parameters}

\item{covparms}{Covariance parameters}

\item{covfun_name}{Name of covariance function}

\item{m}{Number of nearest neighbors to use}

\item{reorder}{TRUE/FALSE for whether reordering should be done. This should
generally be kept at TRUE, unless testing out the effect of
reordering.}

\item{st_scale}{amount by which to scale the spatial and temporal
dimensions for the purpose of selecting neighbors. We recommend setting
this manually when using a spatial-temporal covariance function. When 
\code{lonlat = TRUE}, spatial scale is in radians (earth radius = 1).}
}
\description{
With the prediction locations ordered after the observation locations,
an approximation for the inverse Cholesky of the covariance matrix
is computed, and standard formulas are applied to obtain
the conditional expectation.
}
\details{
We can specify either a GpGp_fit object (the result of 
\code{\link{fit_model}}), OR manually enter the covariance function and
parameters, the observations, observation locations, and design matrix. We 
must specify the prediction locations and the prediction design matrix.
}
