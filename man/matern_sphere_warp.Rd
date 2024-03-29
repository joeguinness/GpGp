% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{matern_sphere_warp}
\alias{matern_sphere_warp}
\alias{d_matern_sphere_warp}
\title{Deformed Matern covariance function on sphere}
\usage{
matern_sphere_warp(covparms, lonlat)

d_matern_sphere_warp(covparms, lonlat)
}
\arguments{
\item{covparms}{A vector with covariance parameters
in the form (variance, range, smoothness, nugget, <5 warping parameters>). 
Range parameter assumes that the sphere has radius 1 (units are radians).}

\item{lonlat}{A matrix with \code{n} rows and one column with longitudes in (-180,180)
and one column of latitudes in (-90,90).
Each row of lonlat describes a point on the sphere.}
}
\value{
A matrix with \code{n} rows and \code{n} columns, with the i,j entry
containing the covariance between observations at \code{lonlat[i,]} and
\code{lonlat[j,]}.
}
\description{
From a matrix of longitudes and latitudes and a vector covariance parameters of the form
(variance, range, smoothness, nugget, <5 warping parameters>), return the square matrix of
all pairwise covariances.
}
\section{Functions}{
\itemize{
\item \code{d_matern_sphere_warp()}: Derivatives with respect to parameters.

}}
\section{Warpings}{

The function first calculates the (x,y,z) 3D coordinates, and then "warps"
the locations to \eqn{(x,y,z) + \Phi(x,y,z)}, where \eqn{\Phi} is a warping
function composed of gradients of spherical harmonic functions of degree 2.
See Guinness (2019, "Gaussian Process Learning via Fisher Scoring of 
Vecchia's Approximation") for details.
The warped locations are input into \code{matern_isotropic}.
}

