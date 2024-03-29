% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{exponential_anisotropic3D}
\alias{exponential_anisotropic3D}
\alias{d_exponential_anisotropic3D}
\title{Geometrically anisotropic exponential covariance function (three dimensions)}
\usage{
exponential_anisotropic3D(covparms, locs)

d_exponential_anisotropic3D(covparms, locs)
}
\arguments{
\item{covparms}{A vector with covariance parameters
in the form (variance, L11, L21, L22, L31, L32, L33, nugget)}

\item{locs}{A matrix with \code{n} rows and \code{3} columns.
Each row of locs is a point in R^3.}
}
\value{
A matrix with \code{n} rows and \code{n} columns, with the i,j entry
containing the covariance between observations at \code{locs[i,]} and
\code{locs[j,]}.
}
\description{
From a matrix of locations and covariance parameters of the form
(variance, L11, L21, L22, L31, L32, L33, nugget), return the square matrix of
all pairwise covariances.
}
\section{Functions}{
\itemize{
\item \code{d_exponential_anisotropic3D()}: Derivatives of anisotropic exponential covariance

}}
\section{Parameterization}{

The covariance parameter vector is (variance, L11, L21, L22, L31, L32, L33, nugget)
where L11, L21, L22, L31, L32, L33 are the six non-zero entries of a lower-triangular
matrix L. The covariances are 
\deqn{ M(x,y) = \sigma^2 exp(-|| L x - L y || ) }
This means that L11 is interpreted as an inverse range parameter in the
first dimension.
The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }.
}

