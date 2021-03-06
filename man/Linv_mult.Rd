% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Linv_mult}
\alias{Linv_mult}
\title{Multiply approximate inverse Cholesky by a vector}
\usage{
Linv_mult(Linv, z, NNarray)
}
\arguments{
\item{Linv}{Entries of the sparse inverse Cholesky factor,
usually the output from \code{\link{vecchia_Linv}}.}

\item{z}{the vector to be multiplied}

\item{NNarray}{A matrix of indices, usually the output from \code{\link{find_ordered_nn}}.
Row \code{i} contains the indices
of the observations that observation \code{i} conditions on. By convention,
the first element of row \code{i} is \code{i}.}
}
\value{
the product of the sparse inverse Cholesky factor with a vector
}
\description{
Vecchia's approximation implies a sparse approximation to the
inverse Cholesky factor of the covariance matrix. This function
returns the result of multiplying that matrix by a vector.
}
\examples{
n <- 2000
locs <- matrix( runif(2*n), n, 2 )
covparms <- c(2, 0.2, 0.75, 0.1)
ord <- order_maxmin(locs)
NNarray <- find_ordered_nn(locs,20)
Linv <- vecchia_Linv( covparms, "matern_isotropic", locs, NNarray )
z1 <- rnorm(n)
y <- fast_Gp_sim_Linv(Linv,NNarray,z1)
z2 <- Linv_mult(Linv, y, NNarray)
print( sum( (z1-z2)^2 ) )
}
