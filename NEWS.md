
# GpGp 0.5.1

Fixed minor issue in documentation.

# GpGp 0.5.0

Bug fixes related to automatic coercion of matrices to vectors
One-dimensional inputs and single-location predictions work now
Workaround for change to roxygen that removed package-alias
Removed maptools from Suggests field of DESCRIPTION

# GpGp 0.4.0

Speedups of evaluations of the likelihood function.  Modifications to
the Fisher scoring algorithm to make it a little more stable.

# GpGp 0.3.2

The major purpose of this release is to fix an error
on CRAN checks for Oracle Developer Studio. We have also
parallelized the computation of the inverse cholesky
factor, and included safeguards against large smoothness
parameters in matern_isotropic.

Many thanks for Youssef Fahmy for finally cracking
the solaris error.

# GpGp 0.3.1

Added matern_anisotropic3D_alt covariance
Fixed some problems with #includes in src files

# GpGp 0.3.0

Now uses OpenMP for parallel computations of the likelihood.

Updated behavior of Fisher scoring algorithm when information
matrix is ill-conditioned (simple regularization of info matrix).


# GpGp 0.2.2

Fixed bug in "fit_model" when missing values are present.
Updated behavior of Fisher scoring algorithm when information matrix ill-conditioned
Allow user to fix a subset of parameters in "fit_model"
Allow user to specify maximum number of iterations in "fit_model"
New faster computational algorithm for predictions
Several new covariance functions, including

  matern15_isotropic
  matern25_isotropic
  matern35_isotropic
  matern45_isotropic
  matern15_scaledim  
  matern25_scaledim  
  matern35_scaledim  
  matern45_scaledim
  

# GpGp 0.2.1

Bug fix for overloaded use of 'pow' function in 'basis.h'

# GpGp 0.2.0

This update includes an implementation of the Fisher Scoring
algorithm described in this paper <https://arxiv.org/abs/1905.08374>,
computed in a single pass through the data.

Much of the C++ code has been rewritten and reorganized,
making use of the Armadillo C++ linear algebra library,
with the help of RcppArmadillo.

There are also several new covariance functions. The complete list of
covariance functions is now:

matern_isotropic
exponential_isotropic
matern_spacetime
exponential_spacetime
matern_scaledim
exponential_scaledim
matern_anisotropic2D
exponential_anisotropic2D
exponential_anisotropic3D
matern_nonstat_var
exponential_nonstat_var
matern_sphere
exponential_sphere
matern_spheretime
exponential_spheretime
matern_sphere_warp
exponential_sphere_warp
matern_spheretime_warp
exponential_spheretime_warp


# GpGp 0.1.1

This is a minor release fixing numerical stability problems
that arise during optimization of the likelihood.

* Added a check in each of the proflik_mean* functions
  to avoid inverting the information matrix when it
  is numerically singular.

* Changed the default number of Nelder-Mead iterations
  in fit_model to 100
