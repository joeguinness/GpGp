
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
