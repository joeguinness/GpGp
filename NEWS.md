# GpGp 0.1.1

This is a minor release fixing numerical stability problems 
that arise during optimization of the likelihood.

* Added a check in each of the proflik_mean* functions 
  to avoid inverting the information matrix when it
  is numerically singular.
  
* Changed the default number of Nelder-Mead iterations 
  in fit_model to 100