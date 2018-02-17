


#' Estimate mean and covariance parameters
#' 
#' Given a response, set of locations, (optionally) a design matrix, 
#' and a specified covariance function, return the maximum
#' approximate likelihood estimates, using Vecchia's 
#' likelihood approximation.
#' 
#' @param y response vector
#' @param locs matrix of locations. Each row is a single spatial or spatial-temporal
#' location. If using one of the "matern_sphere" covariance functions,
#' the locations should be longitudes and latitudes (in that order) in degrees.
#' @param X design matrix. Each row contains covariates for the corresponding
#' observation in \code{y}. If not specified, the function sets \code{X} to be a 
#' matrix with a single column of ones, that is, a constant mean function.
#' @param covfun_name string name of a covariance function. Currently supported
#' are "matern_isotropic", "matern_sphere", and "matern_sphere_time".
#' @param silent TRUE/FALSE for whether to print some information during fitting.
#' @param group TRUE/FALSE for whether to use the grouped version of 
#' the approximation (Guinness, 2018) or not.  The grouped version 
#' is used by default.
#' @param reorder TRUE/FALSE indicating whether maxmin ordering should be used
#' (TRUE) or whether no reordering should be done before fitting (FALSE).
#' @return A list object containing covariance parameter estimates,
#' mean parameter estimates, and covariance matrix for mean parameter estimates.
#' @details The \code{fit_model} is a user-friendly model fitting function
#' that automatically performs many of the auxiliary tasks needed for 
#' using Vecchia's approximation, including reordering, computing
#' nearest neighbors, grouping, and optimization. Optimization proceeds
#' in several steps, using increasingly accurate versions of the approximation. 
#' The first step using 5 neighbors, then 15 neighbors, and the last step
#' uses 30 neighbors. The actual number of neighbors in the grouped
#' version is guaranteed to be larger, though depends on the ordering and
#' the configuration of the locations. We recommend always using
#' \code{group = TRUE} since the grouping is guaranteed to improve
#' the approximation.
#'
#' The Jason-3 windspeed vignette is a useful source for a 
#' use-case of the \code{fit_model} function for data on sphere. The example below
#' shows a very small example with a simulated dataset in 2d.
#' 
#' @examples
#' n1 <- 10
#' n2 <- 10
#' n <- n1*n2
#' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
#' covparms <- c(2,0.1,1/2,0)
#' y <- 7 + fast_Gp_sim(covparms, "matern_isotropic", locs)
#' X <- as.matrix( rep(1,n) )
#' fit <- fit_model(y, locs, X, "matern_isotropic")
#' fit
#' 
#' 
#' @export
fit_model <- function(y, locs, X = NULL, covfun_name = "matern_isotropic",
    silent = FALSE, group = TRUE, reorder = TRUE){
    
    n <- length(y)

    # check that length of observation vector same as
    # number of locations
    if( nrow(locs) != n ){
        stop("length of observation vector y not equal
              to the number of locations (rows in locs)")
    }
    
    # check if design matrix is specified
    if( is.null(X) ){
        if(!silent) cat("Design matrix not specified, using constant mean \n")
        X <- rep(1,n) 
    }
    X <- as.matrix(X)
    
    # check if one of the allowed covariance functions is chosen
    if( ! covfun_name %in% c("matern_isotropic","matern_sphere",
                            "matern_sphere_time","matern_space_time") ) {
        stop("unrecognized covariance function name `covfun_name'. Choose from
             'matern_isotropic', 'matern_sphere', 'matern_sphere_time', or 'matern_space_time' " )
    }
    
    # missing values???
    
    # starting values and covariance-specific settings
    start_var <- stats::var(y)
    start_smooth <- 0.8
    start_nug <- 0.1
    
    randinds <- sample(1:n, min(n,200))    
    if(covfun_name == "matern_isotropic"){
        lonlat <- FALSE
        space_time <- FALSE
        dmat <- fields::rdist(locs[randinds,])
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
    }

    if(covfun_name == "matern_space_time"){
        lonlat <- FALSE
        space_time <- TRUE
        d <- ncol(locs)-1
        dmat1 <- fields::rdist(locs[randinds,1:d])
        dmat2 <- fields::rdist(locs[randinds,d+1,drop=FALSE])
        start_range1 <- mean( dmat1 )/4
        start_range2 <- mean( dmat2 )/1
        start_parms <- c(start_var, start_range1, start_range2, start_smooth, start_nug)
    }

    if( covfun_name == "matern_sphere" ){
        lonlat <- TRUE
        space_time <- FALSE
        dmat <- fields::rdist.earth(locs[randinds,], R = 1)
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
        if(!silent) cat("Using 'matern_sphere'. Assuming that argument 'locs' is (longitude,latitude)\n")    
    }
    
    if( covfun_name == "matern_sphere_time" ){
        lonlat <- TRUE
        space_time <- TRUE
        dmat <- fields::rdist.earth(locs[randinds,1:2], R = 1)
        start_range1 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,3])
        start_range2 <- mean( dmat )/4
        start_parms <- c(start_var, start_range1, start_range2, start_smooth, start_nug)
        if(!silent) cat("Using 'matern_sphere_time'. Assuming that argument 'locs' is (longitude,latitude,time)\n")    
    }
    
    nparms <- length(start_parms)
    
    # get an ordering and reorder everything
    if(reorder){
        if(!silent) cat("Reordering...")
        ord <- order_maxmin(locs, lonlat = lonlat, space_time = space_time)
        if(!silent) cat("Done \n")
    } else {
        ord <- 1:n
    }
    yord <- y[ord]
    locsord <- locs[ord,]
    Xord <- as.matrix( X[ord,] )

    # get nearest neighbors    
    if(!silent) cat("Finding nearest neighbors...")
    NNarray <- find_ordered_nn(locsord, m=30, lonlat = lonlat, space_time = space_time)
    if(!silent) cat("Done \n")
    
    fit <- list(par=log(start_parms[2:nparms]))

    # refine the estimates for m = c(15,30,45)
    for(m in c(5,15,30)){
        if(space_time){
            NNarray <- find_ordered_nn(locsord, m=m, lonlat = lonlat, 
                space_time = space_time, st_scale = exp(fit$par[1:2]) )
        }
        NNlist <- group_obs(NNarray[,1:(m+1)])
        funtomax <- function( logparms ){
            parms <- exp(logparms)
            # set a maximum value for smoothness
            parms[length(parms)-1] = min(parms[length(parms)-1], 4)
            loglik <- proflik_mean_variance_grouped(
                parms, covfun_name, yord, Xord, locsord, NNlist )
            return(-loglik)
        }
        if(!silent) cat("Refining estimates...")
        fit <- stats::optim(fit$par,funtomax,
            control=list(trace=0,maxit=25), method = "BFGS") 
        if(!silent) cat("Done          ")
        if(!silent) cat(paste(  paste( round(exp(fit$par),3), collapse = " " ), "\n" ) )
    }
    fitted_model <- proflik_mean_variance_grouped(
            exp(fit$par), covfun_name, yord, Xord, locsord, NNlist, return_parms = TRUE )
        
    return(fitted_model)
}