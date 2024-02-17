# prediction functions



#' Compute Gaussian process predictions using Vecchia's approximations
#' 
#' With the prediction locations ordered after the observation locations,
#' an approximation for the inverse Cholesky of the covariance matrix
#' is computed, and standard formulas are applied to obtain
#' the conditional expectation.
#' @param fit GpGp_fit object, the result of \code{\link{fit_model}}
#' @param covparms Covariance parameters
#' @param covfun_name Name of covariance function
#' @param y_obs Observations associated with locs_obs
#' @param locs_obs observation locations
#' @param locs_pred prediction locations
#' @param X_obs Design matrix for observations
#' @param X_pred Design matrix for predictions
#' @param beta Linear mean parameters
#' @param m Number of nearest neighbors to use
#' @param reorder TRUE/FALSE for whether reordering should be done. This should
#' generally be kept at TRUE, unless testing out the effect of
#' reordering.
#' @param st_scale amount by which to scale the spatial and temporal
#' dimensions for the purpose of selecting neighbors. We recommend setting
#' this manually when using a spatial-temporal covariance function. When 
#' \code{lonlat = TRUE}, spatial scale is in radians (earth radius = 1).
#' @details We can specify either a GpGp_fit object (the result of 
#' \code{\link{fit_model}}), OR manually enter the covariance function and
#' parameters, the observations, observation locations, and design matrix. We 
#' must specify the prediction locations and the prediction design matrix.
#' @export
predictions <- function(fit = NULL, locs_pred, X_pred, 
    y_obs = fit$y, locs_obs = fit$locs, X_obs = fit$X, beta = fit$betahat,    
    covparms = fit$covparms, covfun_name = fit$covfun_name, 
    m = 60, reorder = TRUE, st_scale = NULL){
    
    # make sure things are matrices
    locs_obs <- as.matrix( locs_obs )
    X_obs <- as.matrix( X_obs )

    locs_pred <- as.matrix(locs_pred)
    X_pred <- as.matrix(X_pred)

    # how many observations and predictions 
    n_obs <- nrow(locs_obs)
    n_pred <- nrow(locs_pred)
    
    # get orderings
    if(reorder){

        if( n_obs < 6e4 ){
            ord1 <- order_maxmin(locs_obs)
        } else {
            ord1 <- sample( 1:n_obs )
        }
        
        if( n_pred < 6e4 ){
            ord2 <- order_maxmin(locs_pred)
        } else {
            ord2 <- sample( 1:n_pred )
        }

    } else {
        ord1 <- 1:n_obs
        ord2 <- 1:n_pred
    }

    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    yord_obs  <- y_obs[ord1]
    
    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )
    inds1 <- 1:n_obs
    inds2 <- (n_obs+1):(n_obs+n_pred)
    
    # figure out if lonlat or not
    lonlat <- get_linkfun(covfun_name)$lonlat
    space_time <- get_linkfun(covfun_name)$space_time
    
    # get nearest neighbor array (in space only)
    NNarray_all <- find_ordered_nn(locs_all,m=m, lonlat = lonlat,st_scale=st_scale)
    
    # get entries of Linv for obs locations and pred locations
    Linv_all <- vecchia_Linv(covparms,covfun_name,locs_all, NNarray_all, n_obs+1)
    
    y_withzeros <- c(yord_obs - c(Xord_obs %*% beta), rep(0,n_pred) )
    v1 <- Linv_mult(Linv_all, y_withzeros, NNarray_all )
    v1[inds1] <- 0
    Linv_all[1:n_obs,1] <- 1.0
    v2 <- -L_mult(Linv_all,v1,NNarray_all)

    condexp <- c(v2[inds2] + Xord_pred %*% beta)
    condexp[ord2] <- condexp
    return(condexp)
}



#' Conditional Simulation using Vecchia's approximation
#' 
#' With the prediction locations ordered after the observation locations,
#' an approximation for the inverse Cholesky of the covariance matrix
#' is computed, and standard formulas are applied to obtain
#' a conditional simulation.
#' @param fit GpGp_fit object, the result of \code{\link{fit_model}}
#' @param covparms Covariance parameters
#' @param covfun_name Name of covariance function
#' @param y_obs Observations associated with locs_obs
#' @param locs_obs observation locations
#' @param locs_pred prediction locations
#' @param X_obs Design matrix for observations
#' @param X_pred Design matrix for predictions
#' @param beta Linear mean parameters
#' @param m Number of nearest neighbors to use. Larger \code{m} gives
#' better approximations.
#' @param st_scale amount by which to scale the spatial and temporal
#' dimensions for the purpose of selecting neighbors. We recommend setting
#' this manually when using a spatial-temporal covariance function. When 
#' \code{lonlat = TRUE}, spatial scale is in radians (earth radius = 1).
#' @param nsims Number of conditional simulations to return.
#' @param reorder TRUE/FALSE for whether reordering should be done. This should
#' generally be kept at TRUE, unless testing out the effect of
#' reordering.
#' @details We can specify either a GpGp_fit object (the result of 
#' \code{\link{fit_model}}), OR manually enter the covariance function and
#' parameters, the observations, observation locations, and design matrix. We 
#' must specify the prediction locations and the prediction design matrix.
#' @export
cond_sim <- function(fit = NULL, locs_pred, X_pred, 
    y_obs = fit$y, locs_obs = fit$locs, X_obs = fit$X, beta = fit$betahat,    
    covparms = fit$covparms, covfun_name = fit$covfun_name, 
    m = 60, reorder = TRUE, st_scale = NULL, nsims = 1 ){

    # make sure things are matrices
    locs_obs <- as.matrix( locs_obs )
    X_obs <- as.matrix( X_obs )

    locs_pred <- as.matrix(locs_pred)
    X_pred <- as.matrix(X_pred)

    # how many observations and predictions 
    n_obs <- nrow(locs_obs)
    n_pred <- nrow(locs_pred)
    
    # get orderings
    if(reorder){

        if( n_obs < 6e4 ){
            ord1 <- order_maxmin(locs_obs)
        } else {
            ord1 <- sample( 1:n_obs )
        }
        
        if( n_pred < 6e4 ){
            ord2 <- order_maxmin(locs_pred)
        } else {
            ord2 <- sample( 1:n_pred )
        }

    } else {
        ord1 <- 1:n_obs
        ord2 <- 1:n_pred
    }

    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    yord_obs  <- y_obs[ord1]

    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )
    inds1 <- 1:n_obs
    inds2 <- (n_obs+1):(n_obs+n_pred)
    
    # figure out if lonlat or not
    lonlat <- get_linkfun(covfun_name)$lonlat
    space_time <- get_linkfun(covfun_name)$space_time
    
    # get nearest neighbor array (in space only)
    NNarray_all <- find_ordered_nn(locs_all,m=m,lonlat = lonlat,st_scale=st_scale)

    # get entries of Linv for obs locations and pred locations
    Linv_all <- vecchia_Linv(covparms,covfun_name,locs_all,NNarray_all)
    
    # an unconditional simulation
    condsim <- matrix(NA, n_pred, nsims)
    for(j in 1:nsims){
        z <- L_mult(Linv_all, stats::rnorm(n_obs+n_pred), NNarray_all)
    
        y_withzeros <- c(yord_obs - Xord_obs %*% beta + z[inds1], rep(0,n_pred) )
        v1 <- Linv_mult(Linv_all, y_withzeros, NNarray_all )
        v1[inds1] <- 0
        v2 <- -L_mult(Linv_all,v1,NNarray_all)

        condsim[ord2,j] <- c(v2[inds2] + Xord_pred %*% beta) - z[inds2]
    }
    return(condsim)
}
