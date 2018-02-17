# prediction functions



#' Compute Gaussian process predictions using Vecchia's approximations
#' 
#' With the prediction locations ordered after the observation locations,
#' an approximation for the inverse Cholesky of the covariance matrix
#' is computed, and standard formulas are applied to obtain
#' the conditional expectation.
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
#' @export
predictions <- function(covparms, covfun_name = "matern_isotropic", y_obs, 
    locs_obs, locs_pred, X_obs, X_pred, beta, m = 60, reorder = TRUE){
    
    n_obs <- nrow(locs_obs)
    n_pred <- nrow(locs_pred)
    
    # get orderings
    if(reorder){
        ord1 <- order_maxmin(locs_obs)
        ord2 <- order_maxmin(locs_pred)
    } else {
        ord1 <- 1:n_obs
        ord2 <- 1:n_pred
    }

    # reorder stuff
    Xord_obs  <- as.matrix( X_obs[ord1,] )
    Xord_pred <- as.matrix( X_pred[ord2,])
    yord_obs  <- y_obs[ord1]
    
    
    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,], locs_pred[ord2,] )
    inds1 <- 1:n_obs
    inds2 <- (n_obs+1):(n_obs+n_pred)
    
    # figure out if lonlat or not
    if( covfun_name == "matern_sphere" || covfun_name == "matern_sphere_time" ){
        lonlat <- TRUE
    } else {
        lonlat <- FALSE
    }
    
    # get nearest neighbor array (in space only)
    NNarray_all <- find_ordered_nn(locs_all,m=m,lonlat = lonlat)
    
    # get entries of Linv for obs locations and pred locations
    Linv_all <- vecchia_Linv(covparms,covfun_name,locs_all,NNarray_all)
    
    y_withzeros <- c(yord_obs - Xord_obs %*% beta, rep(0,n_pred) )
    v1 <- Linv_mult(Linv_all, y_withzeros, NNarray_all )
    v1[inds1] <- 0
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
#' @param nsims Number of conditional simulations to return.
#' @param reorder TRUE/FALSE for whether reordering should be done. This should
#' generally be kept at TRUE, unless testing out the effect of
#' reordering.
#' @export
cond_sim <- function(covparms, covfun_name = "matern_isotropic", y_obs, 
    locs_obs, locs_pred, X_obs, X_pred, beta, m = 60, nsims = 1, reorder = TRUE){
    
    n_obs <- nrow(locs_obs)
    n_pred <- nrow(locs_pred)
    
    # get orderings
    ord1 <- order_maxmin(locs_obs)
    ord2 <- order_maxmin(locs_pred)

    # reorder stuff
    Xord_obs  <- as.matrix( X_obs[ord1,] )
    Xord_pred <- as.matrix( X_pred[ord2,])
    yord_obs  <- y_obs[ord1]

    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,], locs_pred[ord2,] )
    inds1 <- 1:n_obs
    inds2 <- (n_obs+1):(n_obs+n_pred)
    
    # figure out if lonlat or not
    if( covfun_name == "matern_sphere" || covfun_name == "matern_sphere_time" ){
        lonlat <- TRUE
    } else {
        lonlat <- FALSE
    }
    
    # get nearest neighbor array (in space only)
    NNarray_all <- find_ordered_nn(locs_all,m=m,lonlat = lonlat)
    
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
