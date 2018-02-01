# prediction functions



#' Compute Gaussian process predictions using Vecchia's approximations
#' 
#' The prediction locations after the observation locations,
#' an approximation for the inverse Cholesky of the covariance matrix
#' is computed, and standard formulas are applied to obtain
#' the conditional expectation.
#' @param covparms Covariance parameters
#' @param covfun_name Name of covariance function
#' @param locs_obs observation locations
#' @param locs_pred prediction locations
#' @param X_obs Design matrix for observations
#' @param X_pred Design matrix for predictions
#' @param beta Linear mean parameters
#' @param y_obs Observations assocaited with locs_obs
#' @param m Number of nearest neighbors to use
#' @param lonlat Flag indicating that locs contain longitudes and latitudes
#' @export
predictions <- function(covparms, covfun_name = "matern_isotropic", locs_obs, locs_pred, 
    X_obs, X_pred, beta, y_obs, m = 30, lonlat = FALSE){
    
    n_obs <- nrow(locs_obs)
    n_pred <- nrow(locs_pred)

    # put all coordinates together
    locs_all <- rbind( locs_obs, locs_pred )
    inds1 <- 1:n_obs
    inds2 <- (n_obs+1):(n_obs+n_pred)
    
    # get nearest neighbor array (in space only)
    NNarray_all <- find_ordered_nn(locs_all,m=m,lonlat = lonlat)
    
    # get entries of Linv for obs locations and pred locations
    Linv_all <- vecchia_Linv(covparms,covfun_name,locs_all,NNarray_all)
    
    y_withzeros <- c(y_obs - X_obs %*% beta, rep(0,n_pred) )
    v1 <- Linv_mult(Linv_all, y_withzeros, NNarray_all )
    v1[inds1] <- 0
    v2 <- -L_mult(Linv_all,v1,NNarray_all)

    condexp <- c(v2[inds2] + X_pred %*% beta)
    return(condexp)
}
