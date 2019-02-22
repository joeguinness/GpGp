# various profile likelioods

#' Profile likelihood (profiling out mean and variance)
#'
#' The profile likelihood is the maximum likelihood over a subset of
#' the parameters, given specified values of the remaining parameters.
#' In Gaussian process models, we can usually profile out linear mean
#' parameters and an overall variance (scale) parameter.
#'
#' @param subparms All parameters except for variance parameter. The specific meaning
#' of each parameter depends on \code{covfun_name}.
#' @param X design matrix, each column of X is a single covariate
#' @param return_parms flag for whether the function should return the loglikelihood 
#' only (\code{return_parms = FALSE}) or to return both the loglikelihood and 
#' all of the parameter values, including mean vector and variance parameter 
#' (\code{return_parms = TRUE}). Usually, we do the optimization using
#' \code{return_parms = FALSE} and then collect the parameter estimates
#' with another call with \code{return_parms = TRUE}.
#' @inheritParams vecchia_loglik
#' @details It is important that the ordering of \code{y} and \code{locs}
#' correspond to the ordering in \code{NNarray}. See example below.
#' @return Either the loglikelihood only (if \code{return_parms = FALSE}) or a list containing the loglikelihood, parameter values, and covariance matrix for linear mean parameters (if \code{return_parms = TRUE}).
#' @examples
#' n1 <- 50
#' n2 <- 50             # size of grid of locations
#' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
#' n <- nrow(locs)
#' covparms = c(3,0.1,1,0)    # variance, range, smoothness, nugget
#' X = as.matrix( rep(1,n) )  # design matrix
#'
#' # simulated response
#' y <- 2*X[,1] + fast_Gp_sim(covparms, "matern_isotropic", locs, m = 30)
#'
#' ord <- order_maxmin(locs)         # ordering of locations
#' yord <- y[ord]                    # reordered response
#' Xord <- as.matrix( X[ord,] )      # reordered design matrix
#' locsord <- locs[ord,]             # reordered locations
#' NNarray <- find_ordered_nn(locsord, m = 30)     # nearest neighbor indices
#'
#' # loglikelihood at true values of parameters
#' vecchia_loglik( covparms, "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNarray )
#' # profile out mean only (likelihood larger than vecchia_loglik)
#' proflik_mean( covparms[1:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNarray, return_parms = FALSE)
#' # profile out variance (likelihood larger than vecchia_loglik)
#' proflik_variance( covparms[2:4], "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNarray, return_parms = FALSE)
#' # profile out mean and variance (likelihood largest)
#' proflik_mean_variance( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNarray, return_parms = FALSE)
#' # get all parameter values 
#' proflik_mean_variance( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNarray, return_parms = TRUE)
#' 
#' 
#' @export
proflik_mean_variance <- function(subparms,covfun_name = "matern_isotropic",
                    y,X,locs,NNarray,return_parms = FALSE){

    n <- length(y)
    covparms1 <- c(1,subparms)
    Linv <- vecchia_Linv(covparms1,covfun_name,locs,NNarray)
    z <- Linv_mult(Linv,y,NNarray)
    B <- array(NA, dim(X))
    for(j in 1:ncol(X)){
        B[,j] <- Linv_mult(Linv,X[,j],NNarray)
    }
    # information matrix
    infomat <- crossprod(B)
    # is it numerically invertible?
    if( min( eigen(infomat)$values ) < 1e-9 ){
        invertible <- FALSE
    } else {
        invertible <- TRUE
    }
    if( invertible ){ 
        beta <- solve( infomat, crossprod(B,z) ) 
    } else {
        beta <- rep(0, ncol(X))
    }
    resids <- y - X %*% beta
    z_resids <- Linv_mult(Linv,resids,NNarray)
    sigmasq <- c( crossprod(z_resids)/n )

    logdet <- -2*sum(log(Linv[,1])) + n*log(sigmasq)
    quadform <- n
    if( invertible ){
        profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    } else {
        profloglik <- -1e6*stats::sd(y)*n
    }
    if( !return_parms ){
        return(profloglik)
    }
    if( return_parms ){
        betacovmat <- sigmasq*solve(infomat)
        return(list(loglik = profloglik, covparms = c(sigmasq,subparms),
                    beta = beta, betacovmat = betacovmat))
    }
}


#' Profile likelihood (profiling out variance in mean-zero model)
#'
#' The profile likelihood is the maximum likelihood over a subset of
#' the parameters, given specified values of the remaining parameters.
#' In Gaussian process models, we can usually profile out linear mean
#' parameters and an overall variance (scale) parameter.
#'
#' @inheritParams vecchia_loglik
#' @inheritParams proflik_mean_variance
#' @details It is important that the ordering of \code{y} and \code{locs}
#' correspond to the ordering in \code{NNarray}. See example below.
#' @return Either the loglikelihood only (if \code{return_parms = FALSE}) or a list containing the loglikelihood, parameter values, and covariance matrix for linear mean parameters (if \code{return_parms = TRUE}).
#' @examples
#' n1 <- 50
#' n2 <- 50             # size of grid of locations
#' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
#' n <- nrow(locs)
#' covparms = c(3,0.1,1,0)    # variance, range, smoothness, nugget
#' X = as.matrix( rep(1,n) )  # design matrix
#'
#' # simulated response
#' y <- 2*X[,1] + fast_Gp_sim(covparms, "matern_isotropic", locs, m = 30)
#'
#' ord <- order_maxmin(locs)         # ordering of locations
#' yord <- y[ord]                    # reordered response
#' Xord <- as.matrix( X[ord,] )      # reordered design matrix
#' locsord <- locs[ord,]             # reordered locations
#' NNarray <- find_ordered_nn(locsord, m = 30)     # nearest neighbor indices
#'
#' # loglikelihood at true values of parameters
#' vecchia_loglik( covparms, "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNarray )
#' # profile out mean only (likelihood larger than vecchia_loglik)
#' proflik_mean( covparms[1:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNarray, return_parms = FALSE)
#' # profile out variance (likelihood larger than vecchia_loglik)
#' proflik_variance( covparms[2:4], "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNarray, return_parms = FALSE)
#' # profile out mean and variance (likelihood largest)
#' proflik_mean_variance( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNarray, return_parms = FALSE)
#' # get all parameter values 
#' proflik_mean_variance( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNarray, return_parms = TRUE)
#' 
#' 
#' @export
proflik_variance <- function(subparms,covfun_name = "matern_isotropic",
                                  y,locs,NNarray,return_parms = FALSE){

    n <- length(y)
    covparms1 <- c(1,subparms)
    Linv <- vecchia_Linv(covparms1,covfun_name,locs,NNarray)
    z <- Linv_mult(Linv,y,NNarray)
    sigmasq <- c( crossprod(z)/n )

    logdet <- -2*sum(log(Linv[,1])) + n*log(sigmasq)
    quadform <- n
    profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    if( !return_parms ){
        return(profloglik)
    }
    if( return_parms ){
        return(list(loglik = profloglik, covparms = c(sigmasq,subparms) ) )
    }
}



#' Profile likelihood (profiling out mean only)
#'
#' The profile likelihood is the maximum likelihood over a subset of
#' the parameters, given specified values of the remaining parameters.
#' In Gaussian process models, we can usually profile out linear mean
#' parameters and an overall variance (scale) parameter.
#'
#' @param parms All parameters.  The specific meaning
#' of each parameter depends on \code{covfun_name}.
#' @inheritParams vecchia_loglik
#' @inheritParams proflik_mean_variance
#' @details It is important that the ordering of \code{y} and \code{locs}
#' correspond to the ordering in \code{NNarray}. See example below.
#' @return Either the loglikelihood only (if \code{return_parms = FALSE}) or a list containing the loglikelihood, parameter values, and covariance matrix for linear mean parameters (if \code{return_parms = TRUE}).
#' @examples
#' n1 <- 50
#' n2 <- 50             # size of grid of locations
#' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
#' n <- nrow(locs)
#' covparms = c(3,0.1,1,0)    # variance, range, smoothness, nugget
#' X = as.matrix( rep(1,n) )  # design matrix
#'
#' # simulated response
#' y <- 2*X[,1] + fast_Gp_sim(covparms, "matern_isotropic", locs, m = 30)
#'
#' ord <- order_maxmin(locs)         # ordering of locations
#' yord <- y[ord]                    # reordered response
#' Xord <- as.matrix( X[ord,] )      # reordered design matrix
#' locsord <- locs[ord,]             # reordered locations
#' NNarray <- find_ordered_nn(locsord, m = 30)     # nearest neighbor indices
#'
#' # loglikelihood at true values of parameters
#' vecchia_loglik( covparms, "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNarray )
#' # profile out mean only (likelihood larger than vecchia_loglik)
#' proflik_mean( covparms[1:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNarray, return_parms = FALSE)
#' # profile out variance (likelihood larger than vecchia_loglik)
#' proflik_variance( covparms[2:4], "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNarray, return_parms = FALSE)
#' # profile out mean and variance (likelihood largest)
#' proflik_mean_variance( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNarray, return_parms = FALSE)
#' # get all parameter values 
#' proflik_mean_variance( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNarray, return_parms = TRUE)
#' 
#' 
#' @export
proflik_mean <- function(parms,covfun_name = "matern_isotropic",
                                  y,X,locs,NNarray,return_parms = FALSE){

    n <- length(y)
    Linv <- vecchia_Linv(parms,covfun_name,locs,NNarray)
    z <- Linv_mult(Linv,y,NNarray)
    B <- array(NA, dim(X))
    for(j in 1:ncol(X)){
        B[,j] <- Linv_mult(Linv,X[,j],NNarray)
    }
    # information matrix
    infomat <- crossprod(B)
    # is it numerically invertible?
    if( min( eigen(infomat)$values ) < 1e-9 ){
        invertible <- FALSE
    } else {
        invertible <- TRUE
    }
    if( invertible ){ 
        beta <- solve( infomat, crossprod(B,z) ) 
    } else {
        beta <- rep(0, ncol(X))
    }
    resids <- y - X %*% beta
    z_resids <- Linv_mult(Linv,resids,NNarray)

    logdet <- -2*sum(log(Linv[,1]))
    quadform <- sum( z_resids^2 )
    if( invertible ){
        profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    } else {
        profloglik <- -1e6*stats::sd(y)*n
    }
    if( !return_parms ){
        return(profloglik)
    }
    if( return_parms ){
        betacovmat <- solve(infomat)
        return(list(loglik = profloglik, covparms = parms,
                    beta = beta, betacovmat = betacovmat))
    }
}






#' Grouped Version of Profile Likelihood (profiling out mean and variance)
#'
#' The profile likelihood is the maximum likelihood over a subset of
#' the parameters, given specified values of the remaining parameters.
#' In Gaussian process models, we can usually profile out linear mean
#' parameters and an overall variance (scale) parameter.
#'
#' @inheritParams vecchia_loglik
#' @inheritParams proflik_mean_variance
#' @param NNlist List object for grouped version of Vecchia's likelihood. Usually the result of \code{group_obs(NNarray)}.
#' @details It is important that the ordering of \code{y} and \code{locs}
#' correspond to the ordering in \code{NNarray}. See example below.
#' @return Either the loglikelihood only (if \code{return_parms = FALSE}) or a list containing the loglikelihood, parameter values, and covariance matrix for linear mean parameters (if \code{return_parms = TRUE}).
#' @examples
#' n1 <- 50
#' n2 <- 50             # size of grid of locations
#' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
#' n <- nrow(locs)
#' covparms = c(3,0.1,1,0)    # variance, range, smoothness, nugget
#' X = as.matrix( rep(1,n) )  # design matrix
#'
#' # simulated response
#' y <- 2*X[,1] + fast_Gp_sim(covparms, "matern_isotropic", locs, m = 30)
#'
#' ord <- order_maxmin(locs)         # ordering of locations
#' yord <- y[ord]                    # reordered response
#' Xord <- as.matrix( X[ord,] )      # reordered design matrix
#' locsord <- locs[ord,]             # reordered locations
#' NNarray <- find_ordered_nn(locsord, m = 30)     # nearest neighbor indices
#' NNlist <- group_obs(NNarray)
#'
#' # loglikelihood at true values of parameters
#' vecchia_loglik_grouped( covparms, "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNlist )
#' # profile out mean only (likelihood larger than vecchia_loglik)
#' proflik_mean_grouped( covparms[1:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNlist, return_parms = FALSE)
#' # profile out variance (likelihood larger than vecchia_loglik)
#' proflik_variance_grouped( covparms[2:4], "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNlist, return_parms = FALSE)
#' # profile out mean and variance (likelihood largest)
#' proflik_mean_variance_grouped( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNlist, return_parms = FALSE)
#' # get all parameter values 
#' proflik_mean_variance_grouped( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNlist, return_parms = TRUE)
#' 
#' 
#' @export
proflik_mean_variance_grouped <- function(subparms,covfun_name = "matern_isotropic",
                                  y,X,locs,NNlist,return_parms = FALSE){

    n <- length(y)
    covparms1 <- c(1,subparms)
    # get the entries of L inverse
    Linv <- vecchia_Linv_grouped(covparms1,covfun_name,locs,NNlist)
    # Linv y
    z <- Linv_mult_grouped(Linv,y,NNlist)
    # B = Linv X
    B <- array(NA, dim(X))
    for(j in 1:ncol(X)){
        B[,j] <- Linv_mult_grouped(Linv,X[,j],NNlist)
    }
    # information matrix
    infomat <- crossprod(B)
    # is it numerically invertible?
    if( min( eigen(infomat)$values ) < 1e-9 ){
        print("*")
        invertible <- FALSE
    } else {
        invertible <- TRUE
    }
    if( invertible ){ 
        beta <- solve( infomat, crossprod(B,z) ) 
    } else {
        beta <- rep(0, ncol(X))
    }
    resids <- y - X %*% beta
    z_resids <- Linv_mult_grouped(Linv,resids,NNlist)
    sigmasq <- c( crossprod(z_resids)/n )

    Linv_diag_indices <- cumsum( NNlist[["local_resp_inds"]] )
    logdet <- -2*sum(log(Linv[ Linv_diag_indices ])) + n*log(sigmasq)
    quadform <- n
    if( invertible ){
        profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    } else {
        profloglik <- -1e6*stats::sd(y)*n
    }
    if( !return_parms ){
        return(profloglik)
    }
    if( return_parms ){
        betacovmat <- sigmasq*solve(infomat)
        return(list(loglik = profloglik, covparms = c(sigmasq,subparms),
                    beta = beta, betacovmat = betacovmat))
    }
}


#' Grouped Version of Profile Likelihood (profiling out variance only)
#'
#' The profile likelihood is the maximum likelihood over a subset of
#' the parameters, given specified values of the remaining parameters.
#' In Gaussian process models, we can usually profile out linear mean
#' parameters and an overall variance (scale) parameter.
#'
#' @inheritParams vecchia_loglik
#' @inheritParams proflik_mean_variance
#' @inheritParams proflik_mean_variance_grouped
#' @details It is important that the ordering of \code{y} and \code{locs}
#' correspond to the ordering in \code{NNarray}. See example below.
#' @return Either the loglikelihood only (if \code{return_parms = FALSE}) or a list containing the loglikelihood, parameter values, and covariance matrix for linear mean parameters (if \code{return_parms = TRUE}).
#' @examples
#' n1 <- 50
#' n2 <- 50             # size of grid of locations
#' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
#' n <- nrow(locs)
#' covparms = c(3,0.1,1,0)    # variance, range, smoothness, nugget
#' X = as.matrix( rep(1,n) )  # design matrix
#'
#' # simulated response
#' y <- 2*X[,1] + fast_Gp_sim(covparms, "matern_isotropic", locs, m = 30)
#'
#' ord <- order_maxmin(locs)         # ordering of locations
#' yord <- y[ord]                    # reordered response
#' Xord <- as.matrix( X[ord,] )      # reordered design matrix
#' locsord <- locs[ord,]             # reordered locations
#' NNarray <- find_ordered_nn(locsord, m = 30)     # nearest neighbor indices
#' NNlist <- group_obs(NNarray)
#'
#' # loglikelihood at true values of parameters
#' vecchia_loglik_grouped( covparms, "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNlist )
#' # profile out mean only (likelihood larger than vecchia_loglik)
#' proflik_mean_grouped( covparms[1:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNlist, return_parms = FALSE)
#' # profile out variance (likelihood larger than vecchia_loglik)
#' proflik_variance_grouped( covparms[2:4], "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNlist, return_parms = FALSE)
#' # profile out mean and variance (likelihood largest)
#' proflik_mean_variance_grouped( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNlist, return_parms = FALSE)
#' # get all parameter values 
#' proflik_mean_variance_grouped( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNlist, return_parms = TRUE)
#' 
#' 
#' @export
proflik_variance_grouped <- function(subparms,covfun_name = "matern_isotropic",
                             y,locs,NNlist,return_parms = FALSE){

    n <- length(y)
    covparms1 <- c(1,subparms)
    Linv <- vecchia_Linv_grouped(covparms1,covfun_name,locs,NNlist)
    z <- Linv_mult_grouped(Linv,y,NNlist)
    sigmasq <- c( crossprod(z)/n )

    Linv_diag_indices <- cumsum( NNlist[["local_resp_inds"]] )
    logdet <- -2*sum(log(Linv[ Linv_diag_indices ])) + n*log(sigmasq)
    quadform <- n
    profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    if( !return_parms ){
        return(profloglik)
    }
    if( return_parms ){
        return(list(loglik = profloglik, covparms = c(sigmasq,subparms) ) )
    }
}



#' Grouped Version of Profile Likelihood (profiling out mean only)
#'
#' The profile likelihood is the maximum likelihood over a subset of
#' the parameters, given specified values of the remaining parameters.
#' In Gaussian process models, we can usually profile out linear mean
#' parameters and an overall variance (scale) parameter.
#'
#' @param parms All parameters. The specific meaning
#' of each parameter depends on \code{covfun_name}.
#' @inheritParams vecchia_loglik
#' @inheritParams proflik_mean_variance
#' @inheritParams proflik_mean_variance_grouped
#' @details It is important that the ordering of \code{y} and \code{locs}
#' correspond to the ordering in \code{NNarray}. See example below.
#' @return Either the loglikelihood only (if \code{return_parms = FALSE}) or a list containing the loglikelihood, parameter values, and covariance matrix for linear mean parameters (if \code{return_parms = TRUE}).
#' @examples
#' n1 <- 50
#' n2 <- 50             # size of grid of locations
#' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
#' n <- nrow(locs)
#' covparms = c(3,0.1,1,0)    # variance, range, smoothness, nugget
#' X = as.matrix( rep(1,n) )  # design matrix
#'
#' # simulated response
#' y <- 2*X[,1] + fast_Gp_sim(covparms, "matern_isotropic", locs, m = 30)
#'
#' ord <- order_maxmin(locs)         # ordering of locations
#' yord <- y[ord]                    # reordered response
#' Xord <- as.matrix( X[ord,] )      # reordered design matrix
#' locsord <- locs[ord,]             # reordered locations
#' NNarray <- find_ordered_nn(locsord, m = 30)     # nearest neighbor indices
#' NNlist <- group_obs(NNarray)
#'
#' # loglikelihood at true values of parameters
#' vecchia_loglik_grouped( covparms, "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNlist )
#' # profile out mean only (likelihood larger than vecchia_loglik)
#' proflik_mean_grouped( covparms[1:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNlist, return_parms = FALSE)
#' # profile out variance (likelihood larger than vecchia_loglik)
#' proflik_variance_grouped( covparms[2:4], "matern_isotropic", 
#'     yord - 2*Xord[,1], locsord, NNlist, return_parms = FALSE)
#' # profile out mean and variance (likelihood largest)
#' proflik_mean_variance_grouped( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNlist, return_parms = FALSE)
#' # get all parameter values 
#' proflik_mean_variance_grouped( covparms[2:4], "matern_isotropic", 
#'     yord, Xord, locsord, NNlist, return_parms = TRUE)
#' 
#' 
#' @export
proflik_mean_grouped <- function(parms,covfun_name = "matern_isotropic",
                         y,X,locs,NNlist,return_parms = FALSE){

    n <- length(y)
    Linv <- vecchia_Linv_grouped(parms,covfun_name,locs,NNlist)
    z <- Linv_mult_grouped(Linv,y,NNlist)
    B <- array(NA, dim(X))
    for(j in 1:ncol(X)){
        B[,j] <- Linv_mult_grouped(Linv,X[,j],NNlist)
    }
    # information matrix
    infomat <- crossprod(B)
    # is it numerically invertible?
    if( min( eigen(infomat)$values ) < 1e-9 ){
        invertible <- FALSE
    } else {
        invertible <- TRUE
    }
    if( invertible ){ 
        beta <- solve( infomat, crossprod(B,z) ) 
    } else {
        beta <- rep(0, ncol(X))
    }
    resids <- y - X %*% beta
    z_resids <- Linv_mult_grouped(Linv,resids,NNlist)

    Linv_diag_indices <- cumsum( NNlist[["local_resp_inds"]] )
    logdet <- -2*sum(log(Linv[ Linv_diag_indices ]))
    quadform <- sum( z_resids^2 )
    if( invertible ){
        profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    } else {
        profloglik <- -1e6*stats::sd(y)*n
    }
    if( !return_parms ){
        return(profloglik)
    }
    if( return_parms ){
        betacovmat <- solve(infomat)
        return(list(loglik = profloglik, covparms = parms,
                    beta = beta, betacovmat = betacovmat))
    }
}
