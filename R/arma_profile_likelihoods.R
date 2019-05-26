
#' @export
arma_proflik_mean_variance <- function(subparms,covfun_name = "arma_matern_isotropic",
                    y,X,locs,NNarray,return_parms = FALSE){

    n <- length(y)
    covparms1 <- c(1,subparms)
    Linv <- arma_vecchia_Linv(covparms1,covfun_name,locs,NNarray)
    z <- Linv_mult(Linv,y,NNarray)
    B <- array(NA, dim(X))
    for(j in 1:ncol(X)){
        B[,j] <- Linv_mult(Linv,X[,j],NNarray)
    }
    # information matrix
    infomat <- crossprod(B)
    # is it numerically invertible?
    if( sum(is.na(infomat)) > 0 || max(infomat) == Inf || min( eigen(infomat)$values ) < 1e-9 ){
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
