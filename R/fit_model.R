


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
fit_model_optim <- function(y, locs, X = NULL, NNarray = NULL,
    covfun_name = "arma_matern_isotropic",
    silent = FALSE, group = TRUE, reorder = TRUE, method = "Nelder-Mead"){

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
    if( ! covfun_name %in% c(
        "arma_exponential_isotropic",
        "arma_matern_isotropic",
        "arma_matern_anisotropic2D",
        "arma_matern_anisotropic3D",
        "arma_matern_nonstat_var",
        "arma_matern_sphere",
        "arma_matern_sphere_warp",
        "arma_matern_spheretime_warp",
        "arma_matern_spheretime",
        "arma_matern_spacetime"   )
    ){
        stop("unrecognized covariance function name `covfun_name'." )
    }

    # missing values???

    # starting values and covariance-specific settings
    start <- get_start_parms_linkfun(y,X,locs,covfun_name)
    start_parms <- start$start_parms
    link <- start$link
    invlink <- start$invlink
    lonlat <- start$lonlat
    space_time <- start$space_time

    nparms <- length(start_parms)

    # get an ordering and reorder everything
    if(reorder){
        if(!silent) cat("Reordering...")
        if( n < 1e5 ){  # maximum ordering if n < 100000
            ord <- order_maxmin(locs, lonlat = lonlat, space_time = space_time)
        } else {        # otherwise random order
            ord <- sample(n)
        }
        if(!silent) cat("Done \n")
    } else {
        ord <- 1:n
    }
    yord <- y[ord]
    locsord <- locs[ord,]
    Xord <- as.matrix( X[ord,] )

    # get neighbor array if not provided
    if( is.null(NNarray) ){
        if(!silent) cat("Finding nearest neighbors...")
        NNarray <- find_ordered_nn(locsord, m=30, lonlat = lonlat, space_time = space_time)
        if(!silent) cat("Done \n")
    }

    fit <- list(par=invlink(start_parms)[2:nparms])

    # refine the estimates for m = c(15,30,45)
    for(m in c(10,30)){
        print(exp(fit$par[1:2]))
        if(group){ NNlist <- group_obs(NNarray[,1:(m+1)]) }
        funtomax <- function( logparms ){
            parms <- link(c(0,logparms))[2:nparms]
            if(group){
                stop("grouped version not implemented yet")
                #loglik <- arma_proflik_mean_variance_grouped(
                #    parms, covfun_name, yord, Xord, locsord, NNlist )
            } else {
                loglik <- arma_proflik_mean_variance( parms, covfun_name,
                    yord, Xord, locsord, NNarray[,1:(m+1)] )
            }
            return(-loglik)
        }
        if(!silent) cat("Refining estimates...")
        fit <- stats::optim(fit$par,funtomax,
            control=list(trace=5,maxit=10000), method = method)
        # set smoothness to maximum smoothness
        parms <- link(c(0,fit$par))
        # need to be more careful about this part
        #parms[length(parms)-1] <- min(parms[length(parms)-1], max_smooth)
        fit$par <- invlink(parms)[2:nparms]

        if(!silent) cat("Done          ")
        if(!silent) cat(paste(  paste( round(link(c(0,fit$par))[2:nparms],3), collapse = " " ), "\n" ) )
    }
    if(group){
        stop("grouped not implemented yet")
        #fitted_model <- arma_proflik_mean_variance_grouped(
        #    link(fit$par), covfun_name, yord, Xord, locsord, NNlist, return_parms = TRUE )
    } else {
        fitted_model <- arma_proflik_mean_variance(
            link(c(0,fit$par))[2:nparms], covfun_name, yord, Xord, locsord, NNarray, return_parms = TRUE )
    }
    fitted_model$loglik <- fitted_model$loglik
    return(fitted_model)
}













#' @export
fit_model <- function(y, locs, X = NULL, NNarray = NULL,
    covfun_name = "matern_isotropic", start_parms = NULL,
    silent = FALSE, group = TRUE, reorder = TRUE,
    m_seq = c(10,30)){

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
    if( ! covfun_name %in%
            c("arma_exponential_isotropic",
              "arma_matern_isotropic",
              "arma_matern_anisotropic2D",
              "arma_exponential_anisotropic3D",
              "arma_matern_anisotropic3D",
              "arma_matern_nonstat_var",
              "arma_matern_sphere",
              "arma_matern_sphere_warp",
              "arma_matern_spheretime_warp",
              "arma_matern_spheretime",
              "arma_matern_spacetime"   ) )
    {
        stop("unrecognized covariance function name `covfun_name'.")
    }

    # missing values???

    # get starting values for parameters
    if(is.null(start_parms)){
        start <- get_start_parms(y,X,locs,covfun_name)
        start_parms <- start$start_parms
    }
    linkfuns <- get_linkfun(y,X,locs,covfun_name)
    link <- linkfuns$link
    dlink <- linkfuns$dlink
    invlink <- linkfuns$invlink
    lonlat <- linkfuns$lonlat
    space_time <- linkfuns$space_time
    penalty <- get_penalty(y,X,locs,covfun_name)
    pen <- penalty$pen
    dpen <- penalty$dpen
    ddpen <- penalty$ddpen

    # get an ordering and reorder everything
    if(reorder){
        if(!silent) cat("Reordering...")
        if( n < 1e5 ){  # maximum ordering if n < 100000
            ord <- order_maxmin(locs, lonlat = lonlat, space_time = space_time)
        } else {        # otherwise random order
            ord <- sample(n)
        }
        if(!silent) cat("Done \n")
    } else {
        ord <- 1:n
    }
    yord <- y[ord]
    locsord <- locs[ord,]
    Xord <- as.matrix( X[ord,] )

    # get neighbor array if not provided
    if( is.null(NNarray) ){
        if(!silent) cat("Finding nearest neighbors...")
        NNarray <- find_ordered_nn(locsord, m=max(m_seq), lonlat = lonlat, space_time = space_time)
        if(!silent) cat("Done \n")
    }

    # refine the estimates for m in m_seq
    for(m in m_seq){
        if(group){
            NNlist <- group_obs(NNarray[,1:(m+1)])
            likfun <- function(logparms){
                likobj <- arma_vecchia_grad_hess_grouped(link(logparms),covfun_name,
                    yord,Xord,locsord,NNlist)
                likobj$loglik <- -likobj$loglik - pen(link(logparms))
                likobj$grad <- -c(likobj$grad)*dlink(logparms) -
                    dpen(link(logparms))*dlink(logparms)
                likobj$info <- likobj$info*outer(dlink(logparms),dlink(logparms)) -
                    ddpen(link(logparms))*outer(dlink(logparms),dlink(logparms))
                return(likobj)
            }
        } else {
            likfun <- function(logparms){
                likobj <- arma_vecchia_grad_hess(link(logparms),covfun_name,
                    yord,Xord,locsord,NNarray[,1:(m+1)])
                likobj$loglik <- -( likobj$loglik + pen(link(logparms)) )
                likobj$grad <- -( c(likobj$grad)*dlink(logparms) +
                    dpen(link(logparms))*dlink(logparms) )
                likobj$info <- likobj$info*outer(dlink(logparms),dlink(logparms)) -
                    ddpen(link(logparms))*outer(dlink(logparms),dlink(logparms))
                return(likobj)
            }
        }
        fit <- fisher_scoring(likfun,invlink(start_parms),link,silent=silent)
        fit$loglik <- -fit$loglik - pen(fit$covparms)
        start_parms <- fit$covparms
    }
    return(fit)
}





get_start_parms <- function(y,X,locs,covfun_name){

    fitlm <- stats::lm(y ~ X - 1 )
    start_var <- summary(fitlm)$sigma^2
    start_smooth <- 0.8
    start_nug <- 0.1
    n <- length(y)

    randinds <- sample(1:n, min(n,200))
    dmat <- fields::rdist(locs[randinds,])

    if(covfun_name == "arma_exponential_isotropic"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_nug)
    }
    if(covfun_name == "arma_matern_isotropic"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
    }
    if(covfun_name == "arma_matern_anisotropic2D"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, 1/start_range, 0, 1/start_range,
            start_smooth, start_nug)
    }
    if(covfun_name == "arma_exponential_anisotropic3D"){
        dmat <- fields::rdist(locs[randinds,1,drop=FALSE])
        start_range1 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,2,drop=FALSE])
        start_range2 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,3,drop=FALSE])
        start_range3 <- mean( dmat )/4
        start_parms <- c(start_var, 1/start_range1, 0, 1/start_range2,
            0, 0, 1/start_range3, start_nug )
    }
    if(covfun_name == "arma_matern_anisotropic3D"){
        dmat <- fields::rdist(locs[randinds,1,drop=FALSE])
        start_range1 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,2,drop=FALSE])
        start_range2 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,3,drop=FALSE])
        start_range3 <- mean( dmat )/4
        start_parms <- c(start_var, 1/start_range1, 0, 1/start_range2,
            0, 0, 1/start_range3, start_smooth, start_nug )
    }
    if(covfun_name == "arma_matern_sphere"){
        start_range <- 0.2
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
    }
    if(covfun_name == "arma_matern_sphere_warp"){
        # hard-code Lmax = 2
        start_range <- 0.2
        start_parms <- c(start_var, start_range, start_smooth, start_nug,
            rep(0,5))
    }
    if(covfun_name == "arma_matern_spheretime_warp"){
        # hard-code Lmax = 2
        start_range <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/20
        start_parms <- c(start_var, start_range, start_range2, start_smooth,
            start_nug, rep(0,5))
    }
    if(covfun_name == "arma_matern_spacetime"){
        d <- ncol(locs)-1
        dmat1 <- fields::rdist(locs[randinds,1:d])
        dmat2 <- fields::rdist(locs[randinds,d+1,drop=FALSE])
        start_range1 <- mean( dmat1 )/4
        start_range2 <- mean( dmat2 )/1
        start_parms <-
            c(start_var, start_range1, start_range2, start_smooth, start_nug)
    }
    if(covfun_name == "arma_matern_spheretime"){
        start_range1 <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/10
        start_parms <-
            c(start_var, start_range1, start_range2, start_smooth, start_nug)
    }
    if(covfun_name == "arma_matern_spheretime_nonstatvar"){
        start_range1 <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/10
        start_parms <-
            c(start_var, start_range1, start_range2, start_smooth, start_nug,
                rep(0,ncol(locs)-3))
    }
    if(covfun_name == "arma_matern_nonstat_var"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug,
            rep(0, ncol(locs)-2) )
    }
    return( list( start_parms = start_parms ) )
}


# starting values and covariance-specific settings
get_linkfun <- function(y,X,locs,covfun_name, silent = FALSE){

    link <- exp
    dlink <- exp
    invlink <- log
    lonlat <- FALSE
    space_time <- FALSE

    if(covfun_name == "arma_matern_anisotropic2D"){
        link <- function(x){    c( exp(x[1:2]), x[3], exp(x[4:6]) )  }
        dlink <- function(x){   c( exp(x[1:2]), 1.0,  exp(x[4:6]) ) }
        ddlink <- function(x){  c( exp(x[1:2]), 0.0,  exp(x[4:6]) ) }
        invlink <- function(x){ c( log(x[1:2]), x[3], log(x[4:6]) ) }
    }
    if(covfun_name == "arma_exponential_anisotropic3D"){
        link <- function(x)
        { c( exp(x[1:2]), x[3], exp(x[4]), x[5:6], exp(x[7:8]) ) }
        dlink <- function(x)
        { c( exp(x[1:2]), 1.0, exp(x[4]), 1.0,1.0, exp(x[7:8]) ) }
        invlink <- function(x)
        { c( log(x[1:2]), x[3], log(x[4]), x[5:6], log(x[7:8]) ) }
    }
    if(covfun_name == "arma_matern_anisotropic3D"){
        link <- function(x)
        { c( exp(x[1:2]), x[3], exp(x[4]), x[5:6], exp(x[7:9]) ) }
        dlink <- function(x)
        { c( exp(x[1:2]), 1.0, exp(x[4]), 1.0,1.0, exp(x[7:9]) ) }
        invlink <- function(x)
        { c( log(x[1:2]), x[3], log(x[4]), x[5:6], log(x[7:9]) ) }
    }
    if(covfun_name == "arma_matern_nonstat_var"){
        link <- function(x){    c( exp(x[1:3]), exp(x[4]), x[5:length(x)]     ) }
        dlink <- function(x){   c( exp(x[1:3]), exp(x[4]), rep(1,length(x)-4) ) }
        ddlink <- function(x){  c( exp(x[1:3]), exp(x[4]), rep(0,length(x)-4) ) }
        invlink <- function(x){ c( log(x[1:3]), log(x[4]), x[5:length(x)]     ) }
    }
    if(covfun_name == "arma_matern_spheretime_nonstatvar"){
        link <- function(x){    c( exp(x[1:5]), x[6:length(x)]     ) }
        dlink <- function(x){   c( exp(x[1:5]), rep(1,length(x)-5) ) }
        ddlink <- function(x){  c( exp(x[1:5]), rep(0,length(x)-5) ) }
        invlink <- function(x){ c( log(x[1:5]), x[6:length(x)]     ) }
    }

    if(covfun_name == "arma_matern_sphere"){ lonlat <- TRUE }
    if(covfun_name == "arma_matern_sphere_warp"){
        lonlat <- TRUE
        link <- function(x){ c(exp(x[1:4]), x[5:length(x)]) }
        dlink <- function(x){c(exp(x[1:4]), rep(1,length(x)-4))}
        invlink <- function(x){ c(log(x[1:4]),x[5:length(x)])}
    }
    if(covfun_name == "arma_matern_spheretime_warp"){
        lonlat <- TRUE
        link <- function(x){ c(exp(x[1:5]), x[6:length(x)]) }
        dlink <- function(x){c(exp(x[1:5]), rep(1,length(x)-5))}
        invlink <- function(x){ c(log(x[1:5]),x[6:length(x)])}
    }
    if(covfun_name == "arma_matern_spacetime"){ space_time <- FALSE }
    if(covfun_name == "arma_matern_spheretime"){
        lonlat <- TRUE
        space_time <- FALSE
    }

    return(list(
        link = link, dlink = dlink, invlink = invlink,
        lonlat = lonlat, space_time = space_time
    ))
}



# penalty functions
get_penalty <- function(y,X,locs,covfun_name, silent = FALSE){

    fitlm <- stats::lm(y ~ X - 1 )
    vv <- summary(fitlm)$sigma^2
    # by default, no penalty
    pen <- function(x) 0.0
    dpen <- function(x) rep(0,length(x))
    ddpen <- function(x) matrix(0,length(x),length(x))
    # nugget penalty
    pen_nug <- function(x,j){ pen_loglo(x[j],.1,log(0.001)) }
    dpen_nug <- function(x,j){
        dpen <- rep(0,length(x))
        dpen[j] <- dpen_loglo(x[j],.1,log(0.001))
        return(dpen)
    }
    ddpen_nug <- function(x,j){
        ddpen <- matrix(0,length(x),length(x))
        ddpen[j,j] <- ddpen_loglo(x[j],.1,log(0.001))
        return(ddpen)
    }
    # smoothness penalty
    pen_sm <- function(x,j){ pen_loglo(x[j],.1,log(0.2)) }
    dpen_sm <- function(x,j){
        dpen <- rep(0,length(x))
        dpen[j] <- dpen_loglo(x[j],.1,log(0.2))
        return(dpen)
    }
    ddpen_sm <- function(x,j){
        ddpen <- matrix(0,length(x),length(x))
        ddpen[j,j] <- ddpen_loglo(x[j],.1,log(0.2))
        return(ddpen)
    }
    # variance penalty
    # dangerous because vv could get redefined
    pen_var <- function(x,j){ pen_hi(x[j]/vv,1,6) }
    dpen_var <- function(x,j){
        dpen <- rep(0,length(x))
        dpen[j] <- 1/vv*dpen_hi(x[j]/vv,1,6)
        return(dpen)
    }
    ddpen_var <- function(x,j){
        ddpen <- matrix(0,length(x),length(x))
        ddpen[j,j] <- 1/vv^2*ddpen_hi(x[j]/vv,1,6)
        return(ddpen)
    }

    if(covfun_name == "arma_exponential_isotropic"){
          pen <- function(x){  pen_nug(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_isotropic"){
          pen <- function(x){  pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_anisotropic2D"){
          pen <- function(x){  pen_nug(x,6) +   pen_sm(x,5) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,6) +  dpen_sm(x,5) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,6) + ddpen_sm(x,5) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_exponential_anisotropic3D"){
          pen <- function(x){  pen_nug(x,8)  +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,8)  +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,8)  + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_anisotropic3D"){
          pen <- function(x){  pen_nug(x,9) +   pen_sm(x,8) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,9) +  dpen_sm(x,8) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,9) + ddpen_sm(x,8) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_sphere"){
          pen <- function(x){  pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_sphere_warp"){
          pen <- function(x){  pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_spheretime_warp"){
          pen <- function(x){  pen_nug(x,5) +   pen_sm(x,4) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,5) +  dpen_sm(x,4) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,5) + ddpen_sm(x,4) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_spacetime"){
          pen <- function(x){  pen_nug(x,5) +   pen_sm(x,4) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,5) +  dpen_sm(x,4) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,5) + ddpen_sm(x,4) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_spheretime"){
          pen <- function(x){  pen_nug(x,5) +   pen_sm(x,4) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,5) +  dpen_sm(x,4) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,5) + ddpen_sm(x,4) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_spheretime_nonstatvar"){
          pen <- function(x){  pen_nug(x,5) +   pen_sm(x,4) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,5) +  dpen_sm(x,4) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,5) + ddpen_sm(x,4) + ddpen_var(x,1) }
    }
    if(covfun_name == "arma_matern_nonstat_var"){
          pen <- function(x){  pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) }
    }
    return( list( pen = pen, dpen = dpen, ddpen = ddpen ) )
}


# starting values and covariance-specific settings
get_start_parms_linkfun <- function(y,X,locs,covfun_name, silent = FALSE){

    fitlm <- stats::lm(y ~ X - 1 )
    start_var <- summary(fitlm)$sigma^2
    start_smooth <- 0.8
    start_nug <- 0.1
    n <- length(y)




    randinds <- sample(1:n, min(n,200))
    if(covfun_name == "arma_exponential_isotropic"){
        lonlat <- FALSE
        space_time <- FALSE
        dmat <- fields::rdist(locs[randinds,])
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_nug)
        link <- exp
        dlink <- exp
        ddlink <- exp
        invlink <- log
    }

    if(covfun_name == "arma_matern_isotropic"){
        lonlat <- FALSE
        space_time <- FALSE
        dmat <- fields::rdist(locs[randinds,])
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
        link <- exp
        dlink <- exp
        ddlink <- exp
        invlink <- log
    }

    if(covfun_name == "arma_matern_anisotropic2D"){
        lonlat <- FALSE
        space_time <- FALSE
        dmat <- fields::rdist(locs[randinds,])
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, 1/start_range, 0, 1/start_range,
            start_smooth, start_nug)
        link <- function(x){    c( exp(x[1:2]), x[3], exp(x[4:6]) )  }
        dlink <- function(x){   c( exp(x[1:2]), 1.0,  exp(x[4:6]) ) }
        ddlink <- function(x){   c( exp(x[1:2]), 0.0,  exp(x[4:6]) ) }
        invlink <- function(x){ c( log(x[1:2]), x[3], log(x[4:6]) ) }
    }

    if(covfun_name == "arma_matern_anisotropic3D"){
        lonlat <- FALSE
        space_time <- TRUE
        dmat <- fields::rdist(locs[randinds,1,drop=FALSE])
        start_range1 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,2,drop=FALSE])
        start_range2 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,3,drop=FALSE])
        start_range3 <- mean( dmat )/4
        start_parms <- c(start_var, 1/start_range1, 0, 1/start_range2,
            0, 0, 1/start_range3, start_smooth, start_nug )
        link <- function(x)
        { c( exp(x[1:2]), x[3], exp(x[4]), x[5:6], exp(x[7:8]), x[9] ) }
        dlink <- function(x)
        { c( exp(x[1:2]), 1, exp(x[4]), 1, 1, exp(x[7:8]), 1 ) }
        invlink <- function(x)
        { c( log(x[1:2]), x[3], log(x[4]), x[5:6], log(x[7:8]), x[9] ) }
    }

    if(covfun_name == "arma_matern_sphere"){
        lonlat <- TRUE
        space_time <- FALSE
        #dmat <- fields::rdist(locs[randinds,])
        start_range <- 0.2
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
        link <- exp
        dlink <- exp
        invlink <- log
    }

    if(covfun_name == "arma_matern_sphere_warp"){
        # Lmax = 2 hard coded
        lonlat <- TRUE
        space_time <- FALSE
        start_range <- 0.2
        start_parms <- c(start_var, start_range, start_smooth, start_nug,
            rep(0,5))
        link <- function(x){ c(exp(x[1:4]), x[5:length(x)]) }
        dlink <- function(x){c(exp(x[1:4]), rep(1,length(x)-4))}
        invlink <- function(x){ c(log(x[1:4]),x[5:length(x)])}
    }
    if(covfun_name == "arma_matern_spheretime_warp"){
        # Lmax = 2 hard coded
        lonlat <- TRUE
        space_time <- TRUE
        start_range <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/20
        start_parms <- c(start_var, start_range, start_range2, start_smooth,
            start_nug, rep(0,5))
        link <- function(x){ c(exp(x[1:5]), x[6:length(x)]) }
        dlink <- function(x){c(exp(x[1:5]), rep(1,length(x)-5))}
        invlink <- function(x){ c(log(x[1:5]),x[6:length(x)])}
    }


    if(covfun_name == "arma_matern_spacetime"){
        lonlat <- FALSE
        space_time <- TRUE
        d <- ncol(locs)-1
        dmat1 <- fields::rdist(locs[randinds,1:d])
        dmat2 <- fields::rdist(locs[randinds,d+1,drop=FALSE])
        start_range1 <- mean( dmat1 )/4
        start_range2 <- mean( dmat2 )/1
        start_parms <-
            c(start_var, start_range1, start_range2, start_smooth, start_nug)
        link <- exp
        dlink <- exp
        ddlink <- exp
        invlink <- log
    }

    if(covfun_name == "arma_matern_spheretime"){
        lonlat <- TRUE
        space_time <- TRUE
        #dmat <- fields::rdist(locs[randinds,])
        start_range1 <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/10
        start_parms <-
            c(start_var, start_range1, start_range2, start_smooth, start_nug)
        link <- exp
        dlink <- exp
        invlink <- log
    }

    if(covfun_name == "arma_matern_spheretime_nonstatvar"){
        lonlat <- TRUE
        space_time <- FALSE
        start_range1 <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/10
        start_parms <-
            c(start_var, start_range1, start_range2, start_smooth, start_nug,
                rep(0,ncol(locs)-3))
        link <- function(x){    c( exp(x[1:5]), x[6:length(x)]     ) }
        dlink <- function(x){   c( exp(x[1:5]), rep(1,length(x)-5) ) }
        ddlink <- function(x){  c( exp(x[1:5]), rep(0,length(x)-5) ) }
        invlink <- function(x){ c( log(x[1:5]), x[6:length(x)]     ) }
    }

    if(covfun_name == "arma_matern_nonstat_var"){
        lonlat <- FALSE
        space_time <- FALSE
        dmat <- fields::rdist(locs[randinds,])
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug,
            rep(0, ncol(locs)-2) )
        link <- function(x){    c( exp(x[1:3]), exp(x[4]), x[5:length(x)]     ) }
        dlink <- function(x){   c( exp(x[1:3]), exp(x[4]), rep(1,length(x)-4) ) }
        ddlink <- function(x){  c( exp(x[1:3]), exp(x[4]), rep(0,length(x)-4) ) }
        invlink <- function(x){ c( log(x[1:3]), log(x[4]), x[5:length(x)]     ) }
    }

    if(covfun_name == "matern_isotropic"){
        lonlat <- FALSE
        space_time <- FALSE
        dmat <- fields::rdist(locs[randinds,])
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
        link <- exp
        dlink <- exp
        invlink <- log
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
        link <- exp
        dlink <- exp
        invlink <- log
    }

    if( covfun_name == "matern_sphere" ){
        lonlat <- TRUE
        space_time <- FALSE
        dmat <- fields::rdist.earth(locs[randinds,], R = 1)
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
        link <- exp
        dlink <- exp
        invlink <- log
        if(!silent){
            cat("Using 'matern_sphere_time'.\n")
            cat("Assuming that argument 'locs' is (longitude,latitude,time)\n")
        }
    }

    if( covfun_name == "matern_sphere_time" ){
        lonlat <- TRUE
        space_time <- TRUE
        dmat <- fields::rdist.earth(locs[randinds,1:2], R = 1)
        start_range1 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,3])
        start_range2 <- mean( dmat )/4
        start_parms <- c(start_var, start_range1, start_range2, start_smooth, start_nug)
        link <- exp
        dlink <- exp
        invlink <- log
        if(!silent){
            cat("Using 'matern_sphere_time'.\n")
            cat("Assuming that argument 'locs' is (longitude,latitude,time)\n")
        }
    }
    return(list(
        start_parms = start_parms,
        link = link, dlink = dlink, invlink = invlink,
        lonlat = lonlat,space_time = space_time
    ))
}
