
#' @export
fit_model <- function(y, locs, X = NULL, covfun_name = "matern_isotropic",
    silent = FALSE, group = TRUE){
    
    n <- length(y)
    
    
    # check that length of observation vector same as
    # number of locations
    if( nrow(locs) != n ){
        stop("length of observation vector y not equal
              to the number of locations (rows in locs)")
    }
    
    # check if design matrix is specified
    if( is.null(X) ){
        if( !silent ){
            if(!silent) cat("Design matrix not specified, using constant mean \n")
        }
        X <- rep(1,n) 
    }
    X <- as.matrix(X)
    
    # check if one of the allowed covariance functions is chosen
    if( ! covfun_name %in% c("matern_isotropic","matern_sphere","matern_sphere_time") ){
        stop("unrecognized covariance function name `covfun_name'. Choose from
             'matern_isotropic', 'matern_sphere', or 'matern_sphere_time' " )
    }
    
    # missing values???
    
    # starting values and covariance-specific settings
    start_var <- var(y)
    start_smooth <- 1
    start_nug <- 0.1
    
    randinds <- sample(1:n,100)    
    if(covfun_name == "matern_isotropic"){
        lonlat <- FALSE
        dmat <- fields::rdist(locs[randinds,])
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
    }

    if( covfun_name == "matern_sphere" ){
        lonlat <- TRUE
        dmat <- fields::rdist.earth(locs[randinds,], R = 1)
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
        if(!silent) cat("Using 'matern_sphere'. Assuming that argument 'locs' is (longitude,latitude)\n")    
    }
    
    if( covfun_name == "matern_sphere_time" ){
        lonlat <- TRUE
        lonlat <- TRUE
        dmat <- fields::rdist.earth(locs[randinds,1:2], R = 1)
        start_range1 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,3])
        start_range2 <- mean( dmat )/4
        start_parms <- c(start_var, start_range1, start_range2, start_smooth, start_nug)
        if(!silent) cat("Using 'matern_sphere_time'. Assuming that argument 'locs' is (longitude,latitude,time)\n")    
    }
    
    nparms <- length(start_parms)
    
    # get an ordering and reorder everything
    if(!silent) cat("Reordering...")
    ord <- order_maxmin(locs, lonlat = lonlat)
    yord <- y[ord]
    locsord <- locs[ord,]
    Xord <- as.matrix( X[ord,] )
    if(!silent) cat("Done \n")
    
    # get nearest neighbors    
    if(!silent) cat("Finding nearest neighbors...")
    NNarray <- find_ordered_nn(locsord, m=45, lonlat = lonlat)
    if(!silent) cat("Done \n")
    
    fit <- list(par=log(start_parms[2:nparms]))
    
    # refine the estimates for m = c(15,30,45)
    for(m in c(5,15,30)){
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
        fit <- optim(fit$par,funtomax,
            control=list(trace=0,maxit=25), method = "BFGS") 
        if(!silent) cat("Done          ")
        if(!silent) cat(paste(  paste( round(exp(fit$par),3), collapse = " " ), "\n" ) )
    }
    fitted_model <- proflik_mean_variance_grouped(
            exp(fit$par), covfun_name, yord, Xord, locsord, NNlist, return_parms = TRUE )
        
    return(fitted_model)
}