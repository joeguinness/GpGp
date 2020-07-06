

#' Estimate mean and covariance parameters
#'
#' Given a response, set of locations, (optionally) a design matrix,
#' and a specified covariance function, return the maximum
#' Vecchia likelihood estimates, obtained with a Fisher scoring algorithm.
#'
#' @param y response vector
#' @param locs matrix of locations. Each row is a single 
#' spatial or spatial-temporal
#' location. If using one of the covariance functions for data on a sphere,
#' the first column should be longitudes (-180,180) and the second column
#' should be latitudes (-90,90). 
#' If using a spatial-temporal covariance function,
#' the last column should contain the times.
#' @param X design matrix. Each row contains covariates for the corresponding
#' observation in \code{y}. If not specified, the function sets \code{X} to be a
#' matrix with a single column of ones, that is, a constant mean function.
#' @param covfun_name string name of a covariance function. 
#' See \code{\link{GpGp}} for information about supported covariance funtions.
#' @param NNarray Optionally specified array of nearest neighbor indices, 
#' usually from the output of \code{\link{find_ordered_nn}}. If \code{NULL}, 
#' fit_model will compute the nearest neighbors. We recommend that the user
#' not specify this unless there is a good reason to (e.g. if doing a comparison
#' study where one wants to control 
#' \code{NNarray} across different approximations).
#' @param start_parms Optionally specified starting values for parameters. 
#' If \code{NULL},
#' fit_model will select default starting values.
#' @param max_iter maximum number of Fisher scoring iterations
#' @param silent TRUE/FALSE for whether to print some information during fitting.
#' @param group TRUE/FALSE for whether to use the grouped version of
#' the approximation (Guinness, 2018) or not.  The grouped version
#' is used by default and is always recommended.
#' @param reorder TRUE/FALSE indicating whether maxmin ordering should be used
#' (TRUE) or whether no reordering should be done before fitting (FALSE). 
#' If you want
#' to use a customized reordering, then manually reorder \code{y}, \code{locs}, 
#' and \code{X},
#' and then set \code{reorder} to \code{FALSE}. A random reordering is used
#' when \code{nrow(locs) > 1e5}.
#' @param m_seq Sequence of values for number of neighbors. By default, 
#' a 10-neighbor
#' approximation is maximized, then a 30-neighbor approximation is 
#' maximized using the
#' 10 neighbor estimates as starting values. 
#' However, one can specify any sequence
#' of numbers of neighbors, e.g. \code{m_seq = c(10,30,60,90)}.
#' @param fixed_parms Indices of covariance parameters you would like to fix
#' at specific values. If you decide to fix any parameters, you must specify
#' their values in \code{start_parms}, along with the starting values for
#' all other parameters. For example, to fix the nugget at zero in 
#' \code{exponential_isotropic}, set \code{fixed_parms} to \code{c(3)}, and set
#' \code{start_parms} to \code{c(4.7,3.1,0)}. The
#' last element of \code{start_parms} (the nugget parameter) is set to zero,
#' while the starting values for the other two parameters are 4.7 and 3.1.
#' @param st_scale Scaling for spatial and temporal ranges. Only applicable for
#' spatial-temporal models, where it is used in distance
#' calculations when selecting neighbors. \code{st_scale} must be specified
#' when \code{covfun_name} is a spatial-temporal covariance. 
#' See Argo vignette for an example.
#' @param convtol Tolerance for exiting the optimization. 
#' Fisher scoring is stopped
#' when the dot product between the step and the gradient 
#' is less than \code{convtol}.
#' @return An object of class \code{GpGp_fit}, which is a list containing
#' covariance parameter estimates, regression coefficients,
#' covariance matrix for mean parameter estimates, as well as some other
#' information relevant to the model fit.
#' @details \code{fit_model} is a user-friendly model fitting function
#' that automatically performs many of the auxiliary tasks needed for
#' using Vecchia's approximation, including reordering, computing
#' nearest neighbors, grouping, and optimization. The likelihoods use a small
#' penalty on small nuggets, large spatial variances, 
#' and small smoothness parameter.
#'
#' The Jason-3 windspeed vignette and the Argo temperature 
#' vignette are useful sources for a
#' use-cases of the \code{fit_model} function for data on sphere. 
#' The example below shows a very small example with a simulated dataset in 2d.
#'
#' @examples
#' n1 <- 20
#' n2 <- 20
#' n <- n1*n2
#' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
#' covparms <- c(2,0.1,1/2,0)
#' y <- 7 + fast_Gp_sim(covparms, "matern_isotropic", locs)
#' X <- as.matrix( rep(1,n) )
#' ## not run
#' # fit <- fit_model(y, locs, X, "matern_isotropic")
#' # fit
#'
#'
#' @export
fit_model <- function(y, locs, X = NULL, covfun_name = "matern_isotropic",
    NNarray = NULL, start_parms = NULL, reorder = TRUE, group = TRUE,
    m_seq = c(10,30), max_iter = 40, fixed_parms = NULL,
    silent = FALSE, st_scale = NULL, convtol = 1e-4){

    n <- length(y)

    # check that length of observation vector same as
    # number of locations
    if( nrow(locs) != n ){
        stop("length of observation vector y not equal
              to the number of locations (rows in locs)")
    }
    locs <- as.matrix(locs)

    # check if design matrix is specified
    if( is.null(X) ){
        if(!silent) cat("Design matrix not specified, using constant mean \n")
        X <- rep(1,n)
    }
    X <- as.matrix(X)

    # check if one of the allowed covariance functions is chosen
    if( ! covfun_name %in%
            c("exponential_isotropic",
              "matern_isotropic",
              "matern15_isotropic",
              "matern25_isotropic",
              "matern35_isotropic",
              "matern45_isotropic",
              "matern_anisotropic2D",
              "exponential_anisotropic2D",
              "exponential_anisotropic3D",
              "matern_anisotropic3D",
              "matern_nonstat_var",
              "exponential_nonstat_var",
              "matern_sphere",
              "exponential_sphere",
              "matern_sphere_warp",
              "exponential_sphere_warp",
              "matern_spheretime_warp",
              "exponential_spheretime_warp",
              "matern_spheretime",
              "exponential_spheretime",
              "matern_spacetime",
              "exponential_spacetime",
              "matern_scaledim",
              "matern15_scaledim",
              "matern25_scaledim",
              "matern35_scaledim",
              "matern45_scaledim",
              "exponential_scaledim" ) )
    {
        stop("unrecognized covariance function name `covfun_name'.")
    }

    # detect and remove missing values
    not_missing <- apply( cbind(y,locs,X), 1,
        function(x){
            if( sum(is.na(x) | is.infinite(x)) > 0 ){
                return(FALSE)
            } else { return(TRUE) }
        }
    )
    if( sum(not_missing) < n ){
        y <- y[not_missing]
        locs <- locs[not_missing,,drop=FALSE]
        X <- X[not_missing,,drop=FALSE]
        cat(paste0( n - sum(not_missing),
            " observations removed due to missingness or Inf\n"))
    }

    # redefine n
    n <- length(y)
    
    # check that start_parms is specified when fixed_parms is
    if( is.null(fixed_parms) ){
        if( is.null(start_parms) ){
            start <- get_start_parms(y,X,locs,covfun_name)
            start_parms <- start$start_parms
        } else {
            # check if start_parms has the right length
            start <- get_start_parms(y,X,locs,covfun_name)
            if(length(start_parms) != length(start$start_parms) ){
                stop(paste0("start_parms not correct length for ",covfun_name))
            }
        }
        # define the parameters we are not fixing
        active <- rep(TRUE, length(start_parms) )
    } else {
        if( is.null(start_parms) ){
            stop("start_parms must be specified whenever fixed_parms is")
        }
        # check if start_parms has the right length
        start <- get_start_parms(y,X,locs,covfun_name)
        if(length(start_parms) != length(start$start_parms) ){
            stop(paste0("start parms not correct length for ",covfun_name))
        }
        # check whether fixed_parms has appropriate values
        if( max( fixed_parms - floor(fixed_parms) ) > 0 ){
            stop("fixed_parms must contain indices of parms you want to fix")
        }
        if( min( fixed_parms < 1 ) || max(fixed_parms) > length(start_parms) ){
            stop("fixed_parms must be between 1 and number of parameters")
        }
        # define the parameters we are not fixing
        active <- rep(TRUE, length(start_parms) )
        active[fixed_parms] <- FALSE
    } 

    # get link functions
    linkfuns <- get_linkfun(covfun_name)
    link <- linkfuns$link
    dlink <- linkfuns$dlink
    invlink <- linkfuns$invlink
    invlink_startparms <- invlink(start_parms)
    lonlat <- linkfuns$lonlat
    if(lonlat){
    cat("Assuming columns 1 and 2 of locs are (longitude,latidue) in degrees\n")
    }
    space_time <- linkfuns$space_time

    penalty <- get_penalty(y,X,locs,covfun_name) 
    pen <- penalty$pen
    dpen <- penalty$dpen
    ddpen <- penalty$ddpen

    # get an ordering and reorder everything
    if(reorder){
        if(!silent) cat("Reordering...")
        if( n < 1e5 ){  # maxmin ordering if n < 100000
            ord <- order_maxmin(locs, lonlat = lonlat, space_time = space_time)
        } else {        # otherwise random order
            ord <- sample(n)
        }
        if(!silent) cat("Done \n")
    } else {
        ord <- 1:n
    }
    yord <- y[ord]
    locsord <- locs[ord,,drop=FALSE]
    Xord <- as.matrix( X[ord,,drop=FALSE] )

    # get neighbor array if not provided
    if( is.null(NNarray) ){
        if(!silent) cat("Finding nearest neighbors...")
        NNarray <- find_ordered_nn(locsord, m=max(m_seq), lonlat = lonlat,
            st_scale = st_scale)
        if(!silent) cat("Done \n")
    }

    # refine the estimates using m in m_seq
    for(i in 1:length(m_seq)){
        m <- m_seq[i]
        if(group){

            NNlist <- group_obs(NNarray[,1:(m+1)])
            likfun <- function(logparms){
                
                lp <- rep(NA,length(start_parms))
                lp[active] <- logparms
                lp[!active] <- invlink_startparms[!active]
                
                likobj <- vecchia_grouped_profbeta_loglik_grad_info(
                    link(lp),covfun_name,yord,Xord,locsord,NNlist)
                likobj$loglik <- -likobj$loglik - pen(link(lp))
                likobj$grad <- -c(likobj$grad)*dlink(lp) -
                    dpen(link(lp))*dlink(lp)
                likobj$info <- likobj$info*outer(dlink(lp),dlink(lp)) -
                    ddpen(link(lp))*outer(dlink(lp),dlink(lp))
                likobj$grad <- likobj$grad[active]
                likobj$info <- likobj$info[active,active]
                return(likobj)
                
            }

        } else {

            likfun <- function(logparms){

                lp <- rep(NA,length(start_parms))
                lp[active] <- logparms
                lp[!active] <- invlink_startparms[!active]
                
                likobj <- vecchia_profbeta_loglik_grad_info(
                    link(lp),covfun_name,yord,Xord,locsord,NNarray[,1:(m+1)])
                likobj$loglik <- -likobj$loglik - pen(link(lp))
                likobj$grad <- -c(likobj$grad)*dlink(lp) -
                    dpen(link(lp))*dlink(lp)
                likobj$info <- likobj$info*outer(dlink(lp),dlink(lp)) -
                    ddpen(link(lp))*outer(dlink(lp),dlink(lp))
                likobj$grad <- likobj$grad[active]
                likobj$info <- likobj$info[active,active]
                return(likobj)
            }
        }
        fit <- fisher_scoring( likfun,invlink(start_parms)[active],
            link,silent=silent, convtol = convtol, max_iter = max_iter )
        invlink_startparms[active] <- fit$logparms
        #start_parms[active] <- fit$covparms
        start_parms <- link(invlink_startparms)
        fit$loglik <- -fit$loglik - pen(start_parms)
        invlink_startparms <- invlink(start_parms)
    }

    # return fit and information used for predictions
    fit$covfun_name <- covfun_name
    #fit$covparms <- start_parms
    lp <- rep(NA,length(start_parms))
    lp[active] <- fit$logparms
    lp[!active] <- invlink_startparms[!active]
    fit$covparms <- link(lp)
    fit$y <- y
    fit$locs <- locs
    fit$X <- X
    class(fit) <- "GpGp_fit"
    return(fit)
}



#' Print summary of GpGp fit
#'
#' @param object Object of class "GpGp_fit", usually the return value from
#' \code{\link{fit_model}}
#' @param ... additional arguments, for compatability with S3 generic 'summary'
#' @export
summary.GpGp_fit <- function(object, ...){
    cat(paste0("Covariance Function: ",object$covfun_name,"\n\n"))
    cat(paste0("Covariance Parameters: \n"))
    cat(paste0(round(object$covparms,4)),"\n\n")
    cat(paste0("Loglikelihood: ",round(object$loglik,4),"\n\n"))
    X <- as.data.frame(object$X)
    df <- data.frame(
        variable = colnames(X),
        estimate = round(object$betahat,4),
        std_error = round(object$sebeta,4),
        t_stat = round(object$tbeta,4)
    )
    rownames(df) <- c()
    cat("Linear Mean Parameters: \n")
    print(df)
    cat("\n")

}

#' get default starting values of covariance parameters
#'
#' @param y response
#' @param X design matrix
#' @param locs locations
#' @param covfun_name string name of covariance function
get_start_parms <- function(y,X,locs,covfun_name){

    fitlm <- stats::lm(y ~ X - 1 )
    start_var <- summary(fitlm)$sigma^2
    start_smooth <- 0.8
    start_nug <- 0.1
    n <- length(y)

    randinds <- sample(1:n, min(n,200))
    dmat <- fields::rdist(locs[randinds,])

    if(covfun_name == "exponential_isotropic"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_nug)
    }
    if(covfun_name == "matern_isotropic"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
    }
    if(covfun_name == "matern15_isotropic"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_nug)
    }
    if(covfun_name == "matern25_isotropic"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_nug)
    }
    if(covfun_name == "matern35_isotropic"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_nug)
    }
    if(covfun_name == "matern45_isotropic"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_nug)
    }
    if(covfun_name == "matern_anisotropic2D"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, 1/start_range, 0, 1/start_range,
            start_smooth, start_nug)
    }
    if(covfun_name == "exponential_anisotropic2D"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, 1/start_range, 0, 1/start_range, start_nug)
    }
    if(covfun_name == "exponential_anisotropic3D"){
        dmat <- fields::rdist(locs[randinds,1,drop=FALSE])
        start_range1 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,2,drop=FALSE])
        start_range2 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,3,drop=FALSE])
        start_range3 <- mean( dmat )/4
        start_parms <- c(start_var, 1/start_range1, 0, 1/start_range2,
            0, 0, 1/start_range3, start_nug )
    }
    if(covfun_name == "matern_anisotropic3D"){
        dmat <- fields::rdist(locs[randinds,1,drop=FALSE])
        start_range1 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,2,drop=FALSE])
        start_range2 <- mean( dmat )/4
        dmat <- fields::rdist(locs[randinds,3,drop=FALSE])
        start_range3 <- mean( dmat )/4
        start_parms <- c(start_var, 1/start_range1, 0, 1/start_range2,
            0, 0, 1/start_range3, start_smooth, start_nug )
    }
    if(covfun_name == "matern_sphere"){
        start_range <- 0.2
        start_parms <- c(start_var, start_range, start_smooth, start_nug)
    }
    if(covfun_name == "exponential_sphere"){
        start_range <- 0.2
        start_parms <- c(start_var, start_range, start_nug)
    }
    if(covfun_name == "matern_sphere_warp"){
        # hard-code Lmax = 2
        start_range <- 0.2
        start_parms <- c(start_var, start_range, start_smooth, start_nug,
            rep(0,5))
    }
    if(covfun_name == "exponential_sphere_warp"){
        # hard-code Lmax = 2
        start_range <- 0.2
        start_parms <- c(start_var, start_range, start_nug, rep(0,5))
    }
    if(covfun_name == "matern_spheretime_warp"){
        # hard-code Lmax = 2
        start_range <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/20
        start_parms <- c(start_var, start_range, start_range2, start_smooth,
            start_nug, rep(0,5))
    }
    if(covfun_name == "exponential_spheretime_warp"){
        # hard-code Lmax = 2
        start_range <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/20
        start_parms <- c(start_var, start_range, start_range2,
            start_nug, rep(0,5))
    }
    if(covfun_name == "matern_scaledim"){
        d <- ncol(locs)
        start_parms <- c(start_var)
        for(j in 1:d){
            dmat <- fields::rdist(locs[randinds,j])
            start_parms <- c(start_parms, stats::median(dmat)/4 )
        }
        start_parms <- c(start_parms, start_smooth, start_nug)
    }
    if(covfun_name == "exponential_scaledim"){
        d <- ncol(locs)
        start_parms <- c(start_var)
        for(j in 1:d){
            dmat <- fields::rdist(locs[randinds,j])
            start_parms <- c(start_parms, stats::median(dmat)/4 )
        }
        start_parms <- c(start_parms, start_nug)
    }
    if(covfun_name %in% 
            c("matern15_scaledim","matern25_scaledim",
                "matern35_scaledim","matern45_scaledim")
    ){
        d <- ncol(locs)
        start_parms <- c(start_var)
        for(j in 1:d){
            dmat <- fields::rdist(locs[randinds,j])
            start_parms <- c(start_parms, stats::median(dmat)/4 )
        }
        start_parms <- c(start_parms, start_nug)
    }
    if(covfun_name == "matern_spacetime"){
        d <- ncol(locs)-1
        dmat1 <- fields::rdist(locs[randinds,1:d])
        dmat2 <- fields::rdist(locs[randinds,d+1,drop=FALSE])
        start_range1 <- mean( dmat1 )/4
        start_range2 <- mean( dmat2 )/1
        start_parms <-
            c(start_var, start_range1, start_range2, start_smooth, start_nug)
    }
    if(covfun_name == "exponential_spacetime"){
        d <- ncol(locs)-1
        dmat1 <- fields::rdist(locs[randinds,1:d])
        dmat2 <- fields::rdist(locs[randinds,d+1,drop=FALSE])
        start_range1 <- mean( dmat1 )/4
        start_range2 <- mean( dmat2 )/1
        start_parms <- c(start_var, start_range1, start_range2, start_nug)
    }
    if(covfun_name == "matern_spheretime"){
        start_range1 <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/10
        start_parms <-
            c(start_var, start_range1, start_range2, start_smooth, start_nug)
    }
    if(covfun_name == "exponential_spheretime"){
        start_range1 <- 0.2
        start_range2 <- abs(diff(range(locs[,3])))/10
        start_parms <-
            c(start_var, start_range1, start_range2, start_nug)
    }
    if(covfun_name == "matern_nonstat_var"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_smooth, start_nug,
            rep(0, ncol(locs)-2) )
    }
    if(covfun_name == "exponential_nonstat_var"){
        start_range <- mean( dmat )/4
        start_parms <- c(start_var, start_range, start_nug, rep(0,ncol(locs)-2))
    }
    return( list( start_parms = start_parms ) )
}


#' get link function, whether locations are lonlat and space time
#'
#' @param covfun_name string name of covariance function
get_linkfun <- function(covfun_name){

    link <- exp
    dlink <- exp
    invlink <- log
    lonlat <- FALSE
    space_time <- FALSE

    if(covfun_name == "matern_anisotropic2D"){
        link <- function(x){    c( exp(x[1:2]), x[3], exp(x[4:6]) ) }
        dlink <- function(x){   c( exp(x[1:2]), 1.0,  exp(x[4:6]) ) }
        ddlink <- function(x){  c( exp(x[1:2]), 0.0,  exp(x[4:6]) ) }
        invlink <- function(x){ c( log(x[1:2]), x[3], log(x[4:6]) ) }
    }
    if(covfun_name == "exponential_anisotropic2D"){
        link <- function(x){    c( exp(x[1:2]), x[3], exp(x[4:5]) )  }
        dlink <- function(x){   c( exp(x[1:2]), 1.0,  exp(x[4:5]) ) }
        ddlink <- function(x){  c( exp(x[1:2]), 0.0,  exp(x[4:5]) ) }
        invlink <- function(x){ c( log(x[1:2]), x[3], log(x[4:5]) ) }
    }
    if(covfun_name == "exponential_anisotropic3D"){
        link <- function(x)
        { c( exp(x[1:2]), x[3], exp(x[4]), x[5:6], exp(x[7:8]) ) }
        dlink <- function(x)
        { c( exp(x[1:2]), 1.0, exp(x[4]), 1.0,1.0, exp(x[7:8]) ) }
        invlink <- function(x)
        { c( log(x[1:2]), x[3], log(x[4]), x[5:6], log(x[7:8]) ) }
    }
    if(covfun_name == "matern_anisotropic3D"){
        link <- function(x)
        { c( exp(x[1:2]), x[3], exp(x[4]), x[5:6], exp(x[7:9]) ) }
        dlink <- function(x)
        { c( exp(x[1:2]), 1.0, exp(x[4]), 1.0,1.0, exp(x[7:9]) ) }
        invlink <- function(x)
        { c( log(x[1:2]), x[3], log(x[4]), x[5:6], log(x[7:9]) ) }
    }
    if(covfun_name == "matern_nonstat_var"){
        link <- function(x){    c( exp(x[1:3]), exp(x[4]), x[5:length(x)]     )}
        dlink <- function(x){   c( exp(x[1:3]), exp(x[4]), rep(1,length(x)-4) )}
        ddlink <- function(x){  c( exp(x[1:3]), exp(x[4]), rep(0,length(x)-4) )}
        invlink <- function(x){ c( log(x[1:3]), log(x[4]), x[5:length(x)]     )}
    }
    if(covfun_name == "exponential_nonstat_var"){
        link <- function(x){    c( exp(x[1:3]), x[4:length(x)]     ) }
        dlink <- function(x){   c( exp(x[1:3]), rep(1,length(x)-3) ) }
        ddlink <- function(x){  c( exp(x[1:3]), rep(0,length(x)-3) ) }
        invlink <- function(x){ c( log(x[1:3]), x[4:length(x)]     ) }
    }
    if(covfun_name == "matern_sphere"){ lonlat <- TRUE }
    if(covfun_name == "exponential_sphere"){ lonlat <- TRUE }
    if(covfun_name == "matern_sphere_warp"){
        lonlat <- TRUE
        link <- function(x){ c(exp(x[1:4]), x[5:length(x)]) }
        dlink <- function(x){c(exp(x[1:4]), rep(1,length(x)-4))}
        invlink <- function(x){ c(log(x[1:4]),x[5:length(x)])}
    }
    if(covfun_name == "exponential_sphere_warp"){
        lonlat <- TRUE
        link <- function(x){ c(exp(x[1:3]), x[4:length(x)]) }
        dlink <- function(x){c(exp(x[1:3]), rep(1,length(x)-3))}
        invlink <- function(x){ c(log(x[1:3]),x[4:length(x)])}
    }
    if(covfun_name == "matern_spheretime_warp"){
        lonlat <- TRUE
        space_time <- TRUE
        link <- function(x){ c(exp(x[1:5]), x[6:length(x)]) }
        dlink <- function(x){c(exp(x[1:5]), rep(1,length(x)-5))}
        invlink <- function(x){ c(log(x[1:5]),x[6:length(x)])}
    }
    if(covfun_name == "exponential_spheretime_warp"){
        lonlat <- TRUE
        space_time <- TRUE
        link <- function(x){ c(exp(x[1:4]), x[5:length(x)]) }
        dlink <- function(x){c(exp(x[1:4]), rep(1,length(x)-4))}
        invlink <- function(x){ c(log(x[1:4]),x[5:length(x)])}
    }
    if(covfun_name == "matern_spacetime"){ space_time <- TRUE }
    if(covfun_name == "exponential_spacetime"){ space_time <- TRUE }
    if(covfun_name == "matern_spheretime"){
        lonlat <- TRUE
        space_time <- TRUE
    }
    if(covfun_name == "exponential_spheretime"){
        lonlat <- TRUE
        space_time <- TRUE
    }

    return(list(
        link = link, dlink = dlink, invlink = invlink,
        lonlat = lonlat, space_time = space_time
    ))
}



#' get penalty function
#'
#' @inheritParams get_start_parms
get_penalty <- function(y,X,locs,covfun_name){

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

    if(covfun_name == "exponential_isotropic"){
          pen <- function(x){  pen_nug(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name %in% c("matern15_isotropic","matern25_isotropic",
        "matern35_isotropic","matern45_isotropic")){
          pen <- function(x){  pen_nug(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_isotropic"){
          pen <- function(x){  pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_anisotropic2D"){
          pen <- function(x){  pen_nug(x,6) +   pen_sm(x,5) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,6) +  dpen_sm(x,5) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,6) + ddpen_sm(x,5) + ddpen_var(x,1) }
    }
    if(covfun_name == "exponential_anisotropic2D"){
          pen <- function(x){  pen_nug(x,5)  +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,5)  +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,5)  + ddpen_var(x,1) }
    }
    if(covfun_name == "exponential_anisotropic3D"){
          pen <- function(x){  pen_nug(x,8)  +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,8)  +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,8)  + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_anisotropic3D"){
          pen <- function(x){  pen_nug(x,9) +   pen_sm(x,8) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,9) +  dpen_sm(x,8) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,9) + ddpen_sm(x,8) + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_sphere"){
          pen <- function(x){  pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name == "exponential_sphere"){
          pen <- function(x){  pen_nug(x,3)  +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,3)  +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,3)  + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_sphere_warp"){
          pen <- function(x){  pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name == "exponential_sphere_warp"){
          pen <- function(x){  pen_nug(x,3)  +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,3)  +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,3)  + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_spheretime_warp"){
          pen <- function(x){  pen_nug(x,5) +   pen_sm(x,4) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,5) +  dpen_sm(x,4) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,5) + ddpen_sm(x,4) + ddpen_var(x,1) }
    }
    if(covfun_name == "exponential_spheretime_warp"){
          pen <- function(x){  pen_nug(x,4)  +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4)  +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4)  + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_scaledim"){
        d <- ncol(locs)
          pen <- function(x){  pen_nug(x,d+3) +   pen_sm(x,d+2) +   pen_var(x,1)}
         dpen <- function(x){ dpen_nug(x,d+3) +  dpen_sm(x,d+2) +  dpen_var(x,1)}
        ddpen <- function(x){ ddpen_nug(x,d+3) + ddpen_sm(x,d+2) + ddpen_var(x,1)}
    }
    if(covfun_name == "exponential_scaledim"){
        d <- ncol(locs)
          pen <- function(x){  pen_nug(x,d+2)  +   pen_var(x,1)   }
         dpen <- function(x){ dpen_nug(x,d+2)  +  dpen_var(x,1)  }
        ddpen <- function(x){ ddpen_nug(x,d+2) + ddpen_var(x,1) }
    }
    if(covfun_name %in% c("matern15_scaledim","matern25_scaledim",
        "matern35_scaledim","matern45_scaledim")
    ){
        d <- ncol(locs)
          pen <- function(x){  pen_nug(x,d+2)  +   pen_var(x,1)   }
         dpen <- function(x){ dpen_nug(x,d+2)  +  dpen_var(x,1)  }
        ddpen <- function(x){ ddpen_nug(x,d+2) + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_spacetime"){
          pen <- function(x){  pen_nug(x,5) +   pen_sm(x,4) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,5) +  dpen_sm(x,4) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,5) + ddpen_sm(x,4) + ddpen_var(x,1) }
    }
    if(covfun_name == "exponential_spacetime"){
          pen <- function(x){  pen_nug(x,4)  +   pen_var(x,1)   }
         dpen <- function(x){ dpen_nug(x,4)  +  dpen_var(x,1)  }
        ddpen <- function(x){ ddpen_nug(x,4) + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_spheretime"){
          pen <- function(x){  pen_nug(x,5) +   pen_sm(x,4) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,5) +  dpen_sm(x,4) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,5) + ddpen_sm(x,4) + ddpen_var(x,1) }
    }
    if(covfun_name == "exponential_spheretime"){
          pen <- function(x){  pen_nug(x,4)  +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4)  +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4)  + ddpen_var(x,1) }
    }
    if(covfun_name == "matern_nonstat_var"){
          pen <- function(x){  pen_nug(x,4) +   pen_sm(x,3) +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,4) +  dpen_sm(x,3) +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,4) + ddpen_sm(x,3) + ddpen_var(x,1) }
    }
    if(covfun_name == "exponential_nonstat_var"){
          pen <- function(x){  pen_nug(x,3)  +   pen_var(x,1)   }
         dpen <- function(x){  dpen_nug(x,3)  +  dpen_var(x,1)  }
        ddpen <- function(x){  ddpen_nug(x,3)  + ddpen_var(x,1) }
    }
    return( list( pen = pen, dpen = dpen, ddpen = ddpen ) )
}
