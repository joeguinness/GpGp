
context("Covariance Functions")

get_test_locs <- function(covfun_name,n){
    
    nside <- round(sqrt(n))
    longrid <- seq(1,360,length.out=nside)
    latgrid <- seq(-80,80,length.out=nside)
    lonlat <- as.matrix(expand.grid(longrid,latgrid))
    if(covfun_name=="arma_exponential_isotropic"){
        locs <- matrix(runif(3*n),n,3)           
    } else if(covfun_name=="arma_matern_isotropic"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="arma_matern_spacetime"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="arma_matern_anisotropic2D"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="arma_exponential_anisotropic3D"){
        locs <- matrix(runif(3*n),n,3)           
    } else if(covfun_name=="arma_matern_anisotropic3D"){
        locs <- matrix(runif(3*n),n,3)           
    } else if(covfun_name=="arma_matern_nonstat_var"){
        locs <- matrix(runif(6*n),n,6)           
    } else if(covfun_name=="arma_matern_sphere"){
        locs <- lonlat
    } else if(covfun_name=="arma_matern_sphere_warp"){
        locs <- lonlat
    } else if(covfun_name=="arma_matern_spheretime"){
        locs <- cbind(lonlat,runif(n))
    } else if(covfun_name=="arma_matern_spheretime_warp"){
        locs <- cbind(lonlat,runif(n))
    } else {
        stop("unrecognized covariance in testing function")
    }
    return(locs)
}

test_that("covariance functions return positive definite matrix", {
        
    covfun_names <- c(
        "arma_exponential_isotropic",
        "arma_matern_isotropic",
        "arma_matern_spacetime", 
        "arma_matern_anisotropic2D",
        "arma_exponential_anisotropic3D",
        "arma_matern_anisotropic3D",
        "arma_matern_nonstat_var",
        "arma_matern_sphere",
        "arma_matern_sphere_warp",
        "arma_matern_spheretime_warp",
        "arma_matern_spheretime"
    )
    
    n <- 100    
    for(j in 1:length(covfun_names)){
        locs <- get_test_locs(covfun_names[j],n)
        covparms <- get_start_parms(rnorm(n),rep(1,n),locs,covfun_names[j])
        covfun <- get( covfun_names[j] )
        covmat <- covfun( covparms$start_parms, locs )
        cholmat <- t(chol(covmat))
        logdet <- 2*sum(log(diag(cholmat)))
        expect_lt( logdet, sum(log(diag(covmat))) )
    }
    
})



test_that("covariance function derivatives match finite differencing", {
        
    covfun_names <- c(
        "arma_exponential_isotropic",
        "arma_matern_isotropic",
        "arma_matern_spacetime", 
        "arma_matern_anisotropic2D",
        "arma_exponential_anisotropic3D",
        "arma_matern_anisotropic3D",
        "arma_matern_nonstat_var",
        "arma_matern_sphere",
        "arma_matern_sphere_warp",
        "arma_matern_spheretime_warp",
        "arma_matern_spheretime"
    )
    
    n <- 100    
    for(j in 1:length(covfun_names)){
        locs <- get_test_locs(covfun_names[j],n)
        covparms <- get_start_parms(rnorm(n),rep(1,n),locs,covfun_names[j])
        covparms <- covparms$start_parms
        nparms <- length(covparms)
        covfun <- get( covfun_names[j] )
        dcovfun <- get(paste0("d_",covfun_names[j]))
        covmat <- covfun( covparms, locs )
        dcovmat <- dcovfun( covparms,locs )
        ddcov <- array(NA, c(n,n,nparms) )
        eps <- 1e-8
        for(k in 1:nparms){
            dcovparms <- covparms
            dcovparms[k] <- covparms[k] + eps
            cov <- covfun( dcovparms, locs )
            ddcov[,,k] <- (cov - covmat)/eps
        }
        denom <- covmat[1,1]
        expect_equal( dcovmat/denom, ddcov/denom, tolerance = 1e-4 )
    }
})
