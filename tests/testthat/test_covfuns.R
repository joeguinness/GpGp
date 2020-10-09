
context("Covariance Functions")

covfun_names <- c(
    "matern_isotropic",
    "exponential_isotropic",
    "matern15_isotropic",
    "matern25_isotropic",
    "matern35_isotropic",
    "matern45_isotropic",
    "matern_scaledim", 
    "matern15_scaledim", 
    "matern25_scaledim", 
    "matern35_scaledim", 
    "matern45_scaledim", 
    "exponential_scaledim", 
    "matern_spacetime", 
    "exponential_spacetime", 
    "matern_anisotropic2D",
    "exponential_anisotropic2D",
    "matern_anisotropic3D",
    "matern_anisotropic3D_alt",
    "exponential_anisotropic3D",
    "matern_nonstat_var",
    "exponential_nonstat_var",
    "matern_sphere",
    "exponential_sphere",
    "matern_spheretime",
    "exponential_spheretime",
    "matern_sphere_warp",
    "exponential_sphere_warp",
    "matern_spheretime_warp",
    "exponential_spheretime_warp"
)

get_test_locs <- function(covfun_name,n){
    
    nside <- round(sqrt(n))
    longrid <- seq(1,360,length.out=nside)
    latgrid <- seq(-80,80,length.out=nside)
    lonlat <- as.matrix(expand.grid(longrid,latgrid))
    if(covfun_name=="exponential_isotropic"){
        locs <- matrix(runif(3*n),n,3)           
    } else if(covfun_name=="matern_isotropic"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="matern15_isotropic"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="matern25_isotropic"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="matern35_isotropic"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="matern45_isotropic"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="matern_scaledim"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="matern15_scaledim"){
        locs <- matrix(runif(3*n),n,3)
    } else if(covfun_name=="matern25_scaledim"){
        locs <- matrix(runif(3*n),n,3)
    } else if(covfun_name=="matern35_scaledim"){
        locs <- matrix(runif(3*n),n,3)
    } else if(covfun_name=="matern45_scaledim"){
        locs <- matrix(runif(3*n),n,3)
    } else if(covfun_name=="exponential_scaledim"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="matern_spacetime"){
        locs <- matrix(runif(3*n),n,3)
    } else if(covfun_name=="exponential_spacetime"){
        locs <- matrix(runif(3*n),n,3)
    } else if(covfun_name=="matern_anisotropic2D"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="exponential_anisotropic2D"){
        locs <- matrix(runif(2*n),n,2)
    } else if(covfun_name=="exponential_anisotropic3D"){
        locs <- matrix(runif(3*n),n,3)           
    } else if(covfun_name %in% c("matern_anisotropic3D","matern_anisotropic3D_alt") ){
        locs <- matrix(runif(3*n),n,3)           
    } else if(covfun_name=="matern_nonstat_var"){
        locs <- matrix(runif(6*n),n,6)           
    } else if(covfun_name=="exponential_nonstat_var"){
        locs <- matrix(runif(6*n),n,6)           
    } else if(covfun_name=="matern_sphere"){
        locs <- lonlat
    } else if(covfun_name=="exponential_sphere"){
        locs <- lonlat
    } else if(covfun_name=="matern_sphere_warp"){
        locs <- lonlat
    } else if(covfun_name=="exponential_sphere_warp"){
        locs <- lonlat
    } else if(covfun_name=="matern_spheretime"){
        locs <- cbind(lonlat,runif(n))
    } else if(covfun_name=="exponential_spheretime"){
        locs <- cbind(lonlat,runif(n))
    } else if(covfun_name=="matern_spheretime_warp"){
        locs <- cbind(lonlat,runif(n))
    } else if(covfun_name=="exponential_spheretime_warp"){
        locs <- cbind(lonlat,runif(n))
    } else {
        stop("unrecognized covariance in testing function")
    }
    return(locs)
}

test_that("covariance functions return positive definite matrix", {
        
    
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
