
#' Approximate GP simulation
#'
#' Calculates an approximation to the inverse Cholesky
#' factor of the covariance matrix using Vecchia's approximation,
#' then the simulation is produced by solving a linear system
#' with a vector of uncorrelated standard normals
#' @param m Number of nearest neighbors to use in approximation
#' @inheritParams vecchia_meanzero_loglik
#' @return vector of simulated values
#' @examples
#' locs <- as.matrix( expand.grid( (1:50)/50, (1:50)/50 ) )
#' y <- fast_Gp_sim(c(4,0.2,0.5,0), "matern_isotropic",  locs, 30 )
#' fields::image.plot( matrix(y,50,50) )
#' @export
fast_Gp_sim <- function( covparms, covfun_name = "matern_isotropic", locs, m = 30 ){
    
    # make sure locs is a matrix
    locs <- as.matrix(locs)
    
    # figure out if lonlat or not
    lonlat <- get_linkfun(covfun_name)$lonlat
    space_time <- get_linkfun(covfun_name)$space_time
    if(space_time){
        st_scale <- covparms[2:3]
    } else {
        st_scale <- NULL
    }

    n <- nrow(locs)
    m <- min(m,n-1)
    ord <- order_maxmin(locs,lonlat=lonlat,space_time=space_time)
    locsord <- locs[ord,,drop=FALSE]
    NNarray <- find_ordered_nn(locsord,m,lonlat=lonlat,st_scale=st_scale)
    Linv <- vecchia_Linv( covparms, covfun_name, locsord, NNarray )
    y <- fast_Gp_sim_Linv( Linv, NNarray )
    y[ord] <- y
    return(y)

}


#' Approximate GP simulation with specified Linverse
#'
#' In situations where we want to do many gaussian process
#' simulations from the same model, we can compute Linverse
#' once and reuse it, rather than recomputing for each identical simulation.
#' This function also allows the user to input the vector of standard normals \code{z}.
#' @param Linv Matrix containing the entries of Linverse, usually the output from
#' \code{vecchia_Linv}.
#' @param NNarray Matrix of nearest neighbor indices, usually the output from \code{\link{find_ordered_nn}}
#' @param z Optional vector of standard normals. If not specified,
#' these are computed within the function.
#' @return vector of simulated values
#' @examples
#' locs <- as.matrix( expand.grid( (1:50)/50, (1:50)/50 ) )
#' ord <- order_maxmin(locs)
#' locsord <- locs[ord,]
#' m <- 10
#' NNarray <- find_ordered_nn(locsord,m)
#' covparms <- c(2, 0.2, 1, 0)
#' Linv <- vecchia_Linv( covparms, "matern_isotropic", locsord, NNarray )
#' y <- fast_Gp_sim_Linv(Linv,NNarray)
#' y[ord] <- y
#' fields::image.plot( matrix(y,50,50) )
#' @export
fast_Gp_sim_Linv <- function( Linv, NNarray, z = NULL ){

    if( is.null(z) ){ z = stats::rnorm(nrow(Linv)) }
    y <- L_mult( Linv, z, NNarray )
    return(y)

}
