
#' Approximate GP simulation
#'
#' Calculates an approximation to the inverse Cholesky
#' factor of the covariance matrix using Vecchia's approximation,
#' then the simulation is produced by solving a linear system
#' with a vector of uncorrelated standard normals
#' @param covfunNames string name of covariance functions. Currently supported
#' covariance functions are "maternIsotropic", "maternSphere", and "maternSphereTime"
#' @param covparms Vector of covariance function parameters. For "maternIsotropic" and
#' "maternSphere", these are (variance, range, smoothness, nugget). For "maternSphereTime",
#' these are (variance, spatial range, temporal range, smoothness, nugget). See documentation
#' of the covariance functions for more details.
#' @param locs Matrix of locations. Each row is an individual location. For "maternSphere",
#' locations are (lon,lat), where lon is in (-180,180), and lat is in (-90,90). For "maternSphereTime"
#' locations are (lon,lat,time).
#' @param m Number of nearest neighbors to use in approximation
#' @return vector of simulated values
#' @examples
#' locs <- as.matrix( expand.grid( (1:100)/100, (1:100)/100 ) )
#' y <- fast_Gp_sim(c(4,0.2,0.5,0), "matern_isotropic",  locs, 30 )
#' fields::image.plot( matrix(y,100,100) )
#' @export
fast_Gp_sim <- function( covparms, covfun_name = "matern_isotropic", locs, m = 30 ){

    ord <- order_maxmin(locs)
    locsord <- locs[ord,]
    NNarray <- find_ordered_nn(locsord,m)
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
#' \code{vecchiaLinv}.
#' @param NNarray Matrix of nearest neighbor indices, usually the output from findOrderedNN
#' @param z Optional vector of standard normals. If not specified,
#' these are computed within the function.
#' @return vector of simulated values
#' @examples
#' locs <- as.matrix( expand.grid( (1:100)/100, (1:100)/100 ) )
#' ord <- order_maxmin(locs)
#' locsord <- locs[ord,]
#' m = 3
#' NNarray <- find_ordered_nn(locsord,m)
#' covparms <- c(2, 0.2, 1, 0)
#' Linv <- vecchia_Linv( covparms, "matern_isotropic", locsord, NNarray )
#' y <- fast_Gp_sim_Linv(Linv,NNarray)
#' y[ord] <- y
#' fields::image.plot( matrix(y,100,100) )
#' @export
fast_Gp_sim_Linv <- function( Linv, NNarray, z = NULL ){

    if( is.null(z) ){ z = rnorm(nrow(Linv)) }
    y <- L_mult( Linv, z, NNarray )
    return(y)

}
