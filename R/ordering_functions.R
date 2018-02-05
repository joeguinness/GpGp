


#' Distance to specified point ordering
#' 
#' Return the ordering of locations increasing in their
#' distance to some specified location
#' 
#' @param locs A matrix of locations. Each row of \code{locs} contains a location, which can
#' be a point in Euclidean space R^d, a point in space-time R^d x T, 
#' a longitude and latitude (in degrees) giving a point on the sphere, 
#' or a longitude, latitude, and time giving a point in the sphere-time domain.
#' @param loc0 A vector containing a single location in R^d.
#' @param lonlat TRUE/FALSE whether locations are longitudes and latitudes.
#' @return A vector of indices giving the ordering, i.e. 
#' the first element of this vector is the index of the location nearest to \code{loc0}.
#' @examples
#' n <- 100             # Number of locations
#' d <- 2               # dimension of domain
#' locs <- matrix( runif(n*d), n, d )
#' loc0 <- c(1/2,1/2)
#' ord <- order_dist_to_point(locs,loc0)
#' @export
order_dist_to_point <- function( locs, loc0, lonlat = FALSE ){
    
    if(lonlat){
        lon <- locs[,1]
        lat <- locs[,2]
        lonrad <- lon*2*pi/360
        latrad <- (lat+90)*2*pi/360
        x <- sin(latrad)*cos(lonrad)
        y <- sin(latrad)*sin(lonrad)
        z <- cos(latrad)
        locs <- cbind(x,y,z)
    }

    
    if (!requireNamespace("fields", quietly = TRUE)) {
        stop("fields package required for this function. Please install it.",
             call. = FALSE)
    }
    d <- ncol(locs)         # input dimension
    if(d != length(loc0)){
        stop("location in loc0 not in the same domain as the locations in locs")
    }    
    loc0 <- matrix(c(loc0),1,d)             # order by distance to loc0
    distvec <- fields::rdist(locs,loc0)
    orderinds <- order(distvec)
    return(orderinds)
}

#' Middle-out ordering
#' 
#' Return the ordering of locations increasing in their
#' distance to the average location 
#' 
#' @inheritParams order_dist_to_point
#' 
#' @return A vector of indices giving the ordering, i.e. 
#' the first element of this vector is the index of the location nearest the center.
#' @examples
#' n <- 100             # Number of locations
#' d <- 2               # dimension of domain
#' locs <- matrix( runif(n*d), n, d )
#' ord <- order_middleout(locs)
#' @export
order_middleout <- function( locs, lonlat = FALSE ){
    
    if(lonlat){
        lon <- locs[,1]
        lat <- locs[,2]
        lonrad <- lon*2*pi/360
        latrad <- (lat+90)*2*pi/360
        x <- sin(latrad)*cos(lonrad)
        y <- sin(latrad)*sin(lonrad)
        z <- cos(latrad)
        locs <- cbind(x,y,z)
    }

    d <- ncol(locs)
    loc0 <- matrix(colMeans(locs),1,d)
    orderinds <- order_dist_to_point(locs,loc0)
    return(orderinds)
}


#' Sorted coordinate ordering
#' 
#' Return the ordering of locations sorted along one of the
#' coordinates or the sum of multiple coordinates
#' 
#' @param coordinate integer or vector of integers in {1,...,d}. If a single integer,
#' coordinates are ordered along that coordinate. If multiple integers,
#' coordinates are ordered according to the sum of specified coordinate values. For example,
#' when \code{d=2}, \code{coordinate = c(1,2)} orders from bottom left to top right.
#' @inheritParams order_dist_to_point
#' @return A vector of indices giving the ordering, i.e. 
#' the first element of this vector is the index of the first location.
#' @examples
#' n <- 100             # Number of locations
#' d <- 2               # dimension of domain
#' locs <- matrix( runif(n*d), n, d )
#' ord1 <- order_coordinate(locs, 1 )
#' ord12 <- order_coordinate(locs, c(1,2) )
#' @export
order_coordinate <- function( locs, coordinate ){
    order(rowSums(locs[,coordinate,drop=FALSE]))
}



#' Maximum minimum distance ordering
#' 
#' Return the indices of an approximation to the maximum minimum distance ordering.
#' A point in the center is chosen first, and then each successive point
#' is chosen to maximize the minimum distance to previously selected points
#'  
#' @inheritParams order_dist_to_point
#' @param space_time TRUE if locations are euclidean space-time locations, 
#' FALSE otherwise. If set to TRUE, temporal dimension is ignored. 
#' If set to FALSE, temporal dimension treated as another spatial dimension (not recommended).
#' @return A vector of indices giving the ordering, i.e. 
#' the first element of this vector is the index of the first location.
#' @examples
#' # planar coordinates
#' nvec <- c(50,50)
#' locs <- as.matrix( expand.grid( 1:nvec[1]/nvec[1], 1:nvec[2]/nvec[2] ) )
#' ord <- order_maxmin(locs)
#' par(mfrow=c(1,3))
#' plot( locs[ord[1:100],1], locs[ord[1:100],2], xlim = c(0,1), ylim = c(0,1) )
#' plot( locs[ord[1:300],1], locs[ord[1:300],2], xlim = c(0,1), ylim = c(0,1) )
#' plot( locs[ord[1:900],1], locs[ord[1:900],2], xlim = c(0,1), ylim = c(0,1) )
#' 
#' # longitude/latitude coordinates (sphere)
#' latvals <- seq(-80, 80, length.out = 40 )
#' lonvals <- seq( 0, 360, length.out = 81 )[1:80]
#' locs <- as.matrix( expand.grid( lonvals, latvals ) )
#' ord <- order_maxmin(locs, lonlat = TRUE)
#' par(mfrow=c(1,3))
#' plot( locs[ord[1:100],1], locs[ord[1:100],2], xlim = c(0,360), ylim = c(-90,90) )
#' plot( locs[ord[1:300],1], locs[ord[1:300],2], xlim = c(0,360), ylim = c(-90,90) )
#' plot( locs[ord[1:900],1], locs[ord[1:900],2], xlim = c(0,360), ylim = c(-90,90) )
#' 
#' @export
order_maxmin <- function(locs, lonlat = FALSE, space_time = FALSE){

    if(lonlat){
        lon <- locs[,1]
        lat <- locs[,2]
        lonrad <- lon*2*pi/360
        latrad <- (lat+90)*2*pi/360
        x <- sin(latrad)*cos(lonrad)
        y <- sin(latrad)*sin(lonrad)
        z <- cos(latrad)
        locs <- cbind(x,y,z)
    }
    
    if(space_time){ # simply ignore the temporal dimension
        # this will give roughly random order in time
        d <- ncol(locs)-1
        locs <- locs[,1:d,drop=FALSE]
    }
    
    # get number of locs
    n <- nrow(locs)
    m <- round(sqrt(n))
    # m is number of neighbors to search over
    # get the past and future nearest neighbors
    NNall <- FNN::get.knn( locs, k = m )$nn.index
    # pick a random ordering
    index_in_position <- c( sample(n), rep(NA,1*n) )
    position_of_index <- order(index_in_position[1:n])
    # loop over the first n/4 locations
    # move an index to the end if it is a
    # near neighbor of a previous location
    curlen <- n
    #curpos <- 1
    nmoved <- 0
    for(j in 2:(2*n) ){
        nneigh <- round( min(m,n/(j-nmoved+1)) )
        neighbors <- NNall[index_in_position[j],1:nneigh]
        if( min( position_of_index[neighbors], na.rm = TRUE ) < j ){
            nmoved <- nmoved+1
            curlen <- curlen + 1
            position_of_index[ index_in_position[j] ] <- curlen
            index_in_position[curlen] <- index_in_position[j]
            index_in_position[j] <- NA
        }
    }
    ord <- index_in_position[ !is.na( index_in_position ) ]

    return(ord)
}
