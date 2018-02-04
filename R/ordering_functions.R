


#' Distance to specified point ordering
#' 
#' Return the ordering of locations increasing in their
#' distance to some specified location
#' 
#' @param locs A matrix of locations in R^d. Each row of \code{locs} contains a location.
#' @param loc0 A vector containing a single location in R^d.
#' @return A vector of indices giving the ordering, i.e. 
#' the first element of this vector is the index of the location nearest to \code{loc0}.
#' @examples
#' n <- 100             # Number of locations
#' d <- 2               # dimension of domain
#' locs <- matrix( runif(n*d), n, d )
#' loc0 <- c(1/2,1/2)
#' ord <- order_dist_to_point(locs,loc0)
#' @export
order_dist_to_point <- function( locs, loc0 ){

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
#' @param locs A matrix of locations in R^d. Each row of \code{locs} contains a location.
#' @return A vector of indices giving the ordering, i.e. 
#' the first element of this vector is the index of the location nearest the center.
#' @examples
#' n <- 100             # Number of locations
#' d <- 2               # dimension of domain
#' locs <- matrix( runif(n*d), n, d )
#' ord <- order_middleout(locs)
#' @export
order_middleout <- function( locs ){
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
#' @param locs A matrix of locations in R^d. Each row of \code{locs} contains a location.
#' @return A vector of indices giving the ordering, 
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
#' @param locs A matrix of locations in R^d. Each row of \code{locs} contains a location.
#' @param lonlat Flag indicating whether the \code{locs} are longitudes and latitudes.
#' @return A vector of indices giving the ordering, 
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
order_maxmin <- function(locs, lonlat = FALSE){

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
    curpos <- 1
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
