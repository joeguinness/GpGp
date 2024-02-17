


# functions to find nearest neighbors and do the grouping operations




#' Naive brute force nearest neighbor finder
#' 
#' @param locs matrix of locations
#' @param m number of neighbors
#' @return An matrix containing the indices of the neighbors. Row \code{i} of the
#' returned matrix contains the indices of the nearest \code{m}
#' locations to the \code{i}'th location. Indices are ordered within a
#' row to be increasing in distance. By convention, we consider a location
#' to neighbor itself, so the first entry of row \code{i} is \code{i}, the
#' second entry is the index of the nearest location, and so on. Because each
#' location neighbors itself, the returned matrix has \code{m+1} columns.
#' @export
find_ordered_nn_brute <- function( locs, m ){
     # find the m+1 nearest neighbors to locs[j,] in locs[1:j,]
     # by convention, this includes locs[j,], which is distance 0
     n <- dim(locs)[1]
     m <- min(m,n-1)
     NNarray <- matrix(NA,n,m+1)
     for(j in 1:n ){
         distvec <- c(fields::rdist(locs[1:j,,drop=FALSE],locs[j,,drop=FALSE]) )
         NNarray[j,1:min(m+1,j)] <- order(distvec)[1:min(m+1,j)]
     }
     return(NNarray)
}






#' Find ordered nearest neighbors.
#'
#' Given a matrix of locations, find the \code{m} nearest neighbors
#' to each location, subject to the neighbors coming
#' previously in the ordering. The algorithm uses the kdtree
#' algorithm in the FNN package, adapted to the setting
#' where the nearest neighbors must come from previous
#' in the ordering.
#'
#' @param m Number of neighbors to return
#' @inheritParams order_dist_to_point
#' @inheritParams order_maxmin
#' @param st_scale factor by which to scale the spatial and temporal coordinates
#' for distance calculations. The function assumes that the last column of
#' the locations is the temporal dimension, and the rest of the columns
#' are spatial dimensions. The spatial dimensions are divided by \code{st_scale[1]},
#' and the temporal dimension is divided by \code{st_scale[2]}, before distances are
#' calculated. If \code{st_scale} is \code{NULL}, no scaling is used. We 
#' recommend setting \code{st_scale} manually so that each observation gets
#' neighbors that hail multiple directions in space and time.
#' @return An matrix containing the indices of the neighbors. Row \code{i} of the
#' returned matrix contains the indices of the nearest \code{m}
#' locations to the \code{i}'th location. Indices are ordered within a
#' row to be increasing in distance. By convention, we consider a location
#' to neighbor itself, so the first entry of row \code{i} is \code{i}, the
#' second entry is the index of the nearest location, and so on. Because each
#' location neighbors itself, the returned matrix has \code{m+1} columns.
#' @examples
#' locs <- as.matrix( expand.grid( (1:40)/40, (1:40)/40 ) )     
#' ord <- order_maxmin(locs)        # calculate an ordering
#' locsord <- locs[ord,]            # reorder locations
#' m <- 20
#' NNarray <- find_ordered_nn(locsord,20)  # find ordered nearest 20 neighbors
#' ind <- 100
#' # plot all locations in gray, first ind locations in black,
#' # ind location with magenta circle, m neighhbors with blue circle
#' plot( locs[,1], locs[,2], pch = 16, col = "gray" )
#' points( locsord[1:ind,1], locsord[1:ind,2], pch = 16 )
#' points( locsord[ind,1], locsord[ind,2], col = "magenta", cex = 1.5 )
#' points( locsord[NNarray[ind,2:(m+1)],1], 
#'     locsord[NNarray[ind,2:(m+1)],2], col = "blue", cex = 1.5 )
#' @export
find_ordered_nn <- function(locs,m, lonlat = FALSE, st_scale = NULL){
    
    # convert locs to matrix (works on vectors and data frames)
    locs <- as.matrix(locs)

    # number of locations
    n <- nrow(locs)
    m <- min(m,n-1)
    mult <- 2
    
    # FNN::get.knnx has strange behavior for exact matches
    # so add a small amount of noise to each location
    ee <- min(apply( locs, 2, stats::sd ))
    locs <- locs + matrix( ee*1e-4*stats::rnorm(n*ncol(locs)), n, ncol(locs) )    
    
    if(lonlat){ # convert lonlattime to xyztime or lonlat to xyz
        lon <- locs[,1]
        lat <- locs[,2]
        lonrad <- lon*2*pi/360
        latrad <- (lat+90)*2*pi/360
        x <- sin(latrad)*cos(lonrad)
        y <- sin(latrad)*sin(lonrad)
        z <- cos(latrad)
        if(ncol(locs)==3){
            time <- locs[,3]
            locs <- cbind(x,y,z,time)
        } else {
            locs <- cbind(x,y,z)
        }
    }

    
    if( !is.null(st_scale) ){ 
        d <- ncol(locs)-1
        locs[ , 1:d] <- locs[ , 1:d]/st_scale[1]
        locs[ , d+1] <- locs[ , d+1]/st_scale[2]
    }

    # to store the nearest neighbor indices
    NNarray <- matrix(NA,n,m+1)

    # to the first mult*m+1 by brutce force
    maxval <- min( mult*m + 1, n )
    NNarray[1:maxval,] <- find_ordered_nn_brute(locs[1:maxval,,drop=FALSE],m)

    query_inds <- min( maxval+1, n):n
    data_inds <- 1:n

    msearch <- m

    while( length(query_inds) > 0 ){
        msearch <- min( max(query_inds), 2*msearch )
        data_inds <- 1:min( max(query_inds), n )
        NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
        less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
        sum_less_than_k <- apply(less_than_k,1,sum)
        ind_less_than_k <- which(sum_less_than_k >= m+1)

        NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))

        NNarray[ query_inds[ind_less_than_k], ] <- NN_m

        query_inds <- query_inds[-ind_less_than_k]

    }

    return(NNarray)
}










#' Automatic grouping (partitioning) of locations
#'
#' Take in an array of nearest neighbors, and automatically partition
#' the array into groups that share neighbors.
#' This is helpful to speed the computations and improve their accuracy.
#' The function returns a list, with each list element containing one or
#' several rows of NNarray. The algorithm attempts to find groupings such that
#' observations within a group share many common neighbors.
#' @param NNarray Matrix of nearest neighbor indices, usually the result of \code{\link{find_ordered_nn}}.
#' @param exponent Within the algorithm, two groups are merged if the number of unique
#' neighbors raised to the \code{exponent} power is less than the sum of the unique numbers
#' raised to the \code{exponent} power from the two groups.
#' @return A list with elements defining the grouping. The list entries are:
#' \itemize{
#'   \item \code{all_inds}: vector of all indices of all blocks.
#'   \item \code{last_ind_of_block}: The \code{i}th entry tells us the
#'       location in \code{all_inds} of the last index of the \code{i}th block. Thus the length
#'       of \code{last_ind_of_block} is the number of blocks, and \code{last_ind_of_block} can
#'       be used to chop \code{all_inds} up into blocks.
#'   \item \code{global_resp_inds}: The \code{i}th entry tells us the
#'       index of the \code{i}th response, as ordered in \code{all_inds}.
#'   \item \code{local_resp_inds}: The \code{i}th entry tells us the location within 
#'       the block of the response index.
#'   \item \code{last_resp_of_block}: The \code{i}th entry tells us the
#'       location within \code{local_resp_inds} and \code{global_resp_inds} of the last
#'       index of the \code{i}th block. \code{last_resp_of_block} is to 
#'       \code{global_resp_inds} and \code{local_resp_inds}
#'       as \code{last_ind_of_block} is to \code{all_inds}.
#' }
#' @examples
#' locs <- matrix( runif(200), 100, 2 )   # generate random locations
#' ord <- order_maxmin(locs)              # calculate an ordering
#' locsord <- locs[ord,]                  # reorder locations
#' m <- 10
#' NNarray <- find_ordered_nn(locsord,m)  # m nearest neighbor indices
#' NNlist2 <- group_obs(NNarray)          # join blocks if joining reduces squares
#' NNlist3 <- group_obs(NNarray,3)        # join blocks if joining reduces cubes
#' object.size(NNarray)
#' object.size(NNlist2)
#' object.size(NNlist3)
#' mean( NNlist2[["local_resp_inds"]] - 1 )   # average number of neighbors (exponent 2)
#' mean( NNlist3[["local_resp_inds"]] - 1 )   # average number of neighbors (exponent 3)
#' 
#' all_inds <- NNlist2$all_inds
#' last_ind_of_block <- NNlist2$last_ind_of_block
#' inds_of_block_2 <- all_inds[ (last_ind_of_block[1] + 1):last_ind_of_block[2] ]
#' 
#' local_resp_inds <- NNlist2$local_resp_inds
#' global_resp_inds <- NNlist2$global_resp_inds
#' last_resp_of_block <- NNlist2$last_resp_of_block
#' local_resp_of_block_2 <- 
#'     local_resp_inds[(last_resp_of_block[1]+1):last_resp_of_block[2]]
#'
#' global_resp_of_block_2 <- 
#'     global_resp_inds[(last_resp_of_block[1]+1):last_resp_of_block[2]]
#' inds_of_block_2[local_resp_of_block_2]
#' # these last two should be the same
#' 
#' @export
group_obs <- function(NNarray, exponent = 2){
    n <- nrow(NNarray)
    m <- ncol(NNarray)-1

    clust <- vector("list",n)
    for(j in 1:n) clust[[j]] <- j
    for( ell in 2:(m+1) ){  # 2:(m+1)?
        sv <- which( NNarray[,1] - NNarray[,ell] < n )
        for(j in sv){
            k <- NNarray[j,ell]
            if( length(clust[[k]]) > 0){
                nunique <- length(unique(c(NNarray[c(clust[[j]],clust[[k]]),])))

                # this is the rule for deciding whether two groups
                # should be combined
                if( nunique^exponent <= length(unique(c(NNarray[clust[[j]],])))^exponent + length(unique(c(NNarray[clust[[k]],])))^exponent ) {
                    clust[[j]] <- c(clust[[j]],clust[[k]])
                    clust[[k]] <- numeric(0)
                }
            }
        }
    }
    zeroinds <- unlist(lapply(clust,length)==0)
    clust[zeroinds] <- NULL
    nb <- length(clust)
    all_inds <- rep(NA,n*(m+1))
    last_ind_of_block <- rep(NA,length(clust))
    last_resp_of_block <- rep(NA,length(clust))
    #num_resp_inblock <- rep(NA,length(clust))
    local_resp_inds <- rep(NA,n)
    global_resp_inds <- rep(NA,n)

    last_ind <- 0
    last_resp <- 0
    for(j in 1:nb){
        resp <- sort(clust[[j]])
        last_resp_of_block[j] <- last_resp + length(resp)
        global_resp_inds[(last_resp+1):last_resp_of_block[j]] <- resp

        inds_inblock <- sort( unique( c(NNarray[resp,]) ) )
        last_ind_of_block[j] <- last_ind + length(inds_inblock)
        all_inds[(last_ind+1):last_ind_of_block[j]] <- inds_inblock
        last_ind <- last_ind_of_block[j]

        local_resp_inds[(last_resp+1):last_resp_of_block[j]] <- which( inds_inblock %in% resp )
        last_resp <- last_resp_of_block[j]
    }
    all_inds <- all_inds[ !is.na(all_inds) ]
    NNlist <- list( all_inds = all_inds, last_ind_of_block = last_ind_of_block,
                    global_resp_inds = global_resp_inds, local_resp_inds = local_resp_inds,
                    last_resp_of_block = last_resp_of_block )
    return(NNlist)
}
