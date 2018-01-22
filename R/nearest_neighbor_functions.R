


# functions to find nearest neighbors and do the grouping operations




# naive brute force nearest neighbor finder
find_ordered_nn_brute <- function( locs, m ){
     # find the m+1 nearest neighbors to locs[j,] in locs[1:j,]
     # by convention, this includes locs[j,], which is distance 0
     n <- dim(locs)[1]
     NNarray <- matrix(NA,n,m+1)
     for(j in 1:n ){
         distvec <- c(fields::rdist(locs[1:j,,drop=FALSE],locs[j,,drop=FALSE]) )
         NNarray[j,1:min(m+1,j)] <- order(distvec)[1:min(m+1,j)]
     }
     return(NNarray)
}






#' Find ordered nearest neighhbors
#' 
#' Given a matrix of reordered locations, find the \code{m} nearest neighbors
#' to each location, subject to the neighbors coming 
#' previously in the ordering. The algorithm uses the kdtree
#' algorithm in the FNN package, adapted to the setting
#' where the nearest neighbors must come from previous
#' in the ordering.
#'
#' @param locs Matrix containing ordered locations. 
#' Each row is contains a single location
#' @param m Nuber of neighbors to return
#' @return An matrix containing the indices of the neighbors. Row \code{i} of the
#' returned matrix contains the indices of the nearest \code{m} 
#' locations to the \code{i}'th location. Indices are ordered within a 
#' row to be increasing in distance. By convention, consider a location
#' to neighbor itself, so the first entry of row \code{i} is \code{i}, the 
#' second neighbor is the nearest other location, and so on. Because each
#' location neighbors itself, the returned matrix has \code{m+1} columns.
#' @examples
#' locs <- as.matrix( expand.grid( (1:40)/40, (1:40)/40 ) )     # grid of locations
#' ord <- order_maxmin(locs)                                     # calculate an ordering
#' locsord <- locs[ord,]                                        # reorder locations
#' m <- 20                               
#' NNarray <- find_ordered_nn(locsord,20)             # find ordered nearest 20 neighbors
#' ind <- 100
#' # plot all locations in gray, first ind locations in black,
#' # ind location with magenta circle, m neighhbors with blue circle
#' plot( locs[,1], locs[,2], pch = 16, col = "gray" )
#' points( locsord[1:ind,1], locsord[1:ind,2], pch = 16 )
#' points( locsord[ind,1], locsord[ind,2], col = "magenta", cex = 1.5 )
#' points( locsord[NNarray[ind,2:(m+1)],1], locsord[NNarray[ind,2:(m+1)],2], col = "blue", cex = 1.5 )
#' @export
find_ordered_nn <- function(locs,m){

    # number of locations
    n <- nrow(locs)
    mult <- 2

    # to store the nearest neighbor indices
    NNarray <- matrix(NA,n,m+1)
    
    # to the first mult*m+1 by brutce force
    NNarray[1:(mult*m+1),] <- find_ordered_nn_brute(locs[1:(mult*m+1),],m)

    query_inds <- (mult*m+2):n
    data_inds <- 1:n

    msearch <- m

    while( length(query_inds) > 0 ){

        msearch <- min( max(query_inds), 2*msearch )
        data_inds <- 1:max(query_inds)
        NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
        less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
        sum_less_than_k <- apply(less_than_k,1,sum)
        ind_less_than_k <- which(sum_less_than_k >= m+1)
        NN_less_than_k <- NN[ind_less_than_k,]

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
#' @param NNarray Matrix of nearest neighbor indices, usually the result of findOrderedNN
#' @param exponent Within the algorithm, two groups are merged if the number of unique
#' neighbors raised to the \code{exponent} power is less than the sum of the unique numbers 
#' raised to the \code{exponent} power from the two groups.
#' @return A list defining the partition of NNarray. Each list element contains the rows
#' of NNarray corresponding to the grouping, and the collection of all of the elements
#' is exhaustive, meaning that the entirety of NNarray can be reconstructed from the list.
#' @examples 
#' locs <- matrix( runif(200), 100, 2 )                         # generate random locations
#' ord <- order_maxmin(locs)                                     # calculate an ordering
#' locsord <- locs[ord,]                                        # reorder locations
#' m <- 10                               
#' NNarray <- find_ordered_nn(locsord,m)  
#' NNlist2 <- group_obs(NNarray)
#' NNlist3 <- group_obs(NNarray,3)
#' NNlist2[1:10]
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
                if( nunique^2 <= length(unique(c(NNarray[clust[[j]],])))^2 + length(unique(c(NNarray[clust[[k]],])))^2 ) {
                    clust[[j]] <- c(clust[[j]],clust[[k]])
                    clust[[k]] <- numeric(0)
                }
            }
        }
    }
    zeroinds <- unlist(lapply(clust,length)==0)
    clust[zeroinds] <- NULL
    NNlist <- lapply(clust,function(inds) NNarray[inds,,drop=FALSE])
    return(NNlist)
}



