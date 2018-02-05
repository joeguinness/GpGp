

# a short vignette demonstrating how to use the functions
library("GpGp")

# grid size for data locations
gsize <- 50
times <- 10
nvec <- c(gsize,gsize,times)
n <- prod(nvec)

# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locstime <- as.matrix(expand.grid(x1,x2,1:times))

# covariance function and parameters
# covfun <- maternIsotropic
covparms <- c(variance = 4, sp_range = 0.5, t_range = 3, smoothness = 1/2, nugget = 0)

# simulate some data
y <- fast_Gp_sim(covparms, "matern_space_time",locs = locstime,m = 100)
y <- fast_Gp_sim(covparms[c(1,2,4,5)], "matern_isotropic",locs = locstime,m = 100)
y_array <- array( y, nvec )
par(mfrow=c(1,3))
image(y_array[,,1], col = viridis::viridis(64))
image(y_array[,,2], col = viridis::viridis(64))
image(y_array[,,3], col = viridis::viridis(64))


# generate an ordering
ord <- order_maxmin(locs)

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]

# find the ordered m nearest neighbors
m <- 30
NNarray <- find_ordered_nn(locsord,m)

# automatically group the observations
NNlist <- group_obs(NNarray, 2)
object.size(NNarray)
object.size(NNlist)
num_cond <- NNlist[["local_resp_inds"]]-1
mean(num_cond)


# get ungrouped and grouped likelihood
system.time( ll1 <- vecchia_loglik(covparms, "matern_isotropic", yord, locsord, NNarray) )
system.time( ll2 <- vecchia_loglik_grouped(covparms, "matern_isotropic", yord, locsord, NNlist) )

# get entries of L^{-1}
system.time( Linv1 <- vecchia_Linv(covparms, "matern_isotropic", locsord, NNarray) )
system.time( Linv2 <- vecchia_Linv_grouped(covparms, "matern_isotropic", locsord, NNlist) )

# do a multiplication
system.time( z1 <- Linv_mult(Linv1,yord,NNarray) )
system.time( z2 <- Linv_mult_grouped(Linv2,yord,NNlist) )

if( n < 8000 ){ # only do this if we can store the covariance matrix
    covmat <- matern_isotropic(covparms,locsord)
    cholmat <- t(chol(covmat))
    z3 <- forwardsolve(cholmat,yord)
    print( sd(z3-z1) )
    print( sd(z3-z2) )
}
