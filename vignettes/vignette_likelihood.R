

# a short vignette demonstrating how to use the functions
library("GpGp")

# grid size for data locations
gsize <- 100
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

# covariance function and parameters
# covfun <- maternIsotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 1.0, nugget = 0)

# simulate some data
y <- fast_Gp_sim(covparms, "matern_isotropic",locs,40)

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

if( n < 6000 ){ # only do this if we can store the covariance matrix
    covmat <- matern_isotropic(covparms,locsord)
    cholmat <- t(chol(covmat))
    z3 <- forwardsolve(cholmat,yord)
    print( sd(z3-z1) )
    print( sd(z3-z2) )
}





