# a short vignette demonstrating how to use the functions
# library("GpGp")

# grid size for data locations
gsize <- 70
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

# covariance function and parameters
# covfun <- matern_isotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 0.5, nugget = 0.1)

# simulate some data
y <- fast_Gp_sim(covparms, "exponential_isotropic",locs,20)

# generate an ordering
ord <- order_maxmin(locs)

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]

# find the ordered m nearest neighbors
m <- 30
NNarray <- find_ordered_nn(locs,m)

# automatically group the observations
NNlist <- group_obs(NNarray, 2)
object.size(NNarray)
object.size(NNlist)
num_cond <- NNlist[["local_resp_inds"]]-1
mean(num_cond)

X <- matrix( rep(1,n), n, 1 )
Xord <- X[ord,,drop=FALSE]

# get ungrouped and grouped likelihood
ii <- c(1,2,4)
system.time(ll1 <- vecchia_meanzero_loglik(covparms[ii],"exponential_isotropic",yord,locsord,NNarray))
system.time(
  ll2 <- vecchia_profbeta_loglik(covparms[ii],"exponential_isotropic",yord,Xord,locsord,NNarray))
system.time( ll3 <-
               vecchia_profbeta_loglik_grad_info(covparms[ii],"exponential_isotropic",yord,Xord,locsord,NNarray))



system.time( ll1 <- vecchia_grouped_meanzero_loglik(covparms[c(1,2,4)],"exponential_isotropic",yord,locsord,NNlist ) )
system.time( ll2 <- vecchia_grouped_profbeta_loglik(covparms[c(1,2,4)],"exponential_isotropic",yord,Xord,locsord,NNlist ) )
system.time( ll3 <- vecchia_grouped_profbeta_loglik_grad_info(covparms[c(1,2,4)],"exponential_isotropic",yord,Xord,locsord,NNlist ) )

# get entries of L^{-1}
system.time( Linv1 <- vecchia_Linv(covparms[c(1,2,4)], "exponential_isotropic", locsord, NNarray) )

# do a multiplication
system.time( z1 <- Linv_mult(Linv1,yord,NNarray) )


# fit a model
t1 <- proc.time()
gpfit <- fit_model( y = y, locs = locs, X = X, "exponential_isotropic",
                    fixed_parms = 3, start_parms = c(1,1,0.0) 
)
print(proc.time() - t1)





