
# analysis of argo data
devtools::load_all()
data("argo2016")
names(argo2016)
attach(argo2016)
n <- length(lon)

# plot of data locations and times
par(mfrow=c(1,2))
plot(lon,lat,pch=1,cex=0.3)
hist(day,breaks=60)

# plot the data
par(mfrow=c(1,1))
fields::quilt.plot(lon,lat,temp100,nx=200,ny=200)

# define locations
lonlat <- cbind(lon,lat)
lonlattime <- cbind(lonlat,day)

# define covariates (quadratic in latitude)
X <- cbind(rep(1,n),lat,lat^2)

# subset of data for faster results
inds <- sample(n,10000)
# full analysis
#inds <- 1:n

# use the 100-bar temperature data
# also available: 150-bar and 200-bar
temp <- temp100

# store timing results
timing <- rep(NA,4)

# spatial isotropic
t1 <- proc.time()
fit1 <- fit_model(y = temp[inds], X = X[inds,], locs = lonlat[inds,], 
    covfun_name = "exponential_sphere", group = FALSE )
timing[1] <- (proc.time() - t1)[3]


# spacetime, isotropic in each
t1 <- proc.time()
fit2 <- fit_model(y = temp[inds], X = X[inds,], locs = lonlattime[inds,], 
    covfun_name = "exponential_spheretime", group = FALSE, st_scale = c(0.2,16) )
timing[2] <- (proc.time() - t1)[3]


# spatial warping
t1 <- proc.time()
fit3 <- fit_model(y = temp[inds], X = X[inds,], locs = lonlat[inds,], 
    covfun_name = "exponential_sphere_warp", group = FALSE )
timing[3] <- (proc.time() - t1)[3]


# spacetime, warping in space
t1 <- proc.time()
fit4 <- fit_model(y = temp[inds], X = X[inds,], locs = lonlattime[inds,], 
    covfun_name = "exponential_spheretime_warp", group = FALSE, st_scale = c(0.2,16) )
timing[4] <- (proc.time() - t1)[3]

# prediction locations and design matrix
mediantime <- median(day)
latgrid <- seq( min(lat), max(lat), length.out = 60 )
longrid <- seq( 0, 360, length.out = 121)[1:120] # so no locations repeated
locs_pred <- as.matrix( expand.grid(longrid,latgrid) )
n_pred <- nrow(locs_pred)
locstime_pred <- cbind( locs_pred, rep(mediantime, n_pred) )
X_pred <- as.matrix( cbind( rep(1,n_pred), locs_pred[,2], locs_pred[,2]^2 ) )

# predictions: this one uses only data in fit4 object (subset of all data)
pred <- predictions( fit4, locs_pred = locstime_pred, X_pred = X_pred )

# this one uses all data: overrides y_obs = fit$y, locs_obs = fit$locs, etc.
pred <- predictions( fit4, locs_pred = locstime_pred, X_pred = X_pred,
    y_obs = temp, locs_obs = lonlattime, X_obs = X )

# conditional simulations: only uses data in fit4 object
sim <- cond_sim( fit4, locs_pred = locstime_pred, X_pred = X_pred )

# uses all data
sim <- cond_sim( fit4, locs_pred = locstime_pred, X_pred = X_pred,
    y_obs = temp, locs_obs = lonlattime, X_obs = X )

# plot predictions and cond sim
par(mar=c(4,4,1,1),mfrow=c(1,2))
pred_array <- array( pred, c(length(longrid),length(latgrid)) )
fields::image.plot(longrid,latgrid,pred_array)
# plot conditional simulations
sim_array <- array( sim, c(length(longrid),length(latgrid)) )
fields::image.plot(longrid,latgrid,sim_array)


