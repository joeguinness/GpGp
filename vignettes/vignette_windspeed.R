
# analyze the averaged jason3 wind speed data
#library("GpGp")
devtools::load_all(".")


# loads in object
load("datasets/jason3_windspeed_avg.RData")
head(dframe)
attach(dframe)
n <- length(windspeed)

# quick plot
fields::quilt.plot(lon,lat,windspeed)

# higher resolution to see satellite path
fields::quilt.plot(lon,lat,windspeed,nx=400,ny=200)

# longitude latitude locations
locs <- cbind(lon,lat)
locstime <- cbind(locs,time)

# design matrix for constant mean
X <- as.matrix( rep(1,length(windspeed)) )

# get ordering and reorder 
ord <- order_maxmin(locs, lonlat = TRUE)
plot(locs[ord[1:500],1],locs[ord[1:500],2],pch=16,cex=1/2)

# reorder
locsord <- locs[ord,]
locstimeord <- locstime[ord,]
windspeedord <- windspeed[ord]
Xord <- as.matrix( X[ord,] )

# get nearest neighbors
NNarray <- find_ordered_nn(locsord,m=20, lonlat = TRUE)
NNlist <- group_obs(NNarray)

## first analysis: ignore temporal dimension using matern_sphere covariance
funtomax1 <- function( logparms ){
    parms <- exp(logparms)
    parms[5] <- logparms[5]
    loglik <- vecchia_loglik(parms[1:4],"matern_sphere",windspeedord-parms[5],
                                  cbind(lon[ord],lat[ord]),NNarray)
    return(-loglik)
}
startparms <- c(10,0.1,0.8,0.001,12)
fit1 <- optim(c(log(startparms[1:4]),startparms[5]),funtomax1,control=list(trace=5,maxit=500))
fit1$value


## second analysis: space-time covariance using maternSphereTime
funtomax2 <- function( logparms ){
    parms <- exp(logparms)
    parms[6] <- logparms[6]
    loglik <- vecchia_loglik(parms[1:5],"matern_sphere_time",windspeedord-parms[6],
                                  cbind(lon[ord],lat[ord],time[ord]),NNarray)
    return(-loglik)
}
startparms <- c(10,0.1,6e4,0.8,0.001,12)
fit2 <- optim(c(log(startparms[1:5]),startparms[6]),funtomax2,control=list(trace=5,maxit=500))
fit2$value


# fit3 and fit4 are the same as above, except
# profile out mean and variance parameters for faster fitting

# ignoring temporal dimension
funtomax3 <- function( logparms ){
    parms <- exp(logparms)
    loglik <- proflik_mean_variance(parms,"matern_sphere",windspeedord,
                                  Xord,cbind(lon[ord],lat[ord]),NNarray)
    return(-loglik)
}
startparms <- c(0.2,0.8,0.001)
fit3 <- optim(log(startparms),funtomax3,control=list(trace=5,maxit=500))
fit3$value

# with temporal dimension
funtomax4 <- function( logparms ){
    parms <- exp(logparms)
    parms[3] = min(parms[3],4)
    loglik <- proflik_mean_variance(parms,"matern_sphere_time",windspeedord,
                                    Xord,cbind(lon[ord],lat[ord],time[ord]),NNarray)
    return(-loglik)
}
startparms <- c(0.1,6e4,0.8,0.001)
system.time( fit4 <- optim(log(startparms),funtomax4,control=list(trace=5,maxit=50), method = "BFGS") )
fit4$value
fitted_model <- proflik_mean_variance(exp(fit4$par),"matern_sphere_time",windspeedord,
                               Xord,cbind(lon[ord],lat[ord],time[ord]),NNarray,return_parms = TRUE)
fitted_model



#               #
#  Predictions  #
#               #

# get prediction locations
mediantime <- median(time)
latgrid <- seq( min(lat), max(lat), length.out = 60 )
longrid <- seq( 0, 360, length.out = 121)[1:120] # so no locations repeated
locs_pred <- as.matrix( expand.grid(longrid,latgrid) )
n_pred <- nrow(locs_pred)
locstime_pred <- cbind( locs_pred, rep(mediantime, n_pred) )

# reorder prediction locations
ord_pred <- order_maxmin( locs_pred, lonlat = TRUE )
locstimeord_pred <- locstime_pred[ord_pred,]
Xord_pred <- as.matrix( rep(1,n_pred) )

condexpord <- predictions(fitted_model$covparms, "matern_sphere_time", 
    locstimeord, locstimeord_pred, Xord, Xord_pred, fitted_model$beta, windspeedord, 
    m = 60, lonlat = TRUE )

# undo the ordering
condexp <- rep(NA,length(condexpord))
condexp[ord_pred] <- condexpord
condexp_array <- array( condexp, c(length(longrid),length(latgrid)) )
fields::image.plot(condexp_array)


