
# analysis of argo data
devtools::load_all()
data("argo2016Jan2Mar")
names(argo2016Jan2Mar)
attach(argo2016Jan2Mar)
n <- length(lon)

# plot of data locations and times
par(mfrow=c(1,2))
plot(lon,lat,pch=1,cex=0.3)
hist(day,breaks=60)

# interpolate temperature to chosen pressure level
chosen_pres <- 100
temp <- rep(NA,n)
for(j in 1:n){
    # get first pressure larger than chosen_pres
    up <- which( c(pres_prof[[j]][[1]]) - chosen_pres >= 0 )[1]
    lo <- up-1
    # assumes that chosen pressure is between two pressures in dataset
    if(is.na(up) || up==1){ 
        warning(paste0("pressures in observation ",j," do not surround chosen_pres"))
    } else {
        pres_bounds <- c(pres_prof[[j]][[1]])[c(lo,up)]
        temp_bounds <- c(temp_prof[[j]][[1]])[c(lo,up)]
        temp[j] <- temp_bounds[1] + 
            diff(temp_bounds)*(chosen_pres-pres_bounds[1])/diff(pres_bounds)
    }
}

# filter to floats with temperature measurements at pressures surrounding
# chosen_pres
notna <- !is.na(temp)
lon1 <- lon[notna]
lat1 <- lat[notna]
temp1 <- temp[notna]
day1 <- day[notna]
n1 <- length(lon1)

# plot the data
fields::quilt.plot(lon1,lat1,temp1,nx=200,ny=200)

lonlat1 <- cbind(lon1,lat1)
lonlattime1 <- cbind(lonlat1,day1)

X1 <- cbind(rep(1,n1),lat1,lat1^2)

xyz <- array(NA, c(nrow(lonlat1), 3) )
lonrad <- lon1*2*pi/360
latrad <- (lat1+90)*2*pi/360
xyz[,1] <- sin(latrad)*cos(lonrad)
xyz[,2] <- sin(latrad)*sin(lonrad)
xyz[,3] <- cos(latrad)



# subset of data for faster cran checks
inds <- sample(n1,2000)

#NNarray <- find_ordered_nn(xyz,m=30)
timing <- rep(NA,4)

# spatial isotropic
t1 <- proc.time()
fit1 <- fit_model(y = temp1[inds], X = X1[inds,], locs = lonlat1[inds,], 
    covfun_name = "arma_matern_sphere", group = FALSE)
timing[1] <- (proc.time() - t1)[3]


# spacetime, isotropic in each
t1 <- proc.time()
fit2 <- fit_model(y = temp1[inds], X = X1[inds,], locs = lonlattime1[inds,], 
    covfun_name = "arma_matern_spheretime", group = FALSE )
timing[2] <- (proc.time() - t1)[3]


# spatial warping
t1 <- proc.time()
fit3 <- fit_model(y = temp1[inds], X = X1[inds,], locs = lonlat1[inds,], 
    covfun_name = "arma_matern_sphere_warp", group = FALSE )
timing[3] <- (proc.time() - t1)[3]


# spacetime, warping in space
t1 <- proc.time()
fit4 <- fit_model(y = temp1[inds], X = X1[inds,], locs = lonlattime1[inds,], 
    covfun_name = "arma_matern_spheretime_warp", group = FALSE )
timing[4] <- (proc.time() - t1)[3]

# prediction locations and design matrix
mediantime <- median(day1)
latgrid <- seq( min(lat1), max(lat1), length.out = 60 )
longrid <- seq( 0, 360, length.out = 121)[1:120] # so no locations repeated
locs_pred <- as.matrix( expand.grid(longrid,latgrid) )
n_pred <- nrow(locs_pred)
locstime_pred <- cbind( locs_pred, rep(mediantime, n_pred) )
X_pred <- as.matrix( cbind( rep(1,n_pred), locs_pred[,2], locs_pred[,2]^2 ) )
# predictions
pred <- predictions(fit4$covparms, "arma_matern_spheretime_warp", temp1,
    lonlattime1, locstime_pred, X1, X_pred, fit4$betahat, m = 40)
# plot predictions
par(mar=c(4,4,1,1))
pred_array <- array( pred, c(length(longrid),length(latgrid)) )
fields::image.plot(longrid,latgrid,pred_array)
# conditional simulations
sim <- cond_sim(fit4$covparms, "arma_matern_spheretime_warp", temp1,
    lonlattime1, locstime_pred, X1, X_pred, fit4$betahat, m = 40)
# plot conditional simulations
sim_array <- array( sim, c(length(longrid),length(latgrid)) )
fields::image.plot(longrid,latgrid,sim_array)


