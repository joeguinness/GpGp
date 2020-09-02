

#devtools::install_github("joeguinness/GpGp")
#library("GpGp")


#######################################################
#      Part 1: GpGp core functions                    #
#######################################################


# grid size for data locations
gsize <- 100
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

# covariance function and parameters
covparms <- c(variance = 4, range = 0.1, nugget = 0.1)

# simulate some data
y <- GpGp::fast_Gp_sim(covparms, "exponential_isotropic",locs,20)

# generate an ordering
ord <- GpGp::order_maxmin(locs)

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]
X <- matrix( rep(1,n), n, 1 )
Xord <- X[ord,,drop=FALSE]

# find the ordered m nearest neighbors
m <- 30
NNarray <- GpGp::find_ordered_nn(locsord,m)

# calculate likelihoods
t1 <- proc.time()
ll1 <- GpGp::vecchia_profbeta_loglik(
           covparms,"exponential_isotropic",yord,Xord,locsord,NNarray )
t2 <- proc.time()
ll2 <- GpGp::vecchia_profbeta_loglik_grad_info(
           covparms,"exponential_isotropic",yord,Xord,locsord,NNarray )
t3 <- proc.time()

print(t2-t1)
print(t3-t2)

# fit a the model
t1 <- proc.time()
gpfit <- GpGp::fit_model( y = y, locs = locs, covfun = "exponential_isotropic")
print(proc.time() - t1)




#######################################################
#      Part 2: Heaton et al. Temperature Data         #
#######################################################

library("fields")
library("GpGp")

# load the data and peek at it
load(url("http://guinness.cals.cornell.edu/modis_temps.RData"))
modis_temps[1:30,]
nrow(modis_temps)

# make a default image.plot
temp_dims <- c(500,300)
temp_array <- array( modis_temps$subtemp, temp_dims )
png("temps.png",width=600,height=600)
image.plot(temp_array)
dev.off()

# detect missing values
not_missing <- !is.na(modis_temps$subtemp)

# define the response
y <- modis_temps$subtemp[not_missing]

# locations (Longitude and Latitude)
locs <- cbind( modis_temps$lon[not_missing], modis_temps$lat[not_missing]  ) 

# linear covariates (this just uses lon and lat as covariates)
X <- cbind( rep(1,sum(not_missing)), locs )

# fit an exponential covariance GP model using GpGp::fit_model
t1 <- proc.time()
gpfit <- fit_model( y, locs, X, covfun_name = "exponential_sphere" ) 
t2 <- proc.time()
print(t2-t1)

# print out a summary
summary(gpfit)

# interpretation of range parameter
earth_radius <- 6356    # in kilometers
gpfit$covparms[2]*earth_radius

# set up prediction locations and design matrix
pred_inds <- is.na( modis_temps$subtemp )
locs_pred <- as.matrix(modis_temps[pred_inds, c("lon","lat")])
X_pred <- cbind(rep(1,sum(pred_inds)),locs_pred)

# do prediction of the missing data
t1 <- proc.time()
pred <- predictions( fit = gpfit, locs_pred = locs_pred, X_pred = X_pred, m = 30)
print(proc.time()-t1)

# make a map of the predictions
obs_and_pred <- modis_temps$subtemp
obs_and_pred[pred_inds] <- pred

png("temps_preds.png",width=1200,height=600)
par(mfrow=c(1,2))
image.plot( array(modis_temps$subtemp, temp_dims), axes = FALSE )
mtext("Data")
image.plot( array(obs_and_pred, temp_dims), axes = FALSE )
mtext("Conditional Expectation (Prediction)")
dev.off()

# do conditional simulations
# SIMULATE missing data, but consistent with the observed data
ncondsim <- 30
t1 <- proc.time()
sims <- cond_sim( fit = gpfit, locs_pred = locs_pred, X_pred = X_pred, 
    nsims = ncondsim, m = 30)
print(proc.time()-t1)

# prediction standard deviations
sum_squares <- rep(0, 500*300)
for(j in 1:30){
    obs_and_sim <- modis_temps$subtemp
    obs_and_sim[pred_inds] <- sims[,j]
    sum_squares <- sum_squares + (obs_and_pred - obs_and_sim)^2
}
pred_rmse <- sqrt( 1/30*sum_squares )

# plot of data, predictions, 1 conditional simulation, prediction standard devs
png("temps_preds_sims.png",width=800,height=800)
par(mfrow=c(2,2), mar=c(1,1,3,3))
image.plot( array(modis_temps$subtemp, temp_dims), axes = FALSE )
mtext("Data")
image.plot( array(obs_and_pred, temp_dims), axes = FALSE )
mtext("Conditional Expectation (Prediction)")
image.plot( array(obs_and_sim, temp_dims), axes = FALSE )
mtext("One Conditional Simulation")
image.plot( array(pred_rmse, temp_dims), axes = FALSE )
mtext("Prediction Standard Deviation")
dev.off()










#######################################################
#      Part 3: Spatial-Temporal Windspeed Data        #
#######################################################


library("maps")
library("GpGp")

## Reading in Data
data("jason3")
head(jason3)
lat       <- jason3$lat
lon       <- jason3$lon
windspeed <- jason3$windspeed
time      <- jason3$time/3600
n         <- length(windspeed)

## Visualizing the Data

# plot of data from first 6 hours
pdf("jason1.pdf",width=8,height=6)
par(mar=c(4,4,1,1))
inds <- time < 6
fields::quilt.plot(lon[inds],lat[inds],windspeed[inds],
    nx=400,ny=200,xlab="Lon",ylab="Lat",legend.lab = "windspeed (m/s)")
map("world2",add=TRUE)
dev.off()

# plot of data from first day
pdf("jason2.pdf",width=8,height=6)
par(mar=c(4,4,1,1))
inds <- time < 24
fields::quilt.plot(lon[inds],lat[inds],windspeed[inds],
    nx=400,ny=200,xlab="Lon",ylab="Lat",legend.lab = "windspeed (m/s)")
map("world2",add=TRUE)
dev.off()

# plot of all data
pdf("jason3.pdf",width=8,height=6)
par(mar=c(4,4,1,1))
fields::quilt.plot(lon,lat,windspeed,
    nx=400,ny=200,xlab="Lon",ylab="Lat",legend.lab = "windspeed (m/s)")
map("world2",add=TRUE)
dev.off()


## Preparing Variables for Fitting

# longitude/latitude locations, space-time locations, and design matrix
locs <- cbind(lon,lat)
locstime <- cbind(locs,time)
X <- as.matrix( rep(1,length(windspeed)) )

## Model Fitting

t1 <- proc.time()
fit1 <- fit_model(
    windspeed,
    locstime,
    X,
    "exponential_spheretime",
    st_scale = c(0.1,24)
)
print(proc.time() - t1)
summary(fit1)

# try updating the st_scale based on fit1
fit2 <- fit_model(
    windspeed,
    locstime,
    X,
    "exponential_spheretime",
    st_scale = fit1$covparms[2:3],
    start_parms = fit1$covparms
)

# look at the differences
round( fit1$covparms, 4 )
round( fit2$covparms, 4 )

