
# analysis of MODIS temperature data from Heaton et al (2018)

# Download data from here:
# http://guinness.cals.cornell.edu/paper_spatial_data_comparison.html

# install.packages("GpGp")

# attach fields and GpGp packages
library("fields")
library("GpGp")

# load the data and peek at it
load("modis_temps.RData")
modis_temps[1:30,]
nrow(modis_temps)

# make a default image.plot
temp_dims <- c(500,300)
temp_array <- array( modis_temps$temp, temp_dims )
image.plot(temp_array)

# subsampled data (to make things run faster
image.plot(array(modis_temps$subtemp,temp_dims) )

# detect missing values
not_missing <- !is.na(modis_temps$subtemp)

# define the response
y <- modis_temps$subtemp[not_missing]

# locations (Longitude and Latitude)
locs <- cbind( modis_temps$lon[not_missing], modis_temps$lat[not_missing]  ) 

# linear covariates
X <- cbind( rep(1,sum(not_missing)), locs )

# fit an exponential covariance GP model using GpGp::fit_model
?fit_model
t1 <- proc.time()
gpfit <- fit_model( y = y, locs = locs, X = X, "exponential_sphere" )
print(proc.time() - t1)
# exponential covariance = theta1 * ( exp( -dist/theta2 ) + theta3 * ( dist==0 ) )

# print out a summary
summary(gpfit)

earth_radius <- 6356    # in kilometers
gpfit$covparms[2]*earth_radius

# set up prediction locations and design matrix
pred_inds <- is.na( modis_temps$subtemp )
locs_pred <- as.matrix(modis_temps[pred_inds, c("lon","lat")])
X_pred <- cbind(rep(1,sum(pred_inds)),locs_pred)

# do prediction of the missing data
?predictions
t1 <- proc.time()
pred <- predictions( fit = gpfit, locs_pred = locs_pred, X_pred = X_pred, m = 30)
print(proc.time()-t1)

# make a map of the predictions
obs_and_pred <- modis_temps$subtemp
obs_and_pred[pred_inds] <- pred

par(mfrow=c(1,2))
image.plot( array(modis_temps$subtemp, temp_dims), axes = FALSE )
mtext("Data")
image.plot( array(obs_and_pred, temp_dims), axes = FALSE )
mtext("Conditional Expectation (Prediction)")

# do conditional simulations
# SIMULATE missing data, but consistent with the observed data
ncondsim <- 30
t1 <- proc.time()
sims <- cond_sim( fit = gpfit, locs_pred = locs_pred, X_pred = X_pred, 
    nsims = ncondsim, m = 30)
print(proc.time()-t1)

# plot 5 of the 30 conditional simulations
zlims <- range( c(modis_temps$subtemp, sims[,1:5] ), na.rm = TRUE )
for(j in 1:5){
    obs_and_sim <- modis_temps$subtemp
    obs_and_sim[pred_inds] <- sims[,j]
    par(mfrow=c(1,2))
    image.plot( array(modis_temps$subtemp, temp_dims), axes=FALSE, zlim=zlims )
    mtext("Data")
    image.plot( array(obs_and_sim, temp_dims), axes = FALSE, zlim = zlims )
    mtext(paste("Conditional Simulation",j))
}

# prediction standard deviations
sum_squares <- rep(0, 500*300)
for(j in 1:30){
    obs_and_sim <- modis_temps$subtemp
    obs_and_sim[pred_inds] <- sims[,j]
    sum_squares <- sum_squares + (obs_and_pred - obs_and_sim)^2
}
pred_rmse <- sqrt( 1/30*sum_squares )

# plot of data, predictions, 1 conditional simulation, prediction standard devs
par(mfrow=c(2,2), mar=c(1,1,3,3))
image.plot( array(modis_temps$subtemp, temp_dims), axes = FALSE )
mtext("Data")
image.plot( array(obs_and_pred, temp_dims), axes = FALSE )
mtext("Conditional Expectation (Prediction)")
image.plot( array(obs_and_sim, temp_dims), axes = FALSE )
mtext("One Conditional Simulation")
image.plot( array(pred_rmse, temp_dims), axes = FALSE )
mtext("Prediction Standard Deviation")


# guinness@cornell.edu

