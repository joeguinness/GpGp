

data("wind", package = "gstat" )
library("lubridate")
devtools::load_all()

wind_date <- ymd( paste0("19",wind$year,"-",wind$month,"-",wind$day) )
wind_day <- as.numeric(wind_date) - min(as.numeric(wind_date))

num_days <- nrow(wind)
num_stations <- ncol(wind) - 3
wind_time <- rep( wind_day, num_stations)



# get the longitudes and latitudes
lon0 <- strsplit( as.character(wind.loc$Longitude), "[d'\"]+")
lon0 <- lapply(lon0,as.numeric)
lat0 <- strsplit( as.character(wind.loc$Latitude), "[d'\"]+")
lat0 <- lapply(lat0,as.numeric)
station_lon <- rep(NA,num_stations)
station_lat <- rep(NA,num_stations)
for(j in 1:(num_stations)){
    station_lon[j] <- lon0[[j]][1] + lon0[[j]][2]/60
    station_lat[j] <- lat0[[j]][1] + lat0[[j]][2]/60
}
station_lon[12] <- station_lon[12] + lon0[[12]][3]/3600
station_lat[12] <- station_lat[12] + lat0[[12]][3]/3600

# plot the data
# time series plot
for(j in 1:num_stations){
    plot( wind_day/365.25, sqrt(wind[,j+3]), type="l",
        ylim=c(0,7), main = colnames(wind)[j+3] )
}

# temporal correlations
for(j in 1:num_stations){
    print(cor( sqrt(wind[1:(num_days-1),j+3]), sqrt(wind[2:num_days,j+3])))
}

# spatial plot
plot(-station_lon,station_lat,type="n",xlim=c(-12,-5),ylim=c(51,56))
for(j in 1:num_stations){
    mn <- mean(wind[,j+3])
    text(-station_lon[j], station_lat[j], round(mn,2))
    text(-station_lon[j], station_lat[j]+0.2, colnames(wind)[j+3])
}
library("maps")
library("mapdata")
map(  database = "worldHires", regions = "ireland", add = TRUE,
    col="gray")


# get the vector of locations
wind_locs <- cbind( rep(station_lon,each=num_days), rep(station_lat,each=num_days) )
wind_station_name <- rep( wind.loc$Code, each = num_days )
wind_locstime <- cbind(wind_locs, wind_day)

# get the response
wind_sqrt0_mat <- as.matrix(wind[,4:(num_stations+3)])
for(j in 1:num_stations){
    wind_sqrt0_mat[,j] <- wind_sqrt0_mat[,j] - mean(wind_sqrt0_mat[,j])
}
wind_sqrt0 <- c(wind_sqrt0_mat)
n <- length(wind_sqrt0)

# get covariates
X <- cbind(rep(1,n), sin(2*pi*wind_time/365.25), cos(2*pi*wind_time/365.25))

covfun_fit <- c(
    "arma_exponential_isotropic",
    "arma_matern_isotropic",
    "arma_matern_spacetime",
    "arma_exponential_anisotropic3D" 
)

fits <- vector("list", length = length(covfun_fit) )
start_list <- list(c(4,2,0.1),c(4,2,0.8,0.1),c(4,2,2,0.8,0.1),
    c(4,1/2,0,1/2,0,0,1/2,0.1))
#inds <- sample(n,10000)
inds <- 1:n
for(j in 1:length(covfun_fit)){
    start_parms <- start_list[[j]]
    fits[[j]] <- fit_model_fisher5(wind_sqrt0[inds], wind_locstime[inds,], 
        X[inds,], covfun_name = covfun_fit[j], start_parms = start_parms,
        group=FALSE)
}



# plot the covariances
xgrid <- seq(-5,5,length.out=20)
ygrid <- seq(-5,5,length.out=20)
tt <- seq(0,10)
xyt <- as.matrix(expand.grid(xgrid,ygrid,tt))
nx <- length(xgrid)
ny <- length(ygrid)
covmat <- arma_exponential_anisotropic3D(fits[[4]]$covparms,xyt)



implot <- fields::image.plot
par(mfrow=c(1,3))
for(j in 1:3){
    implot( matrix( covmat[(j-1)*nx*ny+(1:(nx*ny)),210], nx, ny ),zlim=c(0,25) )
}






