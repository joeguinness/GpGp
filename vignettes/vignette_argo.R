
# analysis of argo data
devtools::load_all()
data("argo2016Jan2Mar")
names(argo2016Jan2Mar)
attach(argo2016Jan2Mar)
n <- length(lon)

chosen_pres <- 100
par(mfrow=c(1,2))
plot(lon,lat,pch=1,cex=0.3)
hist(day,breaks=60)

temp <- rep(NA,n)
for(j in 1:n){
    ii <- which.min( abs(c(pres_prof[[j]][[1]]) - chosen_pres ) )
    temp[j] <- temp_prof[[j]][[1]][ii,1]
}


fields::quilt.plot(lon,lat,temp,nx=200,ny=200)

lonlat <- cbind(lon,lat)
lonlattime <- cbind(lonlat,day)

X <- cbind(rep(1,n),lat,lat^2)

xyz <- array(NA, c(nrow(lonlat), 3) )
lonrad <- lon*2*pi/360
latrad <- (lat+90)*2*pi/360
xyz[,1] <- sin(latrad)*cos(lonrad)
xyz[,2] <- sin(latrad)*sin(lonrad)
xyz[,3] <- cos(latrad)



# subset of data for faster cran checks
inds <- sample(n,2000)

#NNarray <- find_ordered_nn(xyz,m=30)
timing <- rep(NA,4)

# spatial isotropic
t1 <- proc.time()
fit1 <- fit_model(y = temp[inds], X = X[inds,], locs = lonlat[inds,], 
    covfun_name = "arma_matern_sphere", group = FALSE )
timing[1] <- (proc.time() - t1)[3]


# spacetime, isotropic in each
t1 <- proc.time()
fit2 <- fit_model(y = temp[inds], X = X[inds,], locs = lonlattime[inds,], 
    covfun_name = "arma_matern_spheretime", group = FALSE )
timing[2] <- (proc.time() - t1)[3]


# spatial warping
t1 <- proc.time()
fit3 <- fit_model(y = temp[inds], X = X[inds,], locs = lonlat[inds,], 
    covfun_name = "arma_matern_sphere_warp", group = FALSE )
timing[3] <- (proc.time() - t1)[3]


# spacetime, warping in space
t1 <- proc.time()
fit4 <- fit_model(y = temp[inds], X = X[inds,], locs = lonlattime[inds,], 
    covfun_name = "arma_matern_spheretime_warp", group = FALSE )
timing[4] <- (proc.time() - t1)[3]





