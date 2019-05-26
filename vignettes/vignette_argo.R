
# analysis of argo data
devtools::load_all()

read_and_save_data <- FALSE
if(read_and_save_data){
argo <- R.matlab::readMat("~/Dropbox/research/fisher_vecchia/data/Argo_data_aggr_02_extended.mat")

inds <- argo$profJulDayAggr > 736200

lon <- argo$profLongAggr[inds]
lat <- argo$profLatAggr[inds]
day <- argo$profJulDayAggr[inds]
n <- length(lon)

temp_prof <- argo$profTempAggr[inds]
pres_prof <- argo$profPresAggr[inds]

# write code to save the damn data!
argo2016Jan2Mar <- list( temp_prof = temp_prof, pres_prof = pres_prof, 
    day = day, lon = lon, lat = lat)
save(argo2016Jan2Mar, 
    file = "~/Dropbox/research/fisher_vecchia/data/argo2016Jan2Mar.RData")
}

load("~/Dropbox/research/fisher_vecchia/data/argo2016Jan2Mar.RData")
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

NNarray <- find_ordered_nn(xyz,m=30)

timing <- rep(NA,8)

# spatial isotropic
t1 <- proc.time()
fit1 <- fit_model_fisher5(y = temp, X = X, locs = lonlat, 
    covfun_name = "arma_matern_sphere", group = FALSE, 
    NNarray = NNarray, reorder = FALSE)
timing[1] <- (proc.time() - t1)[3]


# spacetime, isotropic in each
t1 <- proc.time()
fit2 <- fit_model_fisher5(y = temp, X = X, locs = lonlattime, 
    covfun_name = "arma_matern_spheretime", group = FALSE, 
    NNarray = NNarray, reorder = FALSE)
timing[2] <- (proc.time() - t1)[3]


# spatial warping
t1 <- proc.time()
fit3 <- fit_model_fisher5(y = temp, X = X, locs = lonlat, 
    covfun_name = "arma_matern_sphere_warp", group = FALSE, 
    NNarray = NNarray, reorder = FALSE)
timing[3] <- (proc.time() - t1)[3]


# spacetime, warping in space
t1 <- proc.time()
fit4 <- fit_model_fisher5(y = temp, X = X, locs = lonlattime, 
    covfun_name = "arma_matern_spheretime_warp", group = FALSE, 
    NNarray = NNarray, reorder = FALSE)
timing[4] <- (proc.time() - t1)[3]






# spatial isotropic
t1 <- proc.time()
fit5 <- fit_model_arma(y = temp, X = X, locs = lonlat, 
    covfun_name = "arma_matern_sphere", group = FALSE, 
    NNarray = NNarray, reorder = FALSE)
timing[5] <- (proc.time() - t1)[3]


# spacetime, isotropic in each
t1 <- proc.time()
fit6 <- fit_model_arma(y = temp, X = X, locs = lonlattime, 
    covfun_name = "arma_matern_spheretime", group = FALSE, 
    NNarray = NNarray, reorder = FALSE)
timing[6] <- (proc.time() - t1)[3]


# spatial warping
t1 <- proc.time()
fit7 <- fit_model_arma(y = temp, X = X, locs = lonlat, 
    covfun_name = "arma_matern_sphere_warp", group = FALSE, 
    NNarray = NNarray, reorder = FALSE)
timing[7] <- (proc.time() - t1)[3]


# spacetime, warping in space
t1 <- proc.time()
fit8 <- fit_model_arma(y = temp, X = X, locs = lonlattime, 
    covfun_name = "arma_matern_spheretime_warp", group = FALSE, 
    NNarray = NNarray, reorder = FALSE)
timing[8] <- (proc.time() - t1)[3]
#     user    system   elapsed 
# 14845.617   645.256 14648.249


save(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,timing,
    file="argo_fitting_results.RData")


