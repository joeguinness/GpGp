
# vignette demonstrating importance of optimization
# grid size for data locations
devtools::load_all()
gsize <- 20
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))
# generate basis functions
basgridsize <- 3
basvec <- 0:(basgridsize-1)/(basgridsize-1)
basknots <- as.matrix(expand.grid( basvec, basvec ))
basrange <- 2.0/(basgridsize-1)
Zbas <- matrix(NA, n, nrow(basknots) )
for(j in 1:ncol(Zbas)){
    Zbas[,j] <- exp(-(fields::rdist(locs,basknots[j,,drop=FALSE])/basrange)^2)
}
Z <- cbind(locs,Zbas)
#Z <- locs

# covariance function and parameters
# covfun <- maternIsotropic
#cov1 <- "exponential_isotropic"
#covparms <- c(variance = 4, range = 0.1, nugget = 0.1)
cov1 <- "matern_isotropic"
covparms1 <- c(variance = 4, range = 0.3, smoothness = 0.5, nugget = 0.00)
#cov1 <- "matern_isotropic"
#covparms <- c(variance = 4, range = 0.1, smoothness = 0.75, nugget = 0.01)
#cov1 <- "matern_anisotropic2D"
#cov2 <- paste0("arma_",cov1)
#cov2 <- "arma_matern_nonstat_var"
cov2 <- "arma_matern_isotropic"


# simulate some data
y <- fast_Gp_sim(covparms1, cov1 ,locs,40)
X <- matrix( rep(1,n), n, 1 )
NNarray <- find_ordered_nn(locs,m=30)
NNlist <- group_obs(NNarray)
m <- dim(NNarray)[2]

set.seed( as.numeric(Sys.time()) )
t1 <- proc.time()
m1 <- fit_model_arma(y,Z,X,covfun_name = cov2, 
    NNarray = NNarray, reorder=FALSE, group = FALSE)
t2 <- proc.time()
m2 <- fit_model_fisher(y,Z,X,covfun_name = cov2, 
    NNarray = NNarray, reorder=FALSE, group = FALSE, silent = TRUE)
t3 <- proc.time()
m3 <- fit_model_fisher(y,Z,X,covfun_name = cov2, 
    NNarray = NNarray, reorder=FALSE, group = TRUE, silent = TRUE)
t4 <- proc.time()

#fields::image.plot(array(y,nvec))

print(t2-t1)
print(t3-t2)
print(t4-t3)

lls <- c(m1$loglik,m2$loglik,m3$loglik)

cat(sprintf("%11.9f %11.9f %11.9f \n",max(lls)-lls[1],max(lls)-lls[2],max(lls)-lls[3]))
round(m1$covparms,4)
round(m2$covparms,4)
round(m3$covparms,4)


if(FALSE){

sv <- arma_vecchia_Linv( c(2, 0.5, 2.5, 0), cov1, locs, NNarray  )
    
# testing timing of pieces for grouped version
t1 <- proc.time()
sv00 <- vecchia_loglik(covparms1,cov1,y,locs,NNarray)
t2 <- proc.time()
sv0 <- vecchia_loglik_grouped(covparms1,cov1,y,locs,NNlist)
t3 <- proc.time()
sv2 <- arma_vecchia_grad_hess(covparms1,cov2,y,X,locs,NNarray)
t4 <- proc.time()
sv3 <- arma_vecchia_grad_hess_grouped(covparms1,cov2,y,X,locs,NNlist)
t5 <- proc.time()
sv4 <- arma_vecchia_loglik(covparms1,cov2,y,locs,NNarray)
t6 <- proc.time()

print(t2-t1)
print(t3-t2)
print(t4-t3)
print(t5-t4)
print(t6-t5)

print(sv0 - sv3$loglik)
print(sv00 - sv2$loglik)

    
# testing timing of pieces
t1 <- proc.time()
sv0 <- arma_vecchia_loglik(covparms,cov2,y,locs,NNarray)
t2 <- proc.time()
sv1 <- arma_vecchia_grad_hess(covparms,cov2,y,X,locs,NNarray)
t3 <- proc.time()
sv2 <- arma_vecchia_grad_hess2(covparms,cov2,y,X,locs,NNarray)
t4 <- proc.time()

print(t2-t1)
print(t3-t2)
print(t4-t3)
print(round(
c(sv1$loglik-sv2$loglik,sv1$grad-sv2$grad,c(sv1$info-sv2$info),sv1$betahat-sv2$betahat)
,4))





covparms <- m2$covparms
#vv <- seq(0.1,3,length.out=21)
vv <- exp(seq(-1.0,1.0,length.out=40))
llarr <- array(NA, rep(length(vv),2))
for(j1 in 1:dim(llarr)[1]){ 
    cat(paste(j1," "))
for(j2 in 1:dim(llarr)[2]){
    dcovparms <- covparms
    dcovparms[2] <- dcovparms[2]*vv[j1]
    dcovparms[3] <- dcovparms[3]*vv[j2]
    llarr[j1,j2] <- proflik_mean(dcovparms,cov1,y,X,locs,NNarray)
}}
cat("\n")
par(mar=c(3,3,1,1))
fields::image.plot( log(covparms[2]*vv), log(covparms[3]*vv), llarr,
    zlim = c(m1$loglik-600,m1$loglik))
points(log(m1$covparms[2]), log(m1$covparms[3]), pch = 15, col = "white" )
points(log(m2$covparms[2]), log(m2$covparms[3]), pch = 16, col = "white" )
points(log(m3$covparms[1]), log(m3$covparms[2]), pch = 17, col = "white" )

par(mar=c(3,3,1,1))
contour( log(covparms[2]*vv), log(covparms[3]*vv), llarr,
    levels = c(-7500,-7200,-6900, -6800, -6700, -6600, -6560) )
points(log(m1$covparms[2]), log(m1$covparms[3]), pch = 15, col = "white" )
points(log(m2$covparms[2]), log(m2$covparms[3]), pch = 16, col = "white" )
points(log(m3$covparms[1]), log(m3$covparms[2]), pch = 17, col = "white" )


m1$covparms
m2$covparms

vecchia_loglik(m1$covparms, cov1, y-c(m1$beta), locs, NNarray ) -
proflik_mean_variance(m1$covparms[2:length(covparms)], cov1,y,X,locs,NNarray )

arma_vecchia_loglik(m2$covparms, cov2,y-c(m2$betahat),locs,NNarray) -
arma_vecchia_loglik(m1$covparms, cov2,y-c(m1$beta),locs,NNarray)


subparms <- m3$covparms
ddobj <- arma_vecchia_profile_grad_hess(subparms,cov2,y,X,locs,NNarray)
covparms <- c(ddobj$sigmasq,subparms)
sv1 <- arma_vecchia_grad_hess(covparms,cov2,y,X,locs,NNarray)

dd <- ddobj$grad
dde <- rep(NA, length(dd))
eps <- 1e-8
for(j in 1:3){
    dsubparms <- subparms
    dsubparms[j] <- subparms[j]+eps
    ll <- arma_vecchia_profile_grad_hess(dsubparms,cov2,y,X,locs,NNarray)$loglik
    dde[j] <- (ll - ddobj$loglik)/eps
}




eps <- 1e-8
covparms <- m2$covparms
lloc <- locs[ NNarray[50,] , ]
dd <- d_arma_matern_isotropic(covparms,lloc)
dde <- array(NA, dim(dd))
for(j in 1:4){
    dcovparms <- covparms
    dcovparms[j] <- covparms[j]+eps
    cove <- arma_matern_isotropic(dcovparms,lloc)
    cov <- arma_matern_isotropic(covparms,lloc)
    dde[,,j] <- (cove-cov)/eps
}
max(abs(dde-dd))
for(j in 1:4){
    print(max(abs(dde[,,j]-dd[,,j])))
}


# compare matern anisotropic2D with matern isotropic
devtools::load_all()
cov1 <- "arma_matern_isotropic"
cov2 <- "arma_matern_anisotropic2D"
covparms1 <- c(2,0.2,0.5,0.1)
covparms2 <- c(covparms1[1],1/covparms1[2],0.0,1/covparms1[2],covparms1[3:4])
n <- 10
locs <- matrix( runif(2*n), n, 2 )
covmat1 <- arma_matern_isotropic(covparms1,locs)
covmat2 <- arma_matern_anisotropic2D(covparms2,locs)

covmat1-covmat2

# check derivatives of matern anisotropic2D
n <- 10
locs <- matrix( runif(2*n), n, 2 )
covparms <- c(2,5,0.0,5,0.8,0.1)
dd <- d_arma_matern_anisotropic2D(covparms,locs)
eps <- 1e-8
dde <- array(NA, dim(dd))
cov <- arma_matern_anisotropic2D(covparms,locs)
for(j in 1:6){
    dcovparms <- covparms
    dcovparms[j] <- covparms[j]+eps
    print(dcovparms)
    cove <- arma_matern_anisotropic2D(dcovparms,locs)
    dde[,,j] <- (cove-cov)/eps
}
max(abs(dde-dd))
for(j in 1:6){
    print(max(abs(dde[,,j]-dd[,,j])))
}


# test Linv functions
gsize <- 70
nvec <- c(gsize,gsize)
n <- prod(nvec)
# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))
NNarray <- find_ordered_nn(locs,m=30)

cov1 <- "matern_isotropic"
covparms <- c(2,0.2,0.6,0.1)
cov2 <- paste0("arma_",cov1)
Linv1 <- vecchia_Linv(covparms,cov1,locs,NNarray)
Linv2 <- arma_vecchia_Linv(covparms,cov2,locs,NNarray)


# test nonstat covariance function
covparms <- c(2,0.2,0.75,0.1, rep(.1,ncol(Z)-2))
covmat2 <- arma_matern_nonstat_var(covparms,Z)
covmat1 <- arma_matern_isotropic(covparms[1:4], Z[,1:2])


# check derivatives of matern nonstat var
n <- 10
locs <- matrix( runif(2*n), n, 2 )
# generate basis functions
basgridsize <- 3
basvec <- 0:(basgridsize-1)/(basgridsize-1)
basknots <- as.matrix(expand.grid( basvec, basvec ))
basrange <- 2.0/(basgridsize-1)
Zbas <- matrix(NA, n, nrow(basknots) )
for(j in 1:ncol(Zbas)){
    Zbas[,j] <- exp(-(fields::rdist(locs,basknots[j,,drop=FALSE])/basrange)^2)
}
Z <- cbind(locs,Zbas)

covparms <- c(2,0.2, 0.8, 0.1, rep(0,ncol(Z)-1))
dd <- d_arma_matern_nonstat_var(covparms,Z)
eps <- 1e-8
dde <- array(NA, dim(dd))
cov <- arma_matern_nonstat_var(covparms,locs)
for(j in 1:length(covparms)){
    dcovparms <- covparms
    dcovparms[j] <- covparms[j]+eps
    cove <- arma_matern_nonstat_var(dcovparms,Z)
    dde[,,j] <- (cove-cov)/eps
}
max(abs(dde-dd))
for(j in 1:6){
    print(max(abs(dde[,,j]-dd[,,j])))
}






}