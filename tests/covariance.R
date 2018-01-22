

# test to make sure the covariance function is doing the right thing
covparms <- c(2,1/4,1/2,0.1)
n1 <- 40
n2 <- 40
n <- n1*n2
locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )

cov1 <- covparms[1]*( exp( -fields::rdist(locs)/covparms[2] ) + covparms[4]*diag(n) ) 
cov2 <- maternIsotropic(locs,covparms) 

sum( (cov1-cov2)^2 )



