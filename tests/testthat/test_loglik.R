
context("Likelihood Functions")

test_that("likelihood approximations are exact when m = n-1", {

    n1 <- 10
    n2 <- 10
    n <- n1*n2
    locs <- as.matrix( expand.grid( 1:n1, 1:n2 ) )
    ord <- order_maxmin(locs)
    locsord <- locs[ord,]
    m <- n-1
    NNarray <- find_ordered_nn(locsord,m=m)
    NNlist <- group_obs(NNarray)

    covparms <- c(2,40,0.8,0.01)
    y <- fast_Gp_sim(covparms,"matern_isotropic",locsord)
    expect_equal(length(y),n)

    X <- as.matrix(cbind(rep(1,n),runif(n)))
    covmat <- matern_isotropic(covparms,locsord)
    cholmat <- t(chol(covmat))
    Z <- solve(cholmat,X)
    betahat <- c( solve( crossprod(Z), crossprod(Z,solve(cholmat,y)) ) )
    logdet_exact <- 2*sum(log(diag(cholmat)))
    z <- forwardsolve(cholmat,y - X %*% betahat)
    quadform_exact <- c(crossprod(z))
    ll_exact <- -n/2*log(2*pi) - 1/2*logdet_exact - 1/2*quadform_exact
    
    # ungrouped
    llobj <- vecchia_profbeta_loglik_grad_info(covparms,"matern_isotropic",
        y,X,locsord,NNarray)
    expect_equal( llobj$loglik, ll_exact )
    expect_equal( llobj$betahat, betahat )
    
    # grouped
    llobj <- vecchia_grouped_profbeta_loglik_grad_info(covparms,"matern_isotropic",
        y,X,locsord,NNlist)
    expect_equal( llobj$loglik, ll_exact )
    expect_equal( llobj$betahat, betahat )

})



test_that("gradient of likelihood matches finite differencing", {

    n1 <- 10
    n2 <- 10
    n <- n1*n2
    locs <- as.matrix( expand.grid( 1:n1, 1:n2 ) )
    ord <- order_maxmin(locs)
    locsord <- locs[ord,]
    m <- 20
    NNarray <- find_ordered_nn(locsord,m=m)
    NNlist <- group_obs(NNarray)

    covparms <- c(2,40,0.8,0.01)
    y <- fast_Gp_sim(covparms,"matern_isotropic",locsord)

    X <- as.matrix(cbind(rep(1,n),runif(n)))
    covmat <- matern_isotropic(covparms,locsord)
    cholmat <- t(chol(covmat))
    Z <- solve(cholmat,X)
    betahat <- c( solve( crossprod(Z), crossprod(Z,solve(cholmat,y)) ) )
    logdet_exact <- 2*sum(log(diag(cholmat)))
    z <- forwardsolve(cholmat,y - X %*% betahat)
    quadform_exact <- c(crossprod(z))
    ll_exact <- -n/2*log(2*pi) - 1/2*logdet_exact - 1/2*quadform_exact
    
    # ungrouped
    llobj <- vecchia_profbeta_loglik_grad_info(covparms,"matern_isotropic",
        y,X,locsord,NNarray)
    dll <- rep(NA,length(covparms))
    eps <- 1e-8
    for(j in 1:length(covparms)){
        dcovparms <- covparms
        dcovparms[j] <- covparms[j]+eps
        dllobj <- vecchia_profbeta_loglik_grad_info(dcovparms,"matern_isotropic",
        y,X,locsord,NNarray)
        dll[j] <- (dllobj$loglik - llobj$loglik)/eps
    }
    denom <- abs(llobj$grad)
    denom[ denom < 1 ] <- 1
    expect_equal( llobj$grad/denom, dll/denom, tolerance = 1e-4 )
    chol( llobj$info )

    # grouped
    llobj <- vecchia_grouped_profbeta_loglik_grad_info(covparms,"matern_isotropic",
        y,X,locsord,NNlist)
    dll <- rep(NA,length(covparms))
    eps <- 1e-8
    for(j in 1:length(covparms)){
        dcovparms <- covparms
        dcovparms[j] <- covparms[j]+eps
        dllobj <- vecchia_grouped_profbeta_loglik_grad_info(
            dcovparms,"matern_isotropic",y,X,locsord,NNlist)
        dll[j] <- (dllobj$loglik - llobj$loglik)/eps
    }
    denom <- abs(llobj$grad)
    denom[ denom < 1 ] <- 1
    expect_equal( llobj$grad/denom, dll/denom, tolerance = 1e-4 )
    chol( llobj$info )
    expect_equal( llobj$info , t(llobj$info) )

})