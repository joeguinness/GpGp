
context("Check whether functions work in 1 dimension")


test_that("we can simulate, fit, predict, condsim in 1d", {

    # give locs and X as vectors (functions should convert to matrix)
    n <- 200
    locs <- seq(0,1,length.out=n)
    X <- rep(1,n)
    beta <- 1
    covfun_name <- "matern_isotropic"
    covparms <- c(1,0.2,1.5,0.1)

    # simulate data (keeping X as a vector)
    y <- X*beta + fast_Gp_sim(covparms=covparms, covfun_name=covfun_name, locs=locs )
    expect_equal( length(y), n )
    
    # fit a model
    fit <- fit_model( y, locs, X, covfun_name, silent = TRUE )

    # predictions and cond_sim
    X_pred <- X
    locs_pred <- locs

    pred <- predictions( fit, locs_pred, X_pred )
    csim <- cond_sim( fit, locs_pred, X_pred, nsims = 4 )
    expect_equal( length(pred), n )
    expect_equal( dim(csim), c(n,4) )

})
