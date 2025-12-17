context("Predictions")

test_that("Predictions work as expected when X_pred not specified", {

    set.seed(12345)
    
    n1 <- 20
    n2 <- 20
    n <- n1*n2
    locs <- as.matrix( expand.grid( 1:n1, 1:n2 ) )
    y <- fast_Gp_sim( c(1, n1/3, 0.01), "exponential_isotropic", locs, m = 30 )

    m1 <- fit_model( y, locs, covfun_name = "exponential_isotropic", silent = TRUE )

    locs_pred <- locs + matrix( 0.5, n, 2 )

    # do predictions with and without specifying X_pred, make sure they are similar
    p1 <- predictions( m1, locs_pred )
    p2 <- predictions( m1, locs_pred, X_pred = matrix(1, n, 1) )

    expect_lt( mean( (p1-p2)^2 )/mean(p2^2), 1e-4 )
    
    # fit a model with a covariate. Check that predictions gives warning when X_pred not specified
    x1 <- rnorm(n)
    X <- model.matrix( ~ x1 )
    m2 <- fit_model( y, locs, X, covfun_name = "exponential_isotropic", silent = TRUE )

    expect_error( predictions( m2, locs_pred ) )

})    
    
