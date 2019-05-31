context("Fisher Scoring")

test_that("Fisher scoring converges quickly", {

    set.seed(12345)
    
    n1 <- 20
    n2 <- 20
    n <- n1*n2
    locs <- as.matrix( expand.grid( 1:n1, 1:n2 ) )
    ord <- order_maxmin(locs)
    locsord <- locs[ord,]
    m <- 20
    NNarray <- find_ordered_nn(locsord,m=m)
    NNlist <- group_obs(NNarray)
    Xord <- as.matrix(rep(1,n))
    
    covparms <- c(2,n1/4,0.01)
    yord <- fast_Gp_sim(covparms,"exponential_isotropic",locsord)
    fit1 <- fit_model( yord, locsord, X = Xord, m_seq = c(10,20),
        covfun_name = "exponential_isotropic", group = FALSE, 
        reorder = FALSE, NNarray = NNarray, silent = TRUE, convtol = 1e-6 )
    expect_lt( sum(abs(fit1$grad)), 1e-3 )
    
    fit2 <- fit_model( yord, locsord, X = Xord, m_seq = c(10,20),
        covfun_name = "exponential_isotropic", group = TRUE, 
        reorder = FALSE, NNarray = NNarray, silent = TRUE, convtol = 1e-6 )
    expect_lt( sum(abs(fit2$grad)), 1e-3 )
    
    LL11 <- vecchia_profbeta_loglik( fit1$covparms, "exponential_isotropic",
        yord, Xord, locsord, NNarray )
    expect_equal( fit1$betahat, LL11$betahat )
    expect_equal( fit1$loglik, LL11$loglik )
    
    LL12 <- vecchia_meanzero_loglik( fit1$covparms, "exponential_isotropic",
        c(yord - Xord %*% fit1$betahat), locsord, NNarray )
    expect_equal( fit1$loglik, LL12$loglik )
    
    LL21 <- vecchia_grouped_profbeta_loglik( fit2$covparms, "exponential_isotropic",
        yord, Xord, locsord, NNlist )
    expect_equal( fit2$betahat, LL21$betahat )
    expect_equal( fit2$loglik, LL21$loglik )
    
    LL22 <- vecchia_grouped_meanzero_loglik( fit2$covparms, "exponential_isotropic",
        c(yord - Xord %*% fit2$betahat), locsord, NNlist )
    expect_equal( fit2$loglik, LL22$loglik )



})    
    