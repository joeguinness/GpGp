
#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "onepass.h"
//[[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;


//' Vecchia's loglikelihood, gradient, and Fisher information
//'
//' This function returns Vecchia's (1988) approximation to the Gaussian
//' loglikelihood, profiling out the regression coefficients, and returning
//' the gradient and Fisher information. 
//' Vecchia's approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @inheritParams vecchia_meanzero_loglik
//' @param X Design matrix of covariates. Row \code{i} of \code{X} contains
//' the covariates for the observation at row \code{i} of \code{locs}.
//' @return A list containing 
//' \itemize{
//'     \item \code{loglik}: the loglikelihood
//'     \item \code{grad}: gradient with respect to covariance parameters
//'     \item \code{info}: Fisher information for covariance parameters
//'     \item \code{betahat}: profile likelihood estimate of regression coefs
//'     \item \code{betainfo}: information matrix for \code{betahat}.
//' }
//' The covariance matrix for \code{$betahat} is the inverse of \code{$betainfo}.
//' @examples
//' n1 <- 20
//' n2 <- 20
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' X <- cbind(rep(1,n),locs[,2])
//' covparms <- c(2, 0.2, 0.75, 0)
//' y <- X %*% c(1,2) + fast_Gp_sim(covparms, "matern_isotropic", locs, 50 )
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' loglik <- vecchia_profbeta_loglik_grad_info( covparms, "matern_isotropic", 
//'     y, X, locs, NNarray )
//' @export
// [[Rcpp::export]]
List vecchia_profbeta_loglik_grad_info( 
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    IntegerMatrix NNarray ){
    
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    synthesize(covparms, covfun_name, locs, NNarray, y, X,
        &ll, &betahat, &grad, &info, &betainfo, true, true 
    );
    
    List ret = List::create( Named("loglik") = ll, Named("betahat") = betahat,
        Named("grad") = grad, Named("info") = info, Named("betainfo") = betainfo );
    return ret;
        
}

//' Vecchia's approximation to the Gaussian loglikelihood, with profiled 
//' regression coefficients.
//'
//' This function returns Vecchia's (1988) approximation to the Gaussian
//' loglikelihood, profiling out the regression coefficients. 
//' The approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @inheritParams vecchia_meanzero_loglik
//' @param X Design matrix of covariates. Row \code{i} of \code{X} contains
//' the covariates for the observation at row \code{i} of \code{locs}.
//' @return a list containing
//' \itemize{
//'  \item \code{loglik}: the loglikelihood
//'  \item \code{betahat}: profile likelihood estimate of regression coefficients
//'  \item \code{betainfo}: information matrix for \code{betahat}.
//' }
//' The covariance
//' matrix for \code{$betahat} is the inverse of \code{$betainfo}.
//' @examples
//' n1 <- 20
//' n2 <- 20
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' X <- cbind(rep(1,n),locs[,2])
//' covparms <- c(2, 0.2, 0.75, 0)
//' y <- X %*% c(1,2) + fast_Gp_sim(covparms, "matern_isotropic", locs, 50 )
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' loglik <- vecchia_profbeta_loglik( covparms, "matern_isotropic", y, X, locs, NNarray )
//' @export
// [[Rcpp::export]]
List vecchia_profbeta_loglik( 
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    IntegerMatrix NNarray ){
    
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    synthesize(covparms, covfun_name, locs, NNarray, y, X,
        &ll, &betahat, &grad, &info, &betainfo, true, false 
    );
    
    List ret = List::create( Named("loglik") = ll, Named("betahat") = betahat,
        Named("betainfo") = betainfo );
    return ret;
        
}


//' Vecchia's approximation to the Gaussian loglikelihood, zero mean
//'
//' This function returns Vecchia's (1988) approximation to the Gaussian
//' loglikelihood. The approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @param covparms A vector of covariance parameters appropriate
//' for the specified covariance function
//' @param covfun_name See \code{\link{GpGp}} for information about covariance
//' functions.
//' @param y vector of response values
//' @param locs matrix of locations. Row \code{i} of \code{locs} specifies the location
//' of element \code{i} of \code{y}, and so the length of \code{y} should equal
//' the number of rows of \code{locs}.
//' @param NNarray A matrix of indices, usually the output from \code{\link{find_ordered_nn}}.
//' Row \code{i} contains the indices
//' of the observations that observation \code{i} conditions on. By convention,
//' the first element of row \code{i} is \code{i}.
//' @return a list containing
//' \itemize{
//'  \item \code{loglik}: the loglikelihood
//' }
//' @examples
//' n1 <- 20
//' n2 <- 20
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' covparms <- c(2, 0.2, 0.75, 0)
//' y <- fast_Gp_sim(covparms, "matern_isotropic", locs, 50 )
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' loglik <- vecchia_meanzero_loglik( covparms, "matern_isotropic", y, locs, NNarray )
//' @export
// [[Rcpp::export]]
List vecchia_meanzero_loglik( 
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    const NumericMatrix locs,
    IntegerMatrix NNarray ){
    
    NumericMatrix X(1,1);
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    //NumericVector betahat;
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );
    //NumericMatrix betainfo;

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    synthesize(covparms, covfun_name, locs, NNarray, y, X,
        &ll, &betahat, &grad, &info, &betainfo, false, false 
    );
    
    List ret = List::create( Named("loglik") = ll );
    return ret;
        
}



//' Grouped Vecchia loglikelihood, gradient, Fisher information
//'
//' This function returns a grouped version (Guinness, 2018) of Vecchia's (1988) 
//' approximation to the Gaussian
//' loglikelihood, the gradient, and Fisher information, 
//' and the profile likelihood estimate of the regression
//' coefficients. The approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @inheritParams vecchia_grouped_meanzero_loglik
//' @param X Design matrix of covariates. Row \code{i} of \code{X} contains
//' the covariates for the observation at row \code{i} of \code{locs}.
//' @return a list containing
//' \itemize{
//'     \item \code{loglik}: the loglikelihood
//'     \item \code{grad}: gradient with respect to covariance parameters
//'     \item \code{info}: Fisher information for covariance parameters
//'     \item \code{betahat}: profile likelihood estimate of regression coefs
//'     \item \code{betainfo}: information matrix for \code{betahat}.
//' }
//' The covariance
//' matrix for \code{$betahat} is the inverse of \code{$betainfo}.
//' @examples
//' n1 <- 20
//' n2 <- 20
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' X <- cbind(rep(1,n),locs[,2])
//' covparms <- c(2, 0.2, 0.75, 0)
//' y <- fast_Gp_sim(covparms, "matern_isotropic", locs, 50 )
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' NNlist <- group_obs(NNarray)
//' loglik <- vecchia_grouped_profbeta_loglik_grad_info( 
//'     covparms, "matern_isotropic", y, X, locs, NNlist )
//' @export
// [[Rcpp::export]]
List vecchia_grouped_profbeta_loglik_grad_info( 
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    List NNlist ){
    
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    // maybe the synthesize functions should take in an argument that
    // says which compute_pieces function to use
    synthesize_grouped(covparms, covfun_name, locs, NNlist, y, X,
        &ll, &betahat, &grad, &info, &betainfo, true, true 
    );
    
    
    List ret = List::create( Named("loglik") = ll, Named("betahat") = betahat,
        Named("grad") = grad, Named("info") = info, Named("betainfo") = betainfo );
    return ret;
        
}


//' Grouped Vecchia approximation, profiled regression coefficients
//'
//' This function returns a grouped version (Guinness, 2018) of Vecchia's (1988) 
//' approximation to the Gaussian
//' loglikelihood and the profile likelihood estimate of the regression
//' coefficients. The approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @inheritParams vecchia_grouped_meanzero_loglik
//' @param X Design matrix of covariates. Row \code{i} of \code{X} contains
//' the covariates for the observation at row \code{i} of \code{locs}.
//' @return a list containing
//' \itemize{
//'  \item \code{loglik}: the loglikelihood
//'  \item \code{betahat}: profile likelihood estimate of regression coefficients
//'  \item \code{betainfo}: information matrix for \code{betahat}.
//' }
//' The covariance
//' matrix for \code{$betahat} is the inverse of \code{$betainfo}.
//' @examples
//' n1 <- 20
//' n2 <- 20
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' X <- cbind(rep(1,n),locs[,2])
//' covparms <- c(2, 0.2, 0.75, 0)
//' y <- fast_Gp_sim(covparms, "matern_isotropic", locs, 50 )
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' NNlist <- group_obs(NNarray)
//' loglik <- vecchia_grouped_profbeta_loglik( 
//'     covparms, "matern_isotropic", y, X, locs, NNlist )
//' @export
// [[Rcpp::export]]
List vecchia_grouped_profbeta_loglik( 
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    List NNlist ){
    
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    // maybe the synthesize functions should take in an argument that
    // says which compute_pieces function to use
    synthesize_grouped(covparms, covfun_name, locs, NNlist, y, X,
        &ll, &betahat, &grad, &info, &betainfo, true, false 
    );
    
    
    List ret = List::create( Named("loglik") = ll, Named("betahat") = betahat,
        Named("betainfo") = betainfo );
    return ret;
        
}

//' Grouped Vecchia approximation to the Gaussian loglikelihood, zero mean
//'
//' This function returns a grouped version (Guinness, 2018) of Vecchia's (1988) 
//' approximation to the Gaussian
//' loglikelihood. The approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @param covparms A vector of covariance parameters appropriate
//' for the specified covariance function
//' @param covfun_name See \code{\link{GpGp}} for information about covariance
//' functions.
//' @param y vector of response values
//' @param locs matrix of locations. Row \code{i} of \code{locs} specifies the location
//' of element \code{i} of \code{y}, and so the length of \code{y} should equal
//' the number of rows of \code{locs}.
//' @param NNlist A neighbor list object, the output from \code{\link{group_obs}}.
//' @return a list containing
//' \itemize{
//'  \item \code{loglik}: the loglikelihood
//' }
//' @examples
//' n1 <- 20
//' n2 <- 20
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' covparms <- c(2, 0.2, 0.75, 0)
//' y <- fast_Gp_sim(covparms, "matern_isotropic", locs, 50 )
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' NNlist <- group_obs(NNarray)
//' loglik <- vecchia_grouped_meanzero_loglik( covparms, "matern_isotropic", y, locs, NNlist )
//' @export
// [[Rcpp::export]]
List vecchia_grouped_meanzero_loglik( 
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    const NumericMatrix locs,
    List NNlist ){
    
    NumericMatrix X(1,1);
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    //NumericVector betahat;
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );
    //NumericMatrix betainfo;

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    // maybe the synthesize functions should take in an argument that
    // says which compute_pieces function to use
    synthesize_grouped(covparms, covfun_name, locs, NNlist, y, X,
        &ll, &betahat, &grad, &info, &betainfo, false, false 
    );
    
    
    List ret = List::create( Named("loglik") = ll );
    return ret;
        
}
