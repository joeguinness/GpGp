#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "covfuns.h"
#include "vecchiafun.h"

using namespace std;
using namespace Rcpp;


//' Vecchia's approximation to the Gaussian loglikelihood
//' 
//' This function returns Vecchia's (1988) approximation to the Gaussian
//' loglikelihood. The approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @param covparms A vector of covariance parameters appropriate
//' for the specified covariance function
//' @param covfun_name One of "matern_isotropic", "matern_space_time", "matern_sphere", 
//' or "matern_sphere_time".
//' "matern_isotropic" and "matern_sphere" have four covariance parameters, 
//' (variance, range, smoothness, nugget), while "matern_space_time" and 
//' "matern_sphere_time" have five,
//' (variance, spatial range, temporal range, smoothness, nugget). 
//' For more details, see the documentation 
//' for each of the covariance functions by typing, for example, ?matern_isotropic
//' or ?matern_sphere_time.
//' @param y vector of response values
//' @param locs matrix of locations. Row \code{i} of locs specifies the location
//' of element \code{i} of \code{y}, and so the length of \code{y} should equal
//' the number of rows of \code{locs}.
//' @param NNarray A matrix of indices, usually the output from \code{\link{find_ordered_nn}}. 
//' Row \code{i} contains the indices
//' of the observations that observation \code{i} conditions on. By convention,
//' the first element of row \code{i} is \code{i}.
//' @return the Gaussian loglikelihood
//' @examples
//' n1 <- 40
//' n2 <- 40
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' covparms <- c(2, 0.2, 0.75, 0)
//' y <- fast_Gp_sim(covparms, "matern_isotropic", locs, 50 ) 
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' loglik <- vecchia_loglik( covparms, "matern_isotropic", y, locs, NNarray )
//' @export
// [[Rcpp::export]]
NumericVector vecchia_loglik(NumericVector covparms, StringVector covfun_name,
                                  NumericVector y,
                                  const NumericMatrix locs, IntegerMatrix NNarray) {

    NumericVector ll(1);        // loglikelihood to be returned
    NumericMatrix Linv(1,1);    // Linv not to be returned
    vecchia(covparms, covfun_name, locs, NNarray, y, &Linv, &ll, 1);
    return ll;
}

//' Inverse Cholesky factor implied by Vecchia's approximation
//' 
//' This function returns the entries of the sparse approximation to 
//' the Cholesky factor implied by Vecchia's (1988) approximation to the Gaussian
//' loglikelihood. The approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @inheritParams vecchia_loglik
//' @return matrix containing entries of sparse approximation to inverse Cholesky
//' @examples
//' n1 <- 40
//' n2 <- 40
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' covparms <- c(2, 0.2, 0.75, 0)
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' Linv <- vecchia_Linv( covparms, "matern_isotropic", locs, NNarray )
//' @export
// [[Rcpp::export]]
NumericMatrix vecchia_Linv(NumericVector covparms, StringVector covfun_name,
                            NumericMatrix locs, IntegerMatrix NNarray) {

    NumericVector y(NNarray.nrow());
    NumericVector ll(1);        // loglikelihood not to be returned
    NumericMatrix Linv(NNarray.nrow() , NNarray.ncol());    // Linv to be returned
    vecchia(covparms, covfun_name, locs, NNarray, y, &Linv, &ll, 2);
    return Linv;
}


//' Multiply approximate inverse Cholesky by a vector
//' 
//' Vecchia's approximation implies a sparse approximation to the 
//' inverse Cholesky factor of the covariance matrix. This function
//' returns the result of multiplying that matrix by a vector.
//' @param Linv Entries of the sparse inverse Cholesky factor,
//' usually the output from \code{\link{vecchia_Linv}}.
//' @param z the vector to be multiplied
//' @inheritParams vecchia_loglik
//' @return the product of the sparse inverse Cholesky factor with a vector
//' @examples
//' n <- 2000
//' locs <- matrix( runif(2*n), n, 2 )
//' covparms <- c(2, 0.2, 0.75, 0.1)
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' Linv <- vecchia_Linv( covparms, "matern_isotropic", locs, NNarray )
//' z1 <- rnorm(n)
//' y <- fast_Gp_sim_Linv(Linv,NNarray,z1)
//' z2 <- Linv_mult(Linv, y, NNarray)
//' print( sum( (z1-z2)^2 ) )
//' @export
// [[Rcpp::export]]
NumericVector Linv_mult(NumericMatrix Linv, NumericVector z,
                                  IntegerMatrix NNarray) {

    // return x = Linv * z
    int i;
    int j;

    int n = z.length();
    NumericVector x(n);
    for(j=0;j<n;j++){ x[j] = 0.0; }

    // number of neighbors + 1
    int m = NNarray.ncol();

    // rows 1 though n
    for(i=0; i<n; i++){
        int bsize = min(i+1,m);
        for(j=0; j<bsize; j++){
            x( i ) += z( NNarray(i,j) - 1 )*Linv(i,j);
        }
    }

    return x;
}


//' Multiply approximate Cholesky by a vector
//' 
//' Vecchia's approximation implies a sparse approximation to the 
//' inverse Cholesky factor of the covariance matrix. This function
//' returns the result of multiplying the inverse of that matrix by a vector 
//' (i.e. an approximation to the Cholesky factor).
//' @param Linv Entries of the sparse inverse Cholesky factor,
//' usually the output from \code{\link{vecchia_Linv}}.
//' @param z the vector to be multiplied
//' @inheritParams vecchia_loglik
//' @return the product of the Cholesky factor with a vector
//' @examples
//' n <- 2000
//' locs <- matrix( runif(2*n), n, 2 )
//' covparms <- c(2, 0.2, 0.75, 0.1)
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' Linv <- vecchia_Linv( covparms, "matern_isotropic", locs, NNarray )
//' z <- rnorm(n)
//' y1 <- fast_Gp_sim_Linv(Linv,NNarray,z)
//' y2 <- L_mult(Linv, z, NNarray)
//' print( sum( (y1-y2)^2 ) )
//' @export
// [[Rcpp::export]]
NumericVector L_mult(NumericMatrix Linv, NumericVector z,
                               IntegerMatrix NNarray) {

    // return x = L z
    // by solving (L^{-1})x = z
    int i;
    int j;
    int B;

    int n = z.length();
    NumericVector x(n);

    // number of neighbors + 1
    int m = NNarray.ncol();

    // get entry 0
    x(0) = z(0)/Linv(0,0);

    // get entries 1 through n
    for(i=1; i<n; i++){
        B = min(i+1,m);
        x(i) = z(i);
        for(j=1; j<B; j++){
            x(i) -= Linv(i,j)*x( NNarray(i,j) - 1 );
        }
        x(i) = x(i)/Linv(i,0);
    }

    return x;
}



