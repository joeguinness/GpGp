#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace Rcpp;


//' Multiply approximate inverse Cholesky by a vector
//'
//' Vecchia's approximation implies a sparse approximation to the
//' inverse Cholesky factor of the covariance matrix. This function
//' returns the result of multiplying that matrix by a vector.
//' @param Linv Entries of the sparse inverse Cholesky factor,
//' usually the output from \code{\link{vecchia_Linv}}.
//' @param z the vector to be multiplied
//' @inheritParams vecchia_meanzero_loglik
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
//' @inheritParams vecchia_meanzero_loglik
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


//' Multiply transpose of approximate inverse Cholesky by a vector
//'
//' Vecchia's approximation implies a sparse approximation to the
//' inverse Cholesky factor of the covariance matrix. This function
//' returns the result of multiplying the transpose of that matrix by a vector.
//' @param Linv Entries of the sparse inverse Cholesky factor,
//' usually the output from \code{\link{vecchia_Linv}}.
//' @param z the vector to be multiplied
//' @inheritParams vecchia_meanzero_loglik
//' @return the product of the transpose of the 
//' sparse inverse Cholesky factor with a vector
//' @examples
//' n <- 2000
//' locs <- matrix( runif(2*n), n, 2 )
//' covparms <- c(2, 0.2, 0.75, 0.1)
//' NNarray <- find_ordered_nn(locs,20)
//' Linv <- vecchia_Linv( covparms, "matern_isotropic", locs, NNarray )
//' z1 <- rnorm(n)
//' z2 <- Linv_t_mult(Linv, z1, NNarray)
//' @export
// [[Rcpp::export]]
NumericVector Linv_t_mult(NumericMatrix Linv, NumericVector z,
                                  IntegerMatrix NNarray) {

    // return x = t(Linv) * z
    int n = z.length();
    NumericVector x(n);
    for(int j=0;j<n;j++){ x[j] = 0.0; }

    // number of neighbors + 1
    int m = NNarray.ncol();

    // rows 1 though n
    for(int i=0; i<n; i++){
        int bsize = min(i+1,m);
        for(int j=0; j<bsize; j++){
            x( NNarray(i,j) - 1 ) += z( i )*Linv(i,j);
        }
    }

    return x;
}



//' Multiply transpose of approximate Cholesky by a vector
//'
//' Vecchia's approximation implies a sparse approximation to the
//' inverse Cholesky factor of the covariance matrix. This function
//' returns the result of multiplying the transpose of the
//' inverse of that matrix by a vector
//' (i.e. an approximation to the transpose of the Cholesky factor).
//' @param Linv Entries of the sparse inverse Cholesky factor,
//' usually the output from \code{\link{vecchia_Linv}}.
//' @param z the vector to be multiplied
//' @inheritParams vecchia_meanzero_loglik
//' @return the product of the transpose of the Cholesky factor with a vector
//' @examples
//' n <- 2000
//' locs <- matrix( runif(2*n), n, 2 )
//' covparms <- c(2, 0.2, 0.75, 0.1)
//' NNarray <- find_ordered_nn(locs,20)
//' Linv <- vecchia_Linv( covparms, "matern_isotropic", locs, NNarray )
//' z1 <- rnorm(n)
//' z2 <- L_t_mult(Linv, z1, NNarray)
//' @export
// [[Rcpp::export]]
NumericVector L_t_mult(NumericMatrix Linv, NumericVector z,
                               IntegerMatrix NNarray) {

    // return x = L z
    // by solving (L^{-1})x = z
    int n = z.length();
    NumericVector x(n);

    // number of neighbors + 1
    int m = NNarray.ncol();

    // initialize all the entries
    for(int i=0; i<n; i++){ x(i) = z(i)/Linv(i,0); }

    // update entries n-2 through 0
    for(int i=n-1; i>=0; i--){
        int B = min(i+1,m);
        for(int j=1; j<B; j++){
            x( NNarray(i,j) - 1 ) -= Linv(i,j)*x(i)/Linv(NNarray(i,j)-1,0);
        }
    }

    return x;
}
