#ifndef COVMATRIX_FUNS_12_H
#define COVMATRIX_FUNS_12_H



// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]


//' Isotropic Matern covariance function, smoothness = 4.5
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs is a point in R^d.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range, nugget)
//' = \eqn{(\sigma^2,\alpha,\tau^2)}, and the covariance function is parameterized
//' as
//' \deqn{ M(x,y) = \sigma^2 ( \sum_{j=0}^4 c_j || x - y ||^j/ \alpha^j ) exp( - || x - y ||/ \alpha )}
//' where c_0 = 1, c_1 = 1, c_2 = 3/7, c_3 = 2/21, c_4 = 1/105. 
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern45_isotropic(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 2 );
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    // calculate covariances
    arma::mat covmat(n,n);
    double a0 = 1.0;
    double a1 = 1.0;
    double a2 = 3.0/7.0;
    double a3 = 2.0/21.0;
    double a4 = 1.0/105.0;
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        // calculate distance
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        double pdc = a0 + a1*d + a2*d*d + a3*pow(d,3.0) + a4*pow(d,4.0);
        if( d == 0.0 ){
            covmat(i2,i1) = covparms(0);
        } else {
            // calculate covariance
            covmat(i2,i1) = covparms(0)*pdc*std::exp(-d);
        }
        // add nugget
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
        // fill in opposite entry
        else { covmat(i1,i2) = covmat(i2,i1); }
    }}

    return covmat;
}

//' @describeIn matern35_isotropic Derivatives of isotropic
//' matern covariance function with smoothness 3.5
// [[Rcpp::export]]
arma::cube d_matern45_isotropic(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;
    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 2 );
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    double a0 = 1.0;
    double a1 = 1.0;
    double a2 = 3.0/7.0;
    double a3 = 2.0/21.0;
    double a4 = 1.0/105.0;
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        double dpc = -d/covparms(1);
        double pdc = a0 + a1*d + a2*d*d + a3*pow(d,3.0) + a4*pow(d,4.0);
        double ppdc = a1 + 2*a2*d + 3*a3*d*d + 4*a4*pow(d,3.0);
            
        dcovmat(i1,i2,0) += pdc*std::exp(-d);
        dcovmat(i1,i2,1) += covparms(0)*std::exp(-d)*dpc*(ppdc - pdc);
        
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(2);
            dcovmat(i1,i2,2) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}
    return dcovmat;
}



//' Matern covariance function, smoothess = 3.5, different range parameter for each dimension
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range_1, ..., range_d, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs is a point in R^d.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range_1, ..., range_d, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range_1, ..., range_d, nugget).
//' The covariance function is parameterized as
//' \deqn{ M(x,y) = \sigma^2 ( \sum_{j=0}^4 c_j || D^{-1}(x - y) ||^j  ) exp( - || D^{-1}(x - y) || ) }
//' where c_0 = 1, c_1 = 1, c_2 = 3/7, c_3 = 2/21, c_4 = 1/105. 
//' where D is a diagonal matrix with (range_1, ..., range_d) on the diagonals.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern45_scaledim(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;

    if( covparms.n_elem - 2 != dim ){
        stop("length of covparms does not match dim of locs");
    }
            
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( dim + 1 );

    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1+j);
        }
    }

    // calculate covariances
    arma::mat covmat(n,n);
    double a0 = 1.0;
    double a1 = 1.0;
    double a2 = 3.0/7.0;
    double a3 = 2.0/21.0;
    double a4 = 1.0/105.0;
    for(int i1 = 0; i1 < n; i1++){ for(int i2 = 0; i2 <= i1; i2++){
            
        // calculate distance
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        double pdc = a0 + a1*d + a2*d*d + a3*pow(d,3.0) + a4*pow(d,4.0);
        
        if( d == 0.0 ){
            covmat(i2,i1) = covparms(0);
        } else {
            // calculate covariance            
            covmat(i2,i1) = covparms(0)*pdc*exp(-d);
        }
        // add nugget
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
        // fill in opposite entry
        else { covmat(i1,i2) = covmat(i2,i1); }
    }}
    return covmat;
}

//' @describeIn matern35_scaledim Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_matern45_scaledim(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;
    if( covparms.n_elem - 2 != dim ){
        stop("length of covparms does not match dim of locs");
    }
            
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( dim + 1 );

    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1+j);
        }
    }

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    double a0 = 1.0;
    double a1 = 1.0;
    double a2 = 3.0/7.0;
    double a3 = 2.0/21.0;
    double a4 = 1.0/105.0;
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        double pdc = a0 + a1*d + a2*d*d + a3*pow(d,3.0) + a4*pow(d,4.0);
        double ppdc = a1 + 2*a2*d + 3*a3*d*d + 4*a4*pow(d,3.0);

        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
        } else {
            cov = covparms(0)*pdc*exp(-d);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // range parameters
            for(int j=0; j<dim; j++){
                double dj2 = pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
                double dpc = -dj2/covparms(j+1)/d;
                dcovmat(i1,i2,j+1) += covparms(0)*exp(-d)*dpc*( ppdc - pdc );
            }
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(dim+1);
            dcovmat(i1,i2,dim+1) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}
    return dcovmat;
}





#endif

