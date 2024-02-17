#ifndef COVMATRIX_FUNS_03_H
#define COVMATRIX_FUNS_03_H

// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]


//' Matern covariance function, different range parameter for each dimension
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range_1, ..., range_d, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs is a point in R^d.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range_1, ..., range_d, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range_1, ..., range_d, smoothness, nugget).
//' The covariance function is parameterized as
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) (|| D^{-1}(x - y) || )^\nu K_\nu(|| D^{-1}(x - y) || ) }
//' where D is a diagonal matrix with (range_1, ..., range_d) on the diagonals.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern_scaledim(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;

    if( covparms.n_elem - 3 != dim ){
        stop("length of covparms does not match dim of locs");
    }
            
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( dim + 2 );
    double smooth = covparms(dim+1);
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*boost::math::tgamma(smooth));
    
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1+j);
        }
    }
    // try simply plugging this into the matern_isotropic function
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){ for(int i2 = 0; i2 <= i1; i2++){
            
        // calculate distance
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        if( d == 0.0 ){
            covmat(i2,i1) = covparms(0);
        } else {
            // calculate covariance            
            covmat(i2,i1) = normcon*
                pow( d, smooth )* boost::math::cyl_bessel_k(smooth, d);
        }
        // add nugget
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
        // fill in opposite entry
        else { covmat(i1,i2) = covmat(i2,i1); }
    }}
    return covmat;
}

//' @describeIn matern_scaledim Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_matern_scaledim(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;
    if( covparms.n_elem - 3 != dim ){
        stop("length of covparms does not match dim of locs");
    }
            
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( dim + 2 );
    double smooth = covparms(dim+1);
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*boost::math::tgamma(smooth));
    
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1+j);
        }
    }

    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,smooth+eps-1.0)*boost::math::tgamma(smooth+eps));
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
        } else {
            cov = normcon*pow( d, smooth )* boost::math::cyl_bessel_k(smooth, d);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // range parameters
            for(int j=0; j<dim; j++){
                double dj2 = pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
                dcovmat(i1,i2,j+1) += normcon*pow( d, smooth - 1.0)*
                    boost::math::cyl_bessel_k(smooth - 1.0, d)*dj2/covparms(j+1);
            }
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,dim+1) += 
                ( normconeps*pow(d,smooth+eps)* boost::math::cyl_bessel_k(smooth + eps, d) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(dim+2);
            dcovmat(i1,i2,dim+2) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}
    return dcovmat;
}




//' Exponential covariance function, different range parameter for each dimension
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
//' \deqn{ M(x,y) = \sigma^2 exp( - || D^{-1}(x - y) || ) }
//' where D is a diagonal matrix with (range_1, ..., range_d) on the diagonals.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat exponential_scaledim(arma::vec covparms, arma::mat locs ){

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
    for(int i1 = 0; i1 < n; i1++){ for(int i2 = 0; i2 <= i1; i2++){
            
        // calculate distance
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        if( d == 0.0 ){
            covmat(i2,i1) = covparms(0);
        } else {
            // calculate covariance            
            covmat(i2,i1) = covparms(0)*exp(-d);
        }
        // add nugget
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
        // fill in opposite entry
        else { covmat(i1,i2) = covmat(i2,i1); }
    }}
    return covmat;
}

//' @describeIn exponential_scaledim Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_exponential_scaledim(arma::vec covparms, arma::mat locs ){

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
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
        } else {
            cov = covparms(0)*exp(-d);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // range parameters
            for(int j=0; j<dim; j++){
                double dj2 = pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
                dcovmat(i1,i2,j+1) += cov/d*dj2/covparms(j+1);
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







//' Spatial-Temporal Matern covariance function
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range_1, range_2, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d+1} columns.
//' Each row of locs is a point in R^(d+1). The first \code{d} columns
//' should contain the spatial coordinates. The last column contains the times.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range_1, range_2, smoothness, nugget). range_1 is the
//' spatial range, and range_2 is the temporal range.
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range_1, range_2, smoothness, nugget).
//' The covariance function is parameterized as
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) (|| D^{-1}(x - y) || )^\nu K_\nu(|| D^{-1}(x - y) || ) }
//' where D is a diagonal matrix with (range_1, ..., range_1, range_2) on the diagonals.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern_spacetime(arma::vec covparms, arma::mat locs ){
    
    // number of spatial dimensions
    int dim = locs.n_cols - 1;
    int n = locs.n_rows;

    // create scaled locations
    arma::mat locs_scaled(n,dim+1);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    for(int i=0; i<n; i++){
        locs_scaled(i,dim) = locs(i,dim)/covparms(2);
    }
    
    arma::vec newparms(4);
    newparms(0) = covparms(0);
    newparms(1) = 1.0;
    newparms(2) = covparms(3);
    newparms(3) = covparms(4);
    arma::mat covmat = matern_isotropic( newparms, locs_scaled );

    return covmat;
}

//' @describeIn matern_spacetime Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_matern_spacetime(arma::vec covparms, arma::mat locs ){
    
    // number of spatial dimensions
    int dim = locs.n_cols - 1;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 4 );
    double smooth = covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*boost::math::tgamma(smooth));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,smooth+eps-1.0)*boost::math::tgamma(smooth+eps));

    // create scaled locations
    arma::mat locs_scaled(n,dim+1);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    for(int i=0; i<n; i++){
        locs_scaled(i,dim) = locs(i,dim)/covparms(2);
    }
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        double d = 0.0;
        for(int j=0; j<dim+1; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
        } else {
            cov = normcon*pow( d, smooth )*boost::math::cyl_bessel_k(smooth, d);

            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);

            // range parameters
            double dj2 = 0.0; 
            for(int j=0; j<dim; j++){
                dj2 += pow(locs_scaled(i1,j)-locs_scaled(i2,j),2.0);
            }
            dcovmat(i1,i2,1) += normcon*pow( d, smooth - 1.0)*
                boost::math::cyl_bessel_k(smooth - 1.0, d)*dj2/covparms(1);
            dj2 = pow(locs_scaled(i1,dim)-locs_scaled(i2,dim),2.0);
            dcovmat(i1,i2,2) += normcon*pow( d, smooth - 1.0)*
                boost::math::cyl_bessel_k(smooth - 1.0, d)*dj2/covparms(2);
            
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,3) += 
                ( normconeps*pow(d,smooth+eps)* boost::math::cyl_bessel_k(smooth + eps, d) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(4);
            dcovmat(i1,i2,4) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}
    return dcovmat;
}




//' Spatial-Temporal exponential covariance function
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range_1, range_2, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d+1} columns.
//' Each row of locs is a point in R^(d+1). The first \code{d} columns
//' should contain the spatial coordinates. The last column contains the times.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range_1, range_2, nugget). range_1 is the
//' spatial range, and range_2 is the temporal range.
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range_1, range_2, nugget).
//' The covariance function is parameterized as
//' \deqn{ M(x,y) = \sigma^2 exp( - || D^{-1}(x - y) || ) }
//' where D is a diagonal matrix with (range_1, ..., range_1, range_2) on the diagonals.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat exponential_spacetime(arma::vec covparms, arma::mat locs ){
    
    // number of spatial dimensions
    int dim = locs.n_cols - 1;
    int n = locs.n_rows;

    // create scaled locations
    arma::mat locs_scaled(n,dim+1);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    for(int i=0; i<n; i++){
        locs_scaled(i,dim) = locs(i,dim)/covparms(2);
    }
    
    arma::vec newparms(4);
    newparms(0) = covparms(0);
    newparms(1) = 1.0;
    newparms(2) = covparms(3);
    arma::mat covmat = exponential_isotropic( newparms, locs_scaled );

    return covmat;
}

//' @describeIn exponential_spacetime Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_exponential_spacetime(arma::vec covparms, arma::mat locs ){
    
    // number of spatial dimensions
    int dim = locs.n_cols - 1;
    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 3 );

    // create scaled locations
    arma::mat locs_scaled(n,dim+1);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    for(int i=0; i<n; i++){
        locs_scaled(i,dim) = locs(i,dim)/covparms(2);
    }
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        double d = 0.0;
        for(int j=0; j<dim+1; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
        } else {
            cov = covparms(0)*exp(-d);

            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);

            // range parameters
            double dj2 = 0.0; 
            for(int j=0; j<dim; j++){
                dj2 += pow(locs_scaled(i1,j)-locs_scaled(i2,j),2.0);
            }
            dcovmat(i1,i2,1) += covparms(0)*exp(-d)/d*dj2/covparms(1);
            dj2 = pow(locs_scaled(i1,dim)-locs_scaled(i2,dim),2.0);
            dcovmat(i1,i2,2) += covparms(0)*exp(-d)/d*dj2/covparms(2);
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(3);
            dcovmat(i1,i2,3) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}
    return dcovmat;
}


#endif
