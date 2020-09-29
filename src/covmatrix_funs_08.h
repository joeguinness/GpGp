
#ifndef COVMATRIX_FUNS_nonstatvar_H
#define COVMATRIX_FUNS_nonstatvar_H

// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include "basis.h"
#include "covmatrix_funs_01.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]





//' Isotropic Matern covariance function, nonstationary variances
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range, smoothness, nugget, <nonstat variance parameters>), 
//' return the square matrix of all pairwise covariances.
//' @param Z A matrix with \code{n} rows and \code{2} columns for spatial
//' locations + \code{p} columns describing spatial basis functions.
//' Each row of locs gives a point in R^2 (two dimensions only!) + the value
//' of \code{p} spatial basis functions.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range, smoothness, nugget, <nonstat variance parameters>).
//' The number of nonstationary variance parameters should equal \code{p}.
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' This covariance function multiplies the isotropic Matern covariance
//' by a nonstationary variance function. The form of the covariance is
//' \deqn{ C(x,y) = exp( \phi(x) + \phi(y) ) M(x,y) }
//' where M(x,y) is the isotropic Matern covariance, and 
//' \deqn{ \phi(x) = c_1 \phi_1(x) + ... + c_p \phi_p(x) }
//' where \eqn{\phi_1,...,\phi_p} are the spatial basis functions
//' contained in the last \code{p} columns of \code{Z}, and 
//' \eqn{c_1,...,c_p} are the nonstationary variance parameters.
// [[Rcpp::export]]
arma::mat matern_nonstat_var(arma::vec covparms, arma::mat Z ){
    
    // this is a 2D covariance function!
    // first two columns of Z are spatial locations
    // rest of columns are values of basis functions at each location
    // log variance is linear in the basis functions
    // covparms(0) = overall variance
    // covparms(1) = isotropic range
    // covparms(2) = smoothness
    // covparms(3) = nugget
    // covparms(4) ... covparms(covparms.n_elem-1) = coefficients
    // in log linear variance function multiplying Z(,2) .... Z(,Z.n_cols-1)
    int dim = 2;
    int n = Z.n_rows;
    int nbasis = Z.n_cols - dim;
    int nisoparm = 4;
    double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*boost::math::tgamma(covparms(2)) );
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate scaled distance
            double d = 0.0;
            for(int j=0; j<dim; j++){
                d += pow( (Z(i1,j) - Z(i2,j))/covparms(1), 2.0 );
            }
            d = pow( d, 0.5 );
            
            // calculate nonstationary variance
            double v = 0.0;
            for(int j=0; j<nbasis; j++){
                v += ( Z(i1, j+dim) + Z(i2, j+dim) ) * covparms( j + nisoparm );
            }
            v = std::exp(v);
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0) * v;
            } else {
                // calculate covariance
                covmat(i2,i1) = normcon*v *
                    pow( d, covparms(2) ) * boost::math::cyl_bessel_k(covparms(2), d);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

//' @describeIn matern_nonstat_var Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_matern_nonstat_var(arma::vec covparms, arma::mat Z ){

    int dim = 2;
    int n = Z.n_rows;
    int nbasis = Z.n_cols - dim;
    int nisoparm = 4;
    //double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*boost::math::tgamma(covparms(2)));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,covparms(2)+eps-1.0)*boost::math::tgamma(covparms(2) + eps));
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        // calculate scaled distance
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( (Z(i1,j) - Z(i2,j))/covparms(1), 2.0 );
        }
        d = pow( d, 0.5 );
        
        // calculate nonstationary variance
        double v = 0.0;
        for(int j=0; j<nbasis; j++){
            v += ( Z(i1, j+dim) + Z(i2, j+dim) ) * covparms( j + nisoparm );
        }
        v = std::exp(v);
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0) * v;
            dcovmat(i2,i1,0) += v;
            for(int j=0; j<nbasis; j++){
                dcovmat(i2,i1,j+nisoparm) = cov*( Z(i1,j+dim) + Z(i2,j+dim) );
            }
        } else {
            cov = normcon * v *
                pow( d, covparms(2) ) *boost::math::cyl_bessel_k(covparms(2), d);
            // variance parameter
            dcovmat(i2,i1,0) += cov/covparms(0);
            // range parameter 
            dcovmat(i2,i1,1) += normcon * v * pow(d,covparms(2))*
                boost::math::cyl_bessel_k(covparms(2)-1.0, d)*d/covparms(1);
            // smoothness parameter (finite differencing)
            dcovmat(i2,i1,2) += 
                ( normconeps*v*pow(d,covparms(2)+eps)*
                boost::math::cyl_bessel_k(covparms(2)+eps, d)- cov )/eps;
            // log linear variance parameters
            for(int j=0; j<nbasis; j++){
                dcovmat(i2,i1,j+nisoparm) = cov*( Z(i1,j+dim) + Z(i2,j+dim) );
            }
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(3);
            dcovmat(i1,i2,3) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i1,i2,j) = dcovmat(i2,i1,j);
            }
        }
    }}

    return dcovmat;
}






//' Isotropic exponential covariance function, nonstationary variances
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range, nugget, <nonstat variance parameters>), 
//' return the square matrix of all pairwise covariances.
//' @param Z A matrix with \code{n} rows and \code{2} columns for spatial
//' locations + \code{p} columns describing spatial basis functions.
//' Each row of locs gives a point in R^2 (two dimensions only!) + the value
//' of \code{p} spatial basis functions.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range, nugget, <nonstat variance parameters>).
//' The number of nonstationary variance parameters should equal \code{p}.
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' This covariance function multiplies the isotropic exponential covariance
//' by a nonstationary variance function. The form of the covariance is
//' \deqn{ C(x,y) = exp( \phi(x) + \phi(y) ) M(x,y) }
//' where M(x,y) is the isotropic exponential covariance, and 
//' \deqn{ \phi(x) = c_1 \phi_1(x) + ... + c_p \phi_p(x) }
//' where \eqn{\phi_1,...,\phi_p} are the spatial basis functions
//' contained in the last \code{p} columns of \code{Z}, and 
//' \eqn{c_1,...,c_p} are the nonstationary variance parameters.
// [[Rcpp::export]]
arma::mat exponential_nonstat_var(arma::vec covparms, arma::mat Z ){
    
    // this is a 2D covariance function!
    // first two columns of Z are spatial locations
    // rest of columns are values of basis functions at each location
    // log variance is linear in the basis functions
    // covparms(0) = overall variance
    // covparms(1) = isotropic range
    // covparms(2) = nugget
    // covparms(3) ... covparms(covparms.n_elem-1) = coefficients
    // in log linear variance function multiplying Z(,2) .... Z(,Z.n_cols-1)
    int dim = 2;
    int n = Z.n_rows;
    int nbasis = Z.n_cols - dim;
    int nisoparm = 3;
    double nugget = covparms( 0 )*covparms( 2 );

    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate scaled distance
            double d = 0.0;
            for(int j=0; j<dim; j++){
                d += pow( (Z(i1,j) - Z(i2,j))/covparms(1), 2.0 );
            }
            d = pow( d, 0.5 );
            
            // calculate nonstationary variance
            double v = 0.0;
            for(int j=0; j<nbasis; j++){
                v += ( Z(i1, j+dim) + Z(i2, j+dim) ) * covparms( j + nisoparm );
            }
            v = std::exp(v);
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0) * v;
            } else {
                // calculate covariance
                covmat(i2,i1) = covparms(0) * v * exp(-d);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

//' @describeIn exponential_nonstat_var Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_exponential_nonstat_var(arma::vec covparms, arma::mat Z ){

    int dim = 2;
    int n = Z.n_rows;
    int nbasis = Z.n_cols - dim;
    int nisoparm = 3;
    //double nugget = covparms( 0 )*covparms( 2 );

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        // calculate scaled distance
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( (Z(i1,j) - Z(i2,j))/covparms(1), 2.0 );
        }
        d = pow( d, 0.5 );
        
        // calculate nonstationary variance
        double v = 0.0;
        for(int j=0; j<nbasis; j++){
            v += ( Z(i1, j+dim) + Z(i2, j+dim) ) * covparms( j + nisoparm );
        }
        v = std::exp(v);
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0) * v;
            dcovmat(i2,i1,0) += v;
            for(int j=0; j<nbasis; j++){
                dcovmat(i2,i1,j+nisoparm) = cov*( Z(i1,j+dim) + Z(i2,j+dim) );
            }
        } else {
            cov = covparms(0) * v * exp(-d);
            // variance parameter
            dcovmat(i2,i1,0) += cov/covparms(0);
            // range parameter
            dcovmat(i2,i1,1) += covparms(0) * v * exp(-d) * d / covparms(1);
            // log linear variance parameters
            for(int j=0; j<nbasis; j++){
                dcovmat(i2,i1,j+nisoparm) = cov*( Z(i1,j+dim) + Z(i2,j+dim) );
            }
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(2);
            dcovmat(i1,i2,2) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i1,i2,j) = dcovmat(i2,i1,j);
            }
        }
    }}

    return dcovmat;
}



#endif
