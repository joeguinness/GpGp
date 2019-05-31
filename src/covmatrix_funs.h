#ifndef ARMACOVMATRIX_FUNS_H
#define ARMACOVMATRIX_FUNS_H



// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "basis.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]


struct covfun_t {
    arma::mat    (*p_covfun)(NumericVector, NumericMatrix);
    arma::cube (*p_d_covfun)(NumericVector, NumericMatrix);
} ;


//' Isotropic Matern covariance function
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs gives a point in R^d.
//' @param covparms A vector giving positive-valued covariance parameters
//' in the form (variance, range, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range, smoothness, nugget)
//' = \eqn{(\sigma^2,\alpha,\nu,\tau^2)}, and the covariance function is parameterized
//' as
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) (|| x - y ||/\alpha )^\nu K_\nu(|| x - y ||/\alpha ) }
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. The reason for this choice
//' is for simpler profiling of \eqn{ \sigma^2 }.
// [[Rcpp::export]]
arma::mat matern_isotropic(NumericVector covparms, NumericMatrix locs ){

    int dim = locs.ncol();
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*Rf_gammafn(covparms(2)));
    
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
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
                    pow( d, covparms(2) )*Rf_bessel_k(d,covparms(2),1.0);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}



// [[Rcpp::export]]
arma::cube d_matern_isotropic(NumericVector covparms, NumericMatrix locs ){

    int dim = locs.ncol();
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*Rf_gammafn(covparms(2)));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,covparms(2)+eps-1.0)*Rf_gammafn(covparms(2)+eps));
    
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.length(), fill::zeros);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
            dcovmat(i1,i2,1) += 0.0;
            dcovmat(i1,i2,2) += 0.0;
        } else {
            cov = normcon*pow( d, covparms(2) )*Rf_bessel_k(d,covparms(2),1.0);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // range parameter
            dcovmat(i1,i2,1) += normcon*pow(d,covparms(2))*
                Rf_bessel_k(d,covparms(2)-1.0,1.0)*d/covparms(1);
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,2) += 
                ( normconeps*pow(d,covparms(2)+eps)*Rf_bessel_k(d,covparms(2)+eps,1.0) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(3);
            dcovmat(i1,i2,3) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.length(); j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}




//' Isotropic exponential covariance function
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs gives a point in R^d.
//' @param covparms A vector giving positive-valued covariance parameters
//' in the form (variance, range, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range, nugget)
//' = \eqn{(\sigma^2,\alpha,\tau^2)}, and the covariance function is parameterized
//' as
//' \deqn{ M(x,y) = \sigma^2 \exp( - || x - y ||/ \alpha }
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. The reason for this choice
//' is for simpler profiling of \eqn{ \sigma^2 }.
// [[Rcpp::export]]
arma::mat exponential_isotropic(NumericVector covparms, NumericMatrix locs ){

    int dim = locs.ncol();
    int n = locs.nrow();
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
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
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
            covmat(i2,i1) = covparms(0)*std::exp( -d );
        }
        // add nugget
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
        // fill in opposite entry
        else { covmat(i1,i2) = covmat(i2,i1); }
    }}

    return covmat;
}

// [[Rcpp::export]]
arma::cube d_exponential_isotropic(NumericVector covparms, NumericMatrix locs ){

    int dim = locs.ncol();
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 2 );
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.length(), fill::zeros);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        dcovmat(i1,i2,0) += std::exp(-d);
        dcovmat(i1,i2,1) += covparms(0)*std::exp(-d)*d/covparms(1);
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(2);
            dcovmat(i1,i2,2) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.length(); j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}




// [[Rcpp::export]]
arma::mat matern_scaledim(NumericVector covparms, NumericMatrix locs ){

    int dim = locs.ncol();

    if( covparms.length() - 3 != dim ){
        stop("length of covparms does not match dim of locs");
    }
            
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( dim + 2 );
    double smooth = covparms(dim+1);
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*Rf_gammafn(smooth));
    
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
                pow( d, smooth )*Rf_bessel_k( d, smooth, 1.0 );
        }
        // add nugget
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
        // fill in opposite entry
        else { covmat(i1,i2) = covmat(i2,i1); }
    }}
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_scaledim(NumericVector covparms, NumericMatrix locs ){

    int dim = locs.ncol();
    if( covparms.length() - 3 != dim ){
        stop("length of covparms does not match dim of locs");
    }
            
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( dim + 2 );
    double smooth = covparms(dim+1);
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*Rf_gammafn(smooth));
    
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1+j);
        }
    }

    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,smooth+eps-1.0)*Rf_gammafn(smooth+eps));
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.length(), fill::zeros);
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
            cov = normcon*pow( d, smooth )*Rf_bessel_k(d,smooth,1.0);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // range parameters
            for(int j=0; j<dim; j++){
                double dj2 = pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
                dcovmat(i1,i2,j+1) += normcon*pow( d, smooth - 1.0)*
                    Rf_bessel_k(d,smooth-1.0,1.0)*dj2/covparms(j+1);
            }
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,dim+1) += 
                ( normconeps*pow(d,smooth+eps)*Rf_bessel_k(d,smooth+eps,1.0) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(dim+2);
            dcovmat(i1,i2,dim+2) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.length(); j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}
    return dcovmat;
}




// [[Rcpp::export]]
arma::mat matern_anisotropic2D(NumericVector covparms, NumericMatrix locs ){
    
    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = smoothness
    // covparms(5) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget
    
    //int dim = locs.ncol();
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 5 );
    double normcon = covparms(0)/(pow(2.0,covparms(4)-1.0)*Rf_gammafn(covparms(4)));
    
    double b0 = covparms(1)*covparms(1);
    double b1 = covparms(2)*covparms(2)+covparms(3)*covparms(3);
    double b2 = covparms(2)*covparms(3);
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate rescaled distance
            double h0 = locs(i1,0) - locs(i2,0);
            double h1 = locs(i1,1) - locs(i2,1);
            double d = h0*h0*b0 + h1*h1*b1 + 2*h0*h1*b2;
            d = pow( d, 0.5 );
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0);
            } else {
                // calculate covariance            
                covmat(i2,i1) = normcon*
                    pow( d, covparms(4) )*Rf_bessel_k(d,covparms(4),1.0);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_anisotropic2D(NumericVector covparms, NumericMatrix locs ){

    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = smoothness
    // covparms(5) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget

    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 5 );
    double normcon = covparms(0)/(pow(2.0,covparms(4)-1.0)*Rf_gammafn(covparms(4)));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,covparms(4)+eps-1.0)*Rf_gammafn(covparms(4)+eps));
    
    double b0 = covparms(1)*covparms(1);
    double b1 = covparms(2)*covparms(2)+covparms(3)*covparms(3);
    double b2 = covparms(2)*covparms(3);

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.length(), fill::zeros);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        
        // calculate rescaled distance
        double h0 = locs(i1,0) - locs(i2,0);
        double h1 = locs(i1,1) - locs(i2,1);
        double d = h0*h0*b0 + h1*h1*b1 + 2*h0*h1*b2;
        d = pow( d, 0.5 );
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
        } else {
            cov = normcon*pow( d, covparms(4) )*Rf_bessel_k(d,covparms(4),1.0);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // cholesky parameters
            double cov_nu_m1 = normcon*pow(d,covparms(4)-1.0)*
                Rf_bessel_k(d,covparms(4)-1.0,1.0);  
            dcovmat(i1,i2,1) -= cov_nu_m1*(h0*h0*covparms(1));
            dcovmat(i1,i2,2) -= cov_nu_m1*(h1*h1*covparms(2) + h0*h1*covparms(3));
            dcovmat(i1,i2,3) -= cov_nu_m1*(h1*h1*covparms(3) + h0*h1*covparms(2));
            
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,4) += 
                ( normconeps*pow(d,covparms(4)+eps)*Rf_bessel_k(d,covparms(4)+eps,1.0) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(5);
            dcovmat(i1,i2,5) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.length(); j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}




// [[Rcpp::export]]
arma::mat matern_anisotropic3D(NumericVector covparms, NumericMatrix locs ){
    
    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = L20
    // covparms(5) = L21
    // covparms(6) = L22
    // covparms(7) = smoothness
    // covparms(8) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget
    
    //int dim = locs.ncol();
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 8 );
    double smooth = covparms( 7 );
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*Rf_gammafn(smooth));
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate rescaled distance
            double h0 = locs(i1,0) - locs(i2,0);
            double h1 = locs(i1,1) - locs(i2,1);
            double h2 = locs(i1,2) - locs(i2,2);
            // 3 is hard coded here
            double d = 0.0;
            d += pow( covparms(1)*h0, 2 );
            d += pow( covparms(2)*h0 + covparms(3)*h1, 2 );
            d += pow( covparms(4)*h0 + covparms(5)*h1 + covparms(6)*h2, 2 );
            d = pow( d, 0.5 );
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0);
            } else {
                // calculate covariance            
                covmat(i2,i1) = normcon*
                    pow( d, smooth )*Rf_bessel_k(d,smooth,1.0);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_anisotropic3D(NumericVector covparms, NumericMatrix locs ){

    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = L20
    // covparms(5) = L21
    // covparms(6) = L22
    // covparms(7) = smoothness
    // covparms(8) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget

    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 8 );
    double smooth = covparms( 7 );
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*Rf_gammafn(smooth));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,smooth+eps-1.0)*Rf_gammafn(smooth+eps));
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.length(), fill::zeros);
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        // calculate rescaled distance
        double h0 = locs(i1,0) - locs(i2,0);
        double h1 = locs(i1,1) - locs(i2,1);
        double h2 = locs(i1,2) - locs(i2,2);
        // 2 spatial + 1 time dimen is hard coded here
        double d = 0.0;
        d += pow( covparms(1)*h0, 2 );
        d += pow( covparms(2)*h0 + covparms(3)*h1, 2 );
        d += pow( covparms(4)*h0 + covparms(5)*h1 + covparms(6)*h2, 2 );
        d = pow( d, 0.5 );
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
        } else {
            cov = normcon*pow( d, smooth )*Rf_bessel_k( d, smooth, 1.0 );
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // cholesky parameters
            double cov_nu_m1 = normcon*pow( d, smooth - 1.0 )*
                Rf_bessel_k( d, smooth - 1.0, 1.0 );  
            double Limhm = covparms(1)*h0;
                dcovmat(i1,i2,1) = -cov_nu_m1*Limhm*h0;
            Limhm = covparms(2)*h0 + covparms(3)*h1;
                dcovmat(i1,i2,2) = -cov_nu_m1*Limhm*h0;
                dcovmat(i1,i2,3) = -cov_nu_m1*Limhm*h1;
            Limhm = covparms(4)*h0 + covparms(5)*h1 + covparms(6)*h2;
                dcovmat(i1,i2,4) = -cov_nu_m1*Limhm*h0;
                dcovmat(i1,i2,5) = -cov_nu_m1*Limhm*h1;
                dcovmat(i1,i2,6) = -cov_nu_m1*Limhm*h2;
    
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,7) += 
                ( normconeps*pow(d,smooth+eps)*Rf_bessel_k(d,smooth+eps,1.0) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(8);
            dcovmat(i1,i2,8) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.length(); j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}




// [[Rcpp::export]]
arma::mat exponential_anisotropic3D(NumericVector covparms, NumericMatrix locs ){
    
    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = L20
    // covparms(5) = L21
    // covparms(6) = L22
    // covparms(7) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget
    
    //int dim = locs.ncol();
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 7 );

    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate rescaled distance
            double h0 = locs(i1,0) - locs(i2,0);
            double h1 = locs(i1,1) - locs(i2,1);
            double h2 = locs(i1,2) - locs(i2,2);
            // 3 dims is hard coded here
            double d = 0.0;
            d += pow( covparms(1)*h0, 2 );
            d += pow( covparms(2)*h0 + covparms(3)*h1, 2 );
            d += pow( covparms(4)*h0 + covparms(5)*h1 + covparms(6)*h2, 2 );
            d = pow( d, 0.5 );
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0);
            } else {
                // calculate covariance            
                covmat(i2,i1) = covparms(0)*std::exp( -d );
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_exponential_anisotropic3D(NumericVector covparms, NumericMatrix locs ){

    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = L20
    // covparms(5) = L21
    // covparms(6) = L22
    // covparms(7) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget

    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 7 );

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.length(), fill::zeros);
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        // calculate rescaled distance
        double h0 = locs(i1,0) - locs(i2,0);
        double h1 = locs(i1,1) - locs(i2,1);
        double h2 = locs(i1,2) - locs(i2,2);
        // 2 spatial + 1 time dimen is hard coded here
        double d = 0.0;
        d += pow( covparms(1)*h0, 2 );
        d += pow( covparms(2)*h0 + covparms(3)*h1, 2 );
        d += pow( covparms(4)*h0 + covparms(5)*h1 + covparms(6)*h2, 2 );
        d = pow( d, 0.5 );
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
        } else {
            cov = covparms(0)*std::exp( -d );
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // cholesky parameters
            double dcov = -covparms(0)*exp(-d)/d;
            double Limhm = covparms(1)*h0;
                dcovmat(i1,i2,1) = dcov*Limhm*h0;
            Limhm = covparms(2)*h0 + covparms(3)*h1;
                dcovmat(i1,i2,2) = dcov*Limhm*h0;
                dcovmat(i1,i2,3) = dcov*Limhm*h1;
            Limhm = covparms(4)*h0 + covparms(5)*h1 + covparms(6)*h2;
                dcovmat(i1,i2,4) = dcov*Limhm*h0;
                dcovmat(i1,i2,5) = dcov*Limhm*h1;
                dcovmat(i1,i2,6) = dcov*Limhm*h2;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(7);
            dcovmat(i1,i2,7) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.length(); j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}





// [[Rcpp::export]]
arma::mat matern_spacetime(NumericVector covparms, NumericMatrix locs ){
    
    // number of spatial dimensions
    int dim = locs.ncol() - 1;
    int n = locs.nrow();

    // create scaled locations
    NumericMatrix locs_scaled(n,dim+1);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    for(int i=0; i<n; i++){
        locs_scaled(i,dim) = locs(i,dim)/covparms(2);
    }
    
    NumericVector newparms(4);
    newparms(0) = covparms(0);
    newparms(1) = 1.0;
    newparms(2) = covparms(3);
    newparms(3) = covparms(4);
    arma::mat covmat = matern_isotropic( newparms, locs_scaled );

    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_spacetime(NumericVector covparms, NumericMatrix locs ){
    
    // number of spatial dimensions
    int dim = locs.ncol() - 1;
    int n = locs.nrow();
    double nugget = covparms( 0 )*covparms( 4 );
    double smooth = covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*Rf_gammafn(smooth));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,smooth+eps-1.0)*Rf_gammafn(smooth+eps));

    // create scaled locations
    NumericMatrix locs_scaled(n,dim+1);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    for(int i=0; i<n; i++){
        locs_scaled(i,dim) = locs(i,dim)/covparms(2);
    }
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.length(), fill::zeros);
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
            cov = normcon*pow( d, smooth )*Rf_bessel_k(d,smooth,1.0);

            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);

            // range parameters
            double dj2 = 0.0; 
            for(int j=0; j<dim; j++){
                dj2 += pow(locs_scaled(i1,j)-locs_scaled(i2,j),2.0);
            }
            dcovmat(i1,i2,1) += normcon*pow( d, smooth - 1.0)*
                Rf_bessel_k(d,smooth-1.0,1.0)*dj2/covparms(1);
            dj2 = pow(locs_scaled(i1,dim)-locs_scaled(i2,dim),2.0);
            dcovmat(i1,i2,2) += normcon*pow( d, smooth - 1.0)*
                Rf_bessel_k(d,smooth-1.0,1.0)*dj2/covparms(2);
            
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,3) += 
                ( normconeps*pow(d,smooth+eps)*Rf_bessel_k(d,smooth+eps,1.0) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(4);
            dcovmat(i1,i2,4) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.length(); j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}
    return dcovmat;
}





// [[Rcpp::export]]
arma::mat matern_nonstat_var(NumericVector covparms, NumericMatrix Z ){
    
    // this is a 2D covariance function!
    // first two columns of Z are spatial locations
    // rest of columns are values of basis functions at each location
    // log variance is linear in the basis functions
    // covparms(0) = overall variance
    // covparms(1) = isotropic range
    // covparms(2) = smoothness
    // covparms(3) = nugget
    // covparms(4) ... covparms(covparms.length()-1) = coefficients
    // in log linear variance function multiplying Z(,2) .... Z(,Z.ncol()-1)
    int dim = 2;
    int n = Z.nrow();
    int nbasis = Z.ncol() - dim;
    int nisoparm = 4;
    double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*Rf_gammafn(covparms(2)));
    
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
                    pow( d, covparms(2) ) * Rf_bessel_k( d, covparms(2), 1.0 );
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_nonstat_var(NumericVector covparms, NumericMatrix Z ){

    int dim = 2;
    int n = Z.nrow();
    int nbasis = Z.ncol() - dim;
    int nisoparm = 4;
    double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*Rf_gammafn(covparms(2)));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,covparms(2)+eps-1.0)*Rf_gammafn(covparms(2)+eps));
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.length(), fill::zeros);
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
                pow( d, covparms(2) ) * Rf_bessel_k( d, covparms(2), 1.0 );
            // variance parameter
            dcovmat(i2,i1,0) += cov/covparms(0);
            // range parameter
            dcovmat(i2,i1,1) += normcon * v * pow(d,covparms(2))*
                Rf_bessel_k(d,covparms(2)-1.0,1.0)*d/covparms(1);
            // smoothness parameter (finite differencing)
            dcovmat(i2,i1,2) += 
                ( normconeps*v*pow(d,covparms(2)+eps)*
                  Rf_bessel_k(d,covparms(2)+eps,1.0) - cov )/eps;
            // log linear variance parameters
            for(int j=0; j<nbasis; j++){
                dcovmat(i2,i1,j+nisoparm) = cov*( Z(i1,j+dim) + Z(i2,j+dim) );
            }
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(3);
            dcovmat(i1,i2,3) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.length(); j++){
                dcovmat(i1,i2,j) = dcovmat(i2,i1,j);
            }
        }
    }}

    return dcovmat;
}



//' Isotropic Matern covariance function on sphere
//'
//' From a matrix of longitudes and latitudes and a vector covariance parameters of the form
//' (variance, range, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param lonlat A matrix with \code{n} rows and one column with longitudes in (-180,180)
//' and one column of latitudes in (-90,90).
//' Each row of locs describes a point on the sphere.
//' @param covparms A vector giving positive-valued covariance parameters
//' in the form (variance, range, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlat[i,]} and
//' \code{lonlat[j,]}.
//' @section Matern on Sphere Domain:
//' The function first calculates the (x,y,z) 3D coordinates, and then inputs
//' the resulting locations into \code{maternIsotropic}. This means that we construct
//' covariances on the sphere by embedding the sphere in a 3D space. There has been some
//' concern expressed in the literature that such embeddings may produce distortions.
//' The source and nature of such distortions has never been articulated,
//' and to date, no such distortions have been documented. Guinness and
//' Fuentes (2016) argue that 3D embeddings produce reasonable models for data on spheres.
// [[Rcpp::export]]
arma::mat matern_sphere(NumericVector covparms, NumericMatrix lonlat ){

    int n = lonlat.nrow();
    double lonrad;                                  // longitude
    double latrad;                                  // latitude
    Rcpp::NumericMatrix xyz(n, 3);
    for(int i = 0; i < n; i++){
        lonrad = 2*M_PI*lonlat(i,0)/360;
        latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }

    arma::mat covmat = matern_isotropic( covparms, xyz );
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_sphere(NumericVector covparms, NumericMatrix lonlat ){

    int n = lonlat.nrow();
    Rcpp::NumericMatrix xyz(n, 3);
    for(int i = 0; i < n; i++){
        double lonrad = 2*M_PI*lonlat(i,0)/360;
        double latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }

    arma::cube dcovmat = d_matern_isotropic( covparms, xyz );
    return dcovmat;
}





// [[Rcpp::export]]
arma::mat matern_spheretime(NumericVector covparms, NumericMatrix lonlattime ){
    
    int n = lonlattime.nrow();
    // matrix to hold (x,y,z,t)
    NumericMatrix locs(n, 4);
    // convert lonlat to x,y,z
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattime(i,0)/360;
        double latrad = 2*M_PI*(lonlattime(i,1)+90)/360;
        locs(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        locs(i,1) = sin(latrad)*sin(lonrad);
        locs(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){
        locs(i,3) = lonlattime(i,2);
    }
    arma::mat covmat = matern_spacetime( covparms, locs );
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_spheretime(NumericVector covparms, NumericMatrix lonlattime){
    
    int n = lonlattime.nrow();
    // matrix to hold (x,y,z,t)
    NumericMatrix locs(n, 4);
    // convert lonlat to x,y,z
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattime(i,0)/360;
        double latrad = 2*M_PI*(lonlattime(i,1)+90)/360;
        locs(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        locs(i,1) = sin(latrad)*sin(lonrad);
        locs(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){
        locs(i,3) = lonlattime(i,2);
    }
    arma::cube dcovmat = d_matern_spacetime( covparms, locs );
    return dcovmat;
}




// [[Rcpp::export]]
arma::mat matern_sphere_warp(NumericVector covparms, NumericMatrix lonlat ){

    int n = lonlat.nrow();
    int nisoparms = 4;
    NumericVector isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.length() - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    Rcpp::NumericMatrix xyz(n, 3);
    for(int i = 0; i < n; i++){
        double lonrad = 2*M_PI*lonlat(i,0)/360;
        double latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }
    
    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyz, Lmax );

    // warp locations
    for(int i=0; i<n; i++){
        for(int k=0; k<3; k++){
            for(int j=0; j<nbasis; j++){
                xyz(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
            }
        }
    }
    
    // compute covariances
    arma::mat covmat = matern_isotropic( isoparms, xyz );
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_sphere_warp(NumericVector covparms, NumericMatrix lonlat ){

    int n = lonlat.nrow();
    int nisoparms = 4;
    NumericVector isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.length() - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    Rcpp::NumericMatrix xyz(n, 3);
    for(int i = 0; i < n; i++){
        double lonrad = 2*M_PI*lonlat(i,0)/360;
        double latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }
    
    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyz, Lmax );

    // warp locations
    for(int i=0; i<n; i++){
        for(int k=0; k<3; k++){
            for(int j=0; j<nbasis; j++){
                xyz(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
            }
        }
    }

    arma::cube ddcov = d_matern_isotropic( isoparms, xyz );
    arma::cube dcovmat(n,n,covparms.length());
    for(int i=0; i<nisoparms; i++){ dcovmat.slice(i) = ddcov.slice(i); }
    for(int j=0; j<nbasis; j++){
        for(int i1=0; i1<n; i1++){ for(int i2=i1; i2<n; i2++){
            // compute distance
            double dd = 0.0;
            for(int k=0; k<3; k++){ dd += pow( xyz(i2,k) - xyz(i1,k), 2.0); }
            dd = pow(dd,0.5);
            if(dd==0.0){ 
                dcovmat(i2,i1,nisoparms+j) = 0.0; 
            } else {
                // make use of already computed derivatives wrt covparms(1)
                dcovmat(i2,i1,nisoparms+j) = -dcovmat(i2,i1,1)*covparms(1)/dd;
                // pd wrt to location component * pd wrt to basis
                double pd = 0.0;
                for(int k=0; k<3; k++){ 
                    pd += (xyz(i2,k)-xyz(i1,k))/dd * grad_basis(i2,j,k);
                    pd += (xyz(i1,k)-xyz(i2,k))/dd * grad_basis(i1,j,k);
                }
                dcovmat(i2,i1,nisoparms+j) *= pd;
            }
            // fill in opposite side
            dcovmat(i1,i2,nisoparms+j) = dcovmat(i2,i1,nisoparms+j);
        }}
    }
    return dcovmat;
}







// [[Rcpp::export]]
arma::mat matern_spheretime_warp(NumericVector covparms, NumericMatrix lonlattime ){

    int n = lonlattime.nrow();
    int nisoparms = 5;
    NumericVector isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.length() - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    Rcpp::NumericMatrix xyzt(n, 4);
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattime(i,0)/360;
        double latrad = 2*M_PI*(lonlattime(i,1)+90)/360;
        xyzt(i,0) = sin(latrad)*cos(lonrad); // convert lon,lat to x,y,z
        xyzt(i,1) = sin(latrad)*sin(lonrad);
        xyzt(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){ xyzt(i,3) = lonlattime(i,2); }

    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyzt , Lmax );

    // warp locations
    for(int i=0; i<n; i++){ for(int k=0; k<3; k++){ for(int j=0; j<nbasis; j++){
        xyzt(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
    }}}

    // compute covariances
    arma::mat covmat = matern_spacetime( isoparms, xyzt );
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_spheretime_warp(NumericVector covparms, NumericMatrix lonlattime ){

    int n = lonlattime.nrow();
    int nisoparms = 5;
    NumericVector isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.length() - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    Rcpp::NumericMatrix xyzt(n, 4);
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattime(i,0)/360;
        double latrad = 2*M_PI*(lonlattime(i,1)+90)/360;
        xyzt(i,0) = sin(latrad)*cos(lonrad); // convert lon,lat to x,y,z
        xyzt(i,1) = sin(latrad)*sin(lonrad);
        xyzt(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){ xyzt(i,3) = lonlattime(i,2); }

    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyzt , Lmax );
    
    // warp locations
    for(int i=0; i<n; i++){ for(int k=0; k<3; k++){ for(int j=0; j<nbasis; j++){
        xyzt(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
    }}}
    
    arma::cube ddcov = d_matern_spacetime( isoparms, xyzt );
    arma::cube dcovmat(n,n,covparms.length());
    for(int i=0; i<nisoparms; i++){ dcovmat.slice(i) = ddcov.slice(i); }

    for(int j=0; j<nbasis; j++){
        for(int i1=0; i1<n; i1++){ for(int i2=i1; i2<n; i2++){
            // compute distance
            double dd = 0.0;
            for(int k=0; k<3; k++){ 
                dd += pow( (xyzt(i2,k) - xyzt(i1,k))/covparms(1), 2.0); 
            }
            dd += pow( (xyzt(i2,3)-xyzt(i1,3))/covparms(2), 2.0 );
            dd = pow(dd,0.5);
            if(dd==0.0){ 
                dcovmat(i2,i1,nisoparms+j) = 0.0; 
            } else {
                // make use of already computed derivatives wrt covparms(1)
                double pth = 0.0;
                // derivative of distance wrt to temporal range parameter
                pth += -(xyzt(i1,3)-xyzt(i2,3))/pow(covparms(2),3)/dd*xyzt(i1,3);
                pth += -(xyzt(i2,3)-xyzt(i1,3))/pow(covparms(2),3)/dd*xyzt(i2,3);
                // derivative of distance wrt to spatial range parameter
                for(int k=0; k<3; k++){ 
                    pth += -(xyzt(i1,k)-xyzt(i2,k))/pow(covparms(1),3)/dd * xyzt(i1,k);
                    pth += -(xyzt(i2,k)-xyzt(i1,k))/pow(covparms(1),3)/dd * xyzt(i2,k);
                }
                dcovmat(i2,i1,nisoparms+j) = (dcovmat(i2,i1,1)+dcovmat(i2,i1,2))/pth;
                // pd wrt to location component * pd wrt to basis
                double pd = 0.0;
                for(int k=0; k<3; k++){ 
                    pd += (xyzt(i2,k)-xyzt(i1,k))/pow(covparms(1),2)/dd * grad_basis(i2,j,k);
                    pd += (xyzt(i1,k)-xyzt(i2,k))/pow(covparms(1),2)/dd * grad_basis(i1,j,k);
                }
                dcovmat(i2,i1,nisoparms+j) *= pd;
            }
            // fill in opposite side
            dcovmat(i1,i2,nisoparms+j) = dcovmat(i2,i1,nisoparms+j);

        }}
    }
    return dcovmat;
}









// [[Rcpp::export]]
arma::mat matern_spheretime_nonstatvar(NumericVector covparms, NumericMatrix lonlattimeZ ){
    
    int n = lonlattimeZ.nrow();
    int dim = 3;
    int nbasis = lonlattimeZ.ncol() - dim;
    int nisoparm = 5;
    NumericVector isoparms(nisoparm);
    for(int i=0; i<nisoparm; i++){ isoparms(i) = covparms(i); }

    // matrix to hold (x,y,z,t)
    NumericMatrix locs(n, 4);
    // convert lonlat to x,y,z
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattimeZ(i,0)/360;
        double latrad = 2*M_PI*(lonlattimeZ(i,1)+90)/360;
        locs(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        locs(i,1) = sin(latrad)*sin(lonrad);
        locs(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){
        locs(i,3) = lonlattimeZ(i,2);
    }
    arma::mat covmat = matern_spacetime( isoparms, locs );
    
    arma::mat logvarmat(n,n);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<n; i2++){
        logvarmat(i2,i1) = 0.0;
        for(int k=0; k<nbasis; k++){
            logvarmat(i2,i1) += covparms(nisoparm+k)*
                ( lonlattimeZ(i1,dim+k) + lonlattimeZ(i2,dim+k) );
        }
    }}
    covmat = covmat % arma::exp(logvarmat);
    return covmat;
}

// [[Rcpp::export]]
arma::cube d_matern_spheretime_nonstatvar(NumericVector covparms, NumericMatrix lonlattimeZ){
    
    int n = lonlattimeZ.nrow();
    int dim = 3;
    int nbasis = lonlattimeZ.ncol() - dim;
    int nisoparm = 5;
    NumericVector isoparms(nisoparm);
    for(int i=0; i<nisoparm; i++){ isoparms(i) = covparms(i); }

    // matrix to hold (x,y,z,t)
    NumericMatrix locs(n, 4);
    // convert lonlat to x,y,z
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattimeZ(i,0)/360;
        double latrad = 2*M_PI*(lonlattimeZ(i,1)+90)/360;
        locs(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        locs(i,1) = sin(latrad)*sin(lonrad);
        locs(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){
        locs(i,3) = lonlattimeZ(i,2);
    }
    arma::mat covmat = matern_spacetime( isoparms, locs );
    
    arma::cube ddiso = d_matern_spacetime( isoparms, locs );
    arma::cube dcovmat(n,n,covparms.length());
    for(int k=0; k<nisoparm; k++){ dcovmat.slice(k) = ddiso.slice(k); }
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<n; i2++){
        for(int k=0; k<nbasis; k++){
            dcovmat(i2,i1,nisoparm+k) = 
                ( lonlattimeZ(i1,dim+k) + lonlattimeZ(i2,dim+k) )*covmat(i2,i1);
        }
    }}
    return dcovmat;
}



covfun_t get_covfun(std::string covfun_name_string)
{
    
    covfun_t covstruct;
    if( covfun_name_string.compare("matern_isotropic") == 0 )
    { 
        covstruct.p_covfun = &matern_isotropic; 
        covstruct.p_d_covfun = &d_matern_isotropic;
    } 
    else if( covfun_name_string.compare("matern_anisotropic2D") == 0 )
    { 
        covstruct.p_covfun = &matern_anisotropic2D; 
        covstruct.p_d_covfun = &d_matern_anisotropic2D;
    } 
    else if( covfun_name_string.compare("exponential_anisotropic3D") == 0 )
    { 
        covstruct.p_covfun = &exponential_anisotropic3D; 
        covstruct.p_d_covfun = &d_exponential_anisotropic3D;
    } 
    else if( covfun_name_string.compare("matern_anisotropic3D") == 0 )
    { 
        covstruct.p_covfun = &matern_anisotropic3D; 
        covstruct.p_d_covfun = &d_matern_anisotropic3D;
    } 
    else if( covfun_name_string.compare("matern_nonstat_var") == 0 )
    { 
        covstruct.p_covfun = &matern_nonstat_var; 
        covstruct.p_d_covfun = &d_matern_nonstat_var;
    } 
    else if( covfun_name_string.compare("exponential_isotropic") == 0 )
    { 
        covstruct.p_covfun = &exponential_isotropic; 
        covstruct.p_d_covfun = &d_exponential_isotropic;
    }
    else if( covfun_name_string.compare("matern_sphere") == 0 )
    { 
        covstruct.p_covfun = &matern_sphere; 
        covstruct.p_d_covfun = &d_matern_sphere;
    }
    else if( covfun_name_string.compare("matern_sphere_warp") == 0 )
    { 
        covstruct.p_covfun = &matern_sphere_warp; 
        covstruct.p_d_covfun = &d_matern_sphere_warp;
    }
    else if( covfun_name_string.compare("matern_spheretime_warp") == 0 )
    { 
        covstruct.p_covfun = &matern_spheretime_warp; 
        covstruct.p_d_covfun = &d_matern_spheretime_warp;
    }
    else if( covfun_name_string.compare("matern_spheretime") == 0 )
    { 
        covstruct.p_covfun = &matern_spheretime; 
        covstruct.p_d_covfun = &d_matern_spheretime;
    }
    else if( covfun_name_string.compare("matern_spheretime_nonstatvar") == 0 )
    { 
        covstruct.p_covfun = &matern_spheretime_nonstatvar; 
        covstruct.p_d_covfun = &d_matern_spheretime_nonstatvar;
    }
    else if( covfun_name_string.compare("matern_spacetime") == 0 )
    { 
        covstruct.p_covfun = &matern_spacetime; 
        covstruct.p_d_covfun = &d_matern_spacetime;
    }
    else if( covfun_name_string.compare("matern_scaledim") == 0 )
    { 
        covstruct.p_covfun = &matern_scaledim;
        covstruct.p_d_covfun = &d_matern_scaledim;
    }
    else { // stop the program
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
        assert(0);
    }
    return covstruct;
}


#endif
