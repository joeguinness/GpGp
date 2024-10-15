#ifndef COVMATRIX_FUNS_02_H
#define COVMATRIX_FUNS_02_H

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


//' Geometrically anisotropic Matern covariance function (two dimensions)
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, L11, L21, L22, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{2} columns.
//' Each row of locs is a point in R^2.
//' @param covparms A vector with covariance parameters
//' in the form (variance, L11, L21, L22, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, L11, L21, L22, smoothness, nugget)
//' where L11, L21, L22, are the three non-zero entries of a lower-triangular
//' matrix L. The covariances are 
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) (|| L x - L y || )^\nu K_\nu(|| L x - L y ||) }
//' This means that L11 is interpreted as an inverse range parameter in the
//' first dimension.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern_anisotropic2D(arma::vec covparms, arma::mat locs ){
    
	// fail-safe to prevent large smoothness values
    covparms(4) = std::min( covparms(4), 8.0 );
	
    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = smoothness
    // covparms(5) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget
    
    //int dim = locs.n_cols;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 5 );
    double normcon = covparms(0)/(pow(2.0,covparms(4)-1.0)*boost::math::tgamma(covparms(4) ));
    
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
                    pow( d, covparms(4) )*boost::math::cyl_bessel_k(covparms(4), d);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

//' @describeIn matern_anisotropic2D Derivatives of anisotropic Matern covariance
// [[Rcpp::export]]
arma::cube d_matern_anisotropic2D(arma::vec covparms, arma::mat locs ){

	// fail-safe to prevent large smoothness values
    covparms(4) = std::min( covparms(4), 8.0 );
	
    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = smoothness
    // covparms(5) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget

    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 5 );
    double normcon = covparms(0)/(pow(2.0,covparms(4)-1.0)*boost::math::tgamma(covparms(4) ));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,covparms(4)+eps-1.0)*boost::math::tgamma(covparms(4)+eps ));
    
    double b0 = covparms(1)*covparms(1);
    double b1 = covparms(2)*covparms(2)+covparms(3)*covparms(3);
    double b2 = covparms(2)*covparms(3);

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
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
            cov = normcon*pow( d, covparms(4) )*boost::math::cyl_bessel_k(covparms(4), d);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // cholesky parameters
            double cov_nu_m1 = normcon*pow(d,covparms(4)-1.0)*
                boost::math::cyl_bessel_k(covparms(4)-1.0, d);  
            dcovmat(i1,i2,1) -= cov_nu_m1*(h0*h0*covparms(1));
            dcovmat(i1,i2,2) -= cov_nu_m1*(h1*h1*covparms(2) + h0*h1*covparms(3));
            dcovmat(i1,i2,3) -= cov_nu_m1*(h1*h1*covparms(3) + h0*h1*covparms(2));
            
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,4) += 
                ( normconeps*pow(d,covparms(4)+eps)*boost::math::cyl_bessel_k(covparms(4) + eps, d) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(5);
            dcovmat(i1,i2,5) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}



//' Geometrically anisotropic Matern covariance function (three dimensions)
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, L11, L21, L22, L31, L32, L33, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{3} columns.
//' Each row of locs is a point in R^3.
//' @param covparms A vector with covariance parameters
//' in the form (variance, L11, L21, L22, L31, L32, L33, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, L11, L21, L22, L31, L32, L33, smoothness, nugget)
//' where L11, L21, L22, L31, L32, L33 are the six non-zero entries of a lower-triangular
//' matrix L. The covariances are 
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) (|| L x - L y || )^\nu K_\nu(|| L x - L y ||) }
//' This means that L11 is interpreted as an inverse range parameter in the
//' first dimension.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern_anisotropic3D(arma::vec covparms, arma::mat locs ){
    
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
    
    //int dim = locs.n_cols;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 8 );
    double smooth = covparms( 7 );
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*boost::math::tgamma(smooth));
    
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
                    pow( d, smooth )*boost::math::cyl_bessel_k(smooth, d);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

//' @describeIn matern_anisotropic3D Derivatives of anisotropic Matern covariance
// [[Rcpp::export]]
arma::cube d_matern_anisotropic3D(arma::vec covparms, arma::mat locs ){

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

    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 8 );
    double smooth = covparms( 7 );
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*boost::math::tgamma(smooth));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,smooth+eps-1.0)*boost::math::tgamma(smooth+eps));
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
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
            cov = normcon*pow( d, smooth )*boost::math::cyl_bessel_k(smooth, d);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // cholesky parameters
            double cov_nu_m1 = normcon*pow( d, smooth - 1.0 )*
                boost::math::cyl_bessel_k(smooth - 1.0, d);  
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
                ( normconeps*pow(d,smooth+eps)*boost::math::cyl_bessel_k(smooth + eps, d) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(8);
            dcovmat(i1,i2,8) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}





//' Geometrically anisotropic Matern covariance function (three dimensions, alternate parameterization)
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, B11, B12, B13, B22, B23, B33, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{3} columns.
//' Each row of locs is a point in R^3.
//' @param covparms A vector with covariance parameters
//' in the form (variance, B11, B12, B13, B22, B23, B33, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, B11, B12, B13, B22, B23, B33, smoothness, nugget)
//' where B11, B12, B13, B22, B23, B33, transform the three coordinates as
//' \deqn{ u_1 = B11[ x_1 + B12 x_2 + (B13 + B12 B23) x_3] }
//' \deqn{ u_2 = B22[ x_2 + B23 x_3] }
//' \deqn{ u_3 = B33[ x_3 ] }
//' NOTE: the u_1 transformation is different from previous versions of this function.
//' NOTE: now (B13,B23) can be interpreted as a drift vector in space over time.
//' Assuming x is transformed to u and y transformed to v, the covariances are 
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) (|| u - v || )^\nu K_\nu(|| u - v ||) }
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern_anisotropic3D_alt(arma::vec covparms, arma::mat locs ){
    
    // covparms(0) = sigmasq
    // covparms(1) = B00
    // covparms(2) = B10
    // covparms(3) = B11
    // covparms(4) = B20
    // covparms(5) = B21
    // covparms(6) = B22
    // covparms(7) = smoothness
    // covparms(8) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget
    
    // fail-safe to prevent large smoothness values
    covparms(7) = std::min( covparms(7), 8.0 );
	
    //int dim = locs.n_cols;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 8 );
    double smooth = covparms( 7 );
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*boost::math::tgamma(smooth));
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate rescaled distance
            double h0 = locs(i1,0) - locs(i2,0);
            double h1 = locs(i1,1) - locs(i2,1);
            double h2 = locs(i1,2) - locs(i2,2);
	    double r0 = covparms(1)*( h0 + covparms(2)*h1 + (covparms(3)+covparms(2)*covparms(5))*h2 );
	    double r1 = covparms(4)*( h1 + covparms(5)*h2 );
	    double r2 = covparms(6)*( h2 );
	    
            // 3 is hard coded here
            double d = 0.0;
            d += pow( r0, 2 );
            d += pow( r1, 2 );
            d += pow( r2, 2 );
            d =  pow( d, 0.5 );
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0);
            } else {
                // calculate covariance            
                covmat(i2,i1) = normcon*
                    pow( d, smooth )*boost::math::cyl_bessel_k(smooth, d);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

//' @describeIn matern_anisotropic3D Derivatives of anisotropic Matern covariance
// [[Rcpp::export]]
arma::cube d_matern_anisotropic3D_alt(arma::vec covparms, arma::mat locs ){

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

	// fail-safe to prevent large smoothness values
	covparms(7) = std::min( covparms(7), 8.0 );
	
    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 8 );
    double smooth = covparms( 7 );
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*boost::math::tgamma(smooth));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,smooth+eps-1.0)*boost::math::tgamma(smooth+eps));
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        // calculate rescaled distance
        double h0 = locs(i1,0) - locs(i2,0);
        double h1 = locs(i1,1) - locs(i2,1);
        double h2 = locs(i1,2) - locs(i2,2);
        double r0 = covparms(1)*( h0 + covparms(2)*h1 + (covparms(3)+covparms(2)*covparms(5))*h2 );
        double r1 = covparms(4)*( h1 + covparms(5)*h2 );
        double r2 = covparms(6)*( h2 );
     
        // 3 is hard coded here
        double d = 0.0;
        d += pow( r0, 2 );
        d += pow( r1, 2 );
        d += pow( r2, 2 );
        d =  pow( d, 0.5 );
            
        double cov;        
        if( d == 0.0 ){

            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;

        } else {

	    // derivatives of r0, r1, r2 with respect to covparms
	    double dr0_dc1 = r0/covparms(1);
	    double dr0_dc2 = covparms(1)*(h1+covparms(5)*h2);
	    double dr0_dc3 = covparms(1)*h2;
	    double dr0_dc5 = covparms(1)*covparms(2)*h2;

	    double dr1_dc4 = r1/covparms(4);
	    double dr1_dc5 = covparms(4)*h2;

	    double dr2_dc6 = h2;

	    // derivative of scaled distance with respect to covparms
	    double dd_dc1 = r0*dr0_dc1/d;
	    double dd_dc2 = r0*dr0_dc2/d;
	    double dd_dc3 = r0*dr0_dc3/d;
	    double dd_dc4 = r1*dr1_dc4/d;
	    double dd_dc5 = r0*dr0_dc5/d + r1*dr1_dc5/d;
	    double dd_dc6 = r2*dr2_dc6/d;
	    
	    // derivatives of (d)^nu K_nu(d) with respect to covparms
	    // deriv of (d)^n K_nu(d) with respect to d is negative!
            double cov_nu_m1 = normcon*pow( d, smooth - 1.0 )*
                boost::math::cyl_bessel_k(smooth - 1.0, d);  
	    dcovmat(i1,i2,1) = -d*cov_nu_m1*dd_dc1; 
	    dcovmat(i1,i2,2) = -d*cov_nu_m1*dd_dc2; 
	    dcovmat(i1,i2,3) = -d*cov_nu_m1*dd_dc3; 
	    dcovmat(i1,i2,4) = -d*cov_nu_m1*dd_dc4; 
	    dcovmat(i1,i2,5) = -d*cov_nu_m1*dd_dc5; 
	    dcovmat(i1,i2,6) = -d*cov_nu_m1*dd_dc6; 

            // variance parameter
            cov = normcon*pow( d, smooth )*boost::math::cyl_bessel_k(smooth, d);  
            dcovmat(i1,i2,0) += cov/covparms(0);
    
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,7) += 
                ( normconeps*pow(d,smooth+eps)*boost::math::cyl_bessel_k(smooth + eps, d) -
                  cov )/eps;
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(8);
            dcovmat(i1,i2,8) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}





//' Geometrically anisotropic exponential covariance function (two dimensions)
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, L11, L21, L22, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{2} columns.
//' Each row of locs is a point in R^2.
//' @param covparms A vector with covariance parameters
//' in the form (variance, L11, L21, L22, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, L11, L21, L22, nugget)
//' where L11, L21, L22, are the three non-zero entries of a lower-triangular
//' matrix L. The covariances are 
//' \deqn{ M(x,y) = \sigma^2 exp(-|| L x - L y || ) }
//' This means that L11 is interpreted as an inverse range parameter in the
//' first dimension.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat exponential_anisotropic2D(arma::vec covparms, arma::mat locs ){
    
    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget
    
    //int dim = locs.n_cols;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 4 );

    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate rescaled distance
            double h0 = locs(i1,0) - locs(i2,0);
            double h1 = locs(i1,1) - locs(i2,1);
            // 2 dims is hard coded here
            double d = 0.0;
            d += pow( covparms(1)*h0, 2 );
            d += pow( covparms(2)*h0 + covparms(3)*h1, 2 );
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

//' @describeIn exponential_anisotropic2D Derivatives of anisotropic exponential covariance
// [[Rcpp::export]]
arma::cube d_exponential_anisotropic2D(arma::vec covparms, arma::mat locs ){

    // covparms(0) = sigmasq
    // covparms(1) = L00
    // covparms(2) = L10
    // covparms(3) = L11
    // covparms(4) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget

    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 4 );

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        // calculate rescaled distance
        double h0 = locs(i1,0) - locs(i2,0);
        double h1 = locs(i1,1) - locs(i2,1);
        // 2 spatial + 1 time dimen is hard coded here
        double d = 0.0;
        d += pow( covparms(1)*h0, 2 );
        d += pow( covparms(2)*h0 + covparms(3)*h1, 2 );
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









//' Geometrically anisotropic exponential covariance function (three dimensions)
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, L11, L21, L22, L31, L32, L33, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{3} columns.
//' Each row of locs is a point in R^3.
//' @param covparms A vector with covariance parameters
//' in the form (variance, L11, L21, L22, L31, L32, L33, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, L11, L21, L22, L31, L32, L33, nugget)
//' where L11, L21, L22, L31, L32, L33 are the six non-zero entries of a lower-triangular
//' matrix L. The covariances are 
//' \deqn{ M(x,y) = \sigma^2 exp(-|| L x - L y || ) }
//' This means that L11 is interpreted as an inverse range parameter in the
//' first dimension.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat exponential_anisotropic3D(arma::vec covparms, arma::mat locs ){
    
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
    
    //int dim = locs.n_cols;
    int n = locs.n_rows;
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

//' @describeIn exponential_anisotropic3D Derivatives of anisotropic exponential covariance
// [[Rcpp::export]]
arma::cube d_exponential_anisotropic3D(arma::vec covparms, arma::mat locs ){

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

    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 7 );

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
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
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}










//' Geometrically anisotropic exponential covariance function (three dimensions, alternate parameterization)
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, B11, B12, B13, B22, B23, B33, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{3} columns.
//' Each row of locs is a point in R^3.
//' @param covparms A vector with covariance parameters
//' in the form (variance, B11, B12, B13, B22, B23, B33, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, B11, B12, B13, B22, B23, B33, smoothness, nugget)
//' where B11, B12, B13, B22, B23, B33, transform the three coordinates as
//' \deqn{ u_1 = B11[ x_1 + B12 x_2 + (B13 + B12 B23) x_3] }
//' \deqn{ u_2 = B22[ x_2 + B23 x_3] }
//' \deqn{ u_3 = B33[ x_3 ] }
//' (B13,B23) can be interpreted as a drift vector in space over time
//' if first two dimensions are space and third is time.
//' Assuming x is transformed to u and y transformed to v, the covariances are 
//' \deqn{ M(x,y) = \sigma^2 exp( - || u - v || ) }
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat exponential_anisotropic3D_alt(arma::vec covparms, arma::mat locs ){
    
    // covparms(0) = sigmasq
    // covparms(1) = B00
    // covparms(2) = B10
    // covparms(3) = B11
    // covparms(4) = B20
    // covparms(5) = B21
    // covparms(6) = B22
    // covparms(7) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget
    
    //int dim = locs.n_cols;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 7 );
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate rescaled distance
            double h0 = locs(i1,0) - locs(i2,0);
            double h1 = locs(i1,1) - locs(i2,1);
            double h2 = locs(i1,2) - locs(i2,2);
	    double r0 = covparms(1)*( h0 + covparms(2)*h1 + (covparms(3)+covparms(2)*covparms(5))*h2 );
	    double r1 = covparms(4)*( h1 + covparms(5)*h2 );
	    double r2 = covparms(6)*( h2 );
	    
            // 3 is hard coded here
            double d = 0.0;
            d += pow( r0, 2 );
            d += pow( r1, 2 );
            d += pow( r2, 2 );
            d =  pow( d, 0.5 );
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0);
            } else {
                // calculate covariance            
                covmat(i2,i1) = covparms(0)*std::exp(-d);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}

//' @describeIn exponential_anisotropic3D_alt Derivatives of anisotropic Matern covariance
// [[Rcpp::export]]
arma::cube d_exponential_anisotropic3D_alt(arma::vec covparms, arma::mat locs ){

    // covparms(0) = sigmasq
    // covparms(1) = B00
    // covparms(2) = B01
    // covparms(3) = B02
    // covparms(4) = B11
    // covparms(5) = B12
    // covparms(6) = B22
    // covparms(7) = tausq
    // nugget = sigmasq*tausq
    // overall variance = sigmasq*(1 + tausq) = sigmasq + nugget

    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 7 );
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
        
        // calculate rescaled distance
        double h0 = locs(i1,0) - locs(i2,0);
        double h1 = locs(i1,1) - locs(i2,1);
        double h2 = locs(i1,2) - locs(i2,2);
        double r0 = covparms(1)*( h0 + covparms(2)*h1 + (covparms(3)+covparms(2)*covparms(5))*h2 );
        double r1 = covparms(4)*( h1 + covparms(5)*h2 );
        double r2 = covparms(6)*( h2 );
     
        // 3 is hard coded here
        double d = 0.0;
        d += pow( r0, 2 );
        d += pow( r1, 2 );
        d += pow( r2, 2 );
        d =  pow( d, 0.5 );
            
        double cov;        
        if( d == 0.0 ){

            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;

        } else {

	    cov = covparms(0)*std::exp(-d);

	    // derivatives of r0, r1, r2 with respect to covparms
	    double dr0_dc1 = r0/covparms(1);
	    double dr0_dc2 = covparms(1)*(h1+covparms(5)*h2);
	    double dr0_dc3 = covparms(1)*h2;
	    double dr0_dc5 = covparms(1)*covparms(2)*h2;

	    double dr1_dc4 = r1/covparms(4);
	    double dr1_dc5 = covparms(4)*h2;

	    double dr2_dc6 = h2;

	    // derivative of scaled distance with respect to covparms
	    double dd_dc1 = r0*dr0_dc1/d;
	    double dd_dc2 = r0*dr0_dc2/d;
	    double dd_dc3 = r0*dr0_dc3/d;
	    double dd_dc4 = r1*dr1_dc4/d;
	    double dd_dc5 = r0*dr0_dc5/d + r1*dr1_dc5/d;
	    double dd_dc6 = r2*dr2_dc6/d;
	    
	    // derivatives of covparms(0)*exp(-d)
	    dcovmat(i1,i2,1) = -cov*dd_dc1; 
	    dcovmat(i1,i2,2) = -cov*dd_dc2; 
	    dcovmat(i1,i2,3) = -cov*dd_dc3; 
	    dcovmat(i1,i2,4) = -cov*dd_dc4; 
	    dcovmat(i1,i2,5) = -cov*dd_dc5; 
	    dcovmat(i1,i2,6) = -cov*dd_dc6; 

            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
    
        }
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(7);
            dcovmat(i1,i2,7) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}







#endif
