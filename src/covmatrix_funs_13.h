#ifndef COVMATRIX_FUNS_13_H
#define COVMATRIX_FUNS_13_H



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


//' Isotropic Matern covariance function with random effects for categories
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range, smoothness, category variance, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs gives a point in R^d.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range, smoothness, category variance, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range, smoothness, category variance, nugget)
//' = \eqn{(\sigma^2,\alpha,\nu,c^2,\tau^2)}, and the covariance function is parameterized
//' as
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) (|| x - y ||/\alpha )^\nu K_\nu(|| x - y ||/\alpha ) }
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' The category variance \eqn{c^2} is added if two observation from same category
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern_categorical(arma::vec covparms, arma::mat locs ){

    // locs has dim+1 columns.
    // The first dim columns are coordinates
    // the last column is the category
    
    // parameters are
    // covparms(0) : sill variance
    // covparms(1) : range
    // covparms(2) : smoothness
    // covparms(3) : variance of categorical component
    // covparms(4) : nugget
    
    // fail-safe to prevent large smoothness values
    covparms(2) = std::min( covparms(2), 8.0 );
	
    int dim = locs.n_cols - 1;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 4 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*boost::math::tgamma(covparms(2) ));
	
    // create scaled locations
    mat locs_scaled = locs;
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
                    pow( d, covparms(2) )*boost::math::cyl_bessel_k(covparms(2), d);
            }

            // if same category, add categorical component
            // should we do abs( a - b ) < 1e-6 instead of a == b?
            if( locs_scaled(i1,dim) == locs_scaled(i2,dim) ){
                covmat(i2,i1) += covparms(3);
            }
			
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}



//' @describeIn matern_categorical Derivatives of isotropic Matern covariance
// [[Rcpp::export]]
arma::cube d_matern_categorical(arma::vec covparms, arma::mat locs ){

    // locs has dim+1 columns.
    // The first dim columns are coordinates
    // the last column is the category
    
    // parameters are
    // covparms(0) : sill variance
    // covparms(1) : range
    // covparms(2) : smoothness
    // covparms(3) : variance of categorical component
    // covparms(4) : nugget
    
    // fail-safe to prevent large smoothness values
    covparms(2) = std::min( covparms(2), 8.0 );
	
    int dim = locs.n_cols - 1;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 4 );
	double categ = covparms(0)*covparms(3);
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*boost::math::tgamma(covparms(2) ));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,covparms(2)+eps-1.0)*boost::math::tgamma(covparms(2)+eps ));
    
    // create scaled locations
    mat locs_scaled = locs;
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
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
            cov = normcon*pow( d, covparms(2) )*boost::math::cyl_bessel_k(covparms(2), d);
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);
            // range parameter
            dcovmat(i1,i2,1) += normcon*pow(d,covparms(2))*
                boost::math::cyl_bessel_k(covparms(2)-1.0, d)*d/covparms(1);
            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,2) += 
                ( normconeps*pow(d,covparms(2)+eps)*boost::math::cyl_bessel_k(covparms(2) + eps, d) -
                  cov )/eps;
        }
		// update if same category
		if( locs_scaled(i1,dim) == locs_scaled(i2,dim) ){
			dcovmat(i1,i2,3) += 1.0;
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



//' Space-Time Matern covariance function with random effects for categories
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, spatial range, temporal range, smoothness, category, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs gives a point in R^d.
//' @param covparms A vector with covariance parameters
//' in the form (variance, spatial range, temporal range, smoothness, category, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range, smoothness, category, nugget)
//' = \eqn{(\sigma^2,\alpha_1,\alpha_2,\nu,c^2,\tau^2)}, and the covariance function is parameterized
//' as
//' \deqn{ d = ( || x - y ||^2/\alpha_1 + |s-t|^2/\alpha_2^2 )^{1/2} }
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) (d)^\nu K_\nu(d) }
//' (x,s) and (y,t) are the space-time locations of a pair of observations.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' The category variance \eqn{c^2} is added if two observation from same category
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern_spacetime_categorical(arma::vec covparms, arma::mat locs ){

    // locs has dim+2 columns.
    // The first dim columns are euclidean spatial coordinates
    // next column is time dimension
    // the last column is the category
    
    // parameters are
    // covparms(0) : sill variance
    // covparms(1) : spatial range
    // covparms(2) : temporal range
    // covparms(3) : smoothness
    // covparms(4) : variance of categorical component
    // covparms(5) : nugget
    
    // fail-safe to prevent large smoothness values
    covparms(3) = std::min( covparms(3), 8.0 );
	
    int dim = locs.n_cols - 2;
    int n = locs.n_rows;
    double smooth = covparms(3);
    double nugget = covparms(0)*covparms(5);
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*boost::math::tgamma(smooth ));
	
    // create scaled locations
    mat locs_scaled = locs;
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    for(int i=0; i<n; i++){
        locs_scaled(i,dim) = locs(i,dim)/covparms(2);
    }
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate distance
            double d = 0.0;
            for(int j=0; j<dim+1; j++){
                d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
            }
            d = pow( d, 0.5 );
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0);
            } else {
                // calculate covariance        
                covmat(i2,i1) = normcon*
                    pow( d, smooth )*boost::math::cyl_bessel_k(smooth, d);
            }

            // if same category, add categorical component
            // should we do abs( a - b ) < 1e-6 instead of a == b?
            if( locs_scaled(i1,dim+1) == locs_scaled(i2,dim+1) ){
                covmat(i2,i1) += covparms(4);
            }
			
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}



//' @describeIn matern_spacetime_categorical Derivatives of isotropic Matern covariance
// [[Rcpp::export]]
arma::cube d_matern_spacetime_categorical(arma::vec covparms, arma::mat locs ){

    // locs has dim+2 columns.
    // The first dim columns are euclidean spatial coordinates
    // next column is time dimension
    // the last column is the category
    
    // parameters are
    // covparms(0) : sill variance
    // covparms(1) : spatial range
    // covparms(2) : temporal range
    // covparms(3) : smoothness
    // covparms(4) : variance of categorical component
    // covparms(5) : nugget
    
    // fail-safe to prevent large smoothness values
    covparms(3) = std::min( covparms(3), 8.0 );
	
    int dim = locs.n_cols - 2;
    int n = locs.n_rows;
    double smooth = covparms(3);
    double normcon = covparms(0)/(pow(2.0,smooth-1.0)*boost::math::tgamma(smooth));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,smooth + eps-1.0)*boost::math::tgamma(smooth + eps ));
    
    // create scaled locations
    mat locs_scaled = locs;
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
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){

        double d = 0.0;
        double d1 = 0.0;
        double d2 = 0.0;
        for(int j=0; j<dim; j++){
            d1 += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
	d2 = pow( locs_scaled(i1,dim) - locs_scaled(i2,dim), 2.0 );
        d = pow( d1 + d2, 0.5 );
        
        double cov;        
        if( d == 0.0 ){
            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;
            dcovmat(i1,i2,1) += 0.0;
            dcovmat(i1,i2,2) += 0.0;
        } else {

            cov = normcon*pow( d, smooth )*boost::math::cyl_bessel_k( smooth, d );
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);

	    // derivatives of (d)^nu K_nu(d) with respect to covparms
	    // deriv of (d)^nu K_nu(d) with respect to d is negative!
	    // note that normcon contains covparms(0) 
            double cov_nu_m1 = normcon*pow( d, smooth - 1.0 )*
                boost::math::cyl_bessel_k(smooth - 1.0, d);  
	    double dM_dd = -d*cov_nu_m1;
	    
            double dd_dc1 = -1/d*d1/covparms(1);
            double dd_dc2 = -1/d*d2/covparms(2);

	    dcovmat(i1,i2,1) = dM_dd * dd_dc1;
	    dcovmat(i1,i2,2) = dM_dd * dd_dc2;

            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,3) += 
                ( normconeps*pow(d,smooth+eps)*boost::math::cyl_bessel_k(smooth + eps, d) -
                  cov )/eps;
        }
        // update if same category
        if( locs_scaled(i1,dim+1) == locs_scaled(i2,dim+1) ){
            dcovmat(i1,i2,4) += 1.0;
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







//' Space-Time Matern covariance function with local random effects for categories
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, spatial range, temporal range, smoothness, cat variance, cat spatial range, cat temporal range, cat smoothness, nugget),
//' return the square matrix of
//' all pairwise covariances.
//' This is the covariance for the following model for data from cateogory k
//' \deqn{ Y_k(x_i,t_i) = Z_0(x_i,t_i) + Z_k(x_i,t_i) + e_i }
//' where Z_0 is Matern with parameters (variance,spatial range,temporal range,smoothness)
//' and Z_1,...,Z_K are independent Materns with parameters
//' (cat variance, cat spatial range, cat temporal range, cat smoothness),
//' and e_1, ..., e_n are independent normals with variance (variance * nugget)
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs gives a point in R^d.
//' @param covparms A vector with covariance parameters
//' in the form (variance, spatial range, temporal range, smoothness, category, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range, smoothness, category, nugget)
//' = \eqn{(\sigma^2,\alpha_1,\alpha_2,\nu,c^2,\tau^2)}, and the covariance function is parameterized
//' as
//' \deqn{ d = ( || x - y ||^2/\alpha_1 + |s-t|^2/\alpha_2^2 )^{1/2} }
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) (d)^\nu K_\nu(d) }
//' (x,s) and (y,t) are the space-time locations of a pair of observations.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' The category variance \eqn{c^2} is added if two observation from same category
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern_spacetime_categorical_local(arma::vec covparms, arma::mat locs ){

    // locs has dim+2 columns.
    // The first dim columns are euclidean spatial coordinates
    // next column is time dimension
    // the last column is the category
    
    // parameters are
    // covparms(0) : sill variance
    // covparms(1) : spatial range
    // covparms(2) : temporal range
    // covparms(3) : smoothness
    // covparms(4) : variance of categorical component
    // covparms(5) : spatial range of categorical component
    // covparms(6) : temporal range of categorical component
    // covparms(7) : smoothness of categorical component
    // covparms(8) : nugget
    
    // fail-safe to prevent large smoothness values
    covparms(3) = std::min( covparms(3), 8.0 );
    covparms(7) = std::min( covparms(7), 8.0 );
	
    int dim = locs.n_cols - 2;
    int n = locs.n_rows;
    double smooth1 = covparms(3);
    double smooth2 = covparms(7);
    double nugget = covparms(0)*covparms(8);
    double normcon1 = covparms(0)/(pow(2.0,smooth1-1.0)*boost::math::tgamma(smooth1 ));
    double normcon2 = covparms(4)/(pow(2.0,smooth2-1.0)*boost::math::tgamma(smooth2 ));
	
    // create scaled locations
    mat locs_scaled1 = locs;
    mat locs_scaled2 = locs;
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled1(i,j) = locs(i,j)/covparms(1);
            locs_scaled2(i,j) = locs(i,j)/covparms(5);
        }
    }
    for(int i=0; i<n; i++){
        locs_scaled1(i,dim) = locs(i,dim)/covparms(2);
        locs_scaled2(i,dim) = locs(i,dim)/covparms(6);
    }
    
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            
            // calculate distance
            double d = 0.0;
            for(int j=0; j<dim+1; j++){
                d += pow( locs_scaled1(i1,j) - locs_scaled1(i2,j), 2.0 );
            }
            d = pow( d, 0.5 );
            
            if( d == 0.0 ){
                covmat(i2,i1) = covparms(0);
            } else {
                // calculate covariance        
                covmat(i2,i1) = normcon1*
                    pow( d, smooth1 )*boost::math::cyl_bessel_k(smooth1, d);
            }

            // if same category, add categorical component
            // should we do abs( a - b ) < 1e-6 instead of a == b?
            if( locs(i1,dim+1) == locs(i2,dim+1) ){

                // calculate distance
                double d = 0.0;
                for(int j=0; j<dim+1; j++){
                    d += pow( locs_scaled2(i1,j) - locs_scaled2(i2,j), 2.0 );
                }
                d = pow( d, 0.5 );
                
                if( d == 0.0 ){
                    covmat(i2,i1) += covparms(4);
                } else {
                    // calculate covariance        
                    covmat(i2,i1) += normcon2*
                        pow( d, smooth2 )*boost::math::cyl_bessel_k(smooth2, d);
                }

            }
			
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}



//' @describeIn matern_spacetime_categorical_local Derivatives of isotropic Matern covariance
// [[Rcpp::export]]
arma::cube d_matern_spacetime_categorical_local(arma::vec covparms, arma::mat locs ){

    // locs has dim+2 columns.
    // The first dim columns are euclidean spatial coordinates
    // next column is time dimension
    // the last column is the category
    
    // parameters are
    // covparms(0) : sill variance
    // covparms(1) : spatial range
    // covparms(2) : temporal range
    // covparms(3) : smoothness
    // covparms(4) : variance of categorical component
    // covparms(5) : spatial range of categorical component
    // covparms(6) : temporal range of categorical component
    // covparms(7) : smoothness of categorical component
    // covparms(8) : nugget
    
    // fail-safe to prevent large smoothness values
    covparms(3) = std::min( covparms(3), 8.0 );
    covparms(7) = std::min( covparms(7), 8.0 );
	
    int dim = locs.n_cols - 2;
    int n = locs.n_rows;
    double smooth1 = covparms(3);
    double smooth2 = covparms(7);
    double nugget = covparms(0)*covparms(8);
    double normcon1 = covparms(0)/(pow(2.0,smooth1-1.0)*boost::math::tgamma(smooth1 ));
    double normcon2 = covparms(4)/(pow(2.0,smooth2-1.0)*boost::math::tgamma(smooth2 ));
	
    double eps = 1e-8;
    double normconeps1 = 
        covparms(0)/(pow(2.0,smooth1 + eps-1.0)*boost::math::tgamma(smooth1 + eps ));
    double normconeps2 = 
        covparms(4)/(pow(2.0,smooth2 + eps-1.0)*boost::math::tgamma(smooth2 + eps ));
    
    // create scaled locations
    mat locs_scaled1 = locs;
    mat locs_scaled2 = locs;
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled1(i,j) = locs(i,j)/covparms(1);
            locs_scaled2(i,j) = locs(i,j)/covparms(5);
        }
    }
    for(int i=0; i<n; i++){
        locs_scaled1(i,dim) = locs(i,dim)/covparms(2);
        locs_scaled2(i,dim) = locs(i,dim)/covparms(6);
    }
    
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){

        double d = 0.0;
        double d1 = 0.0;
        double d2 = 0.0;
        for(int j=0; j<dim; j++){
            d1 += pow( locs_scaled1(i1,j) - locs_scaled1(i2,j), 2.0 );
        }
	d2 = pow( locs_scaled1(i1,dim) - locs_scaled1(i2,dim), 2.0 );
        d = pow( d1 + d2, 0.5 );
        
        double cov;        
        if( d == 0.0 ){

            cov = covparms(0);
            dcovmat(i1,i2,0) += 1.0;

        } else {

            cov = normcon1*pow( d, smooth1 )*boost::math::cyl_bessel_k( smooth1, d );
            // variance parameter
            dcovmat(i1,i2,0) += cov/covparms(0);

	    // derivatives of (d)^nu K_nu(d) with respect to covparms
	    // deriv of (d)^nu K_nu(d) with respect to d is negative!
	    // note that normcon contains covparms(0) 
            double cov_nu_m1 = normcon1*pow( d, smooth1 - 1.0 )*
                boost::math::cyl_bessel_k(smooth1 - 1.0, d);  
	    double dM_dd = -d*cov_nu_m1;
	    
            double dd_dc1 = -1/d*d1/covparms(1);
            double dd_dc2 = -1/d*d2/covparms(2);

	    dcovmat(i1,i2,1) = dM_dd * dd_dc1;
	    dcovmat(i1,i2,2) = dM_dd * dd_dc2;

            // smoothness parameter (finite differencing)
            dcovmat(i1,i2,3) += 
                ( normconeps1*pow(d,smooth1+eps)*boost::math::cyl_bessel_k(smooth1 + eps, d) -
                  cov )/eps;
        }

        // update if same category
        if( locs(i1,dim+1) == locs(i2,dim+1) ){

            double d = 0.0;
            double d1 = 0.0;
            double d2 = 0.0;
            for(int j=0; j<dim; j++){
                d1 += pow( locs_scaled2(i1,j) - locs_scaled2(i2,j), 2.0 );
            }
	    d2 = pow( locs_scaled2(i1,dim) - locs_scaled2(i2,dim), 2.0 );
            d = pow( d1 + d2, 0.5 );
            
            double cov;        
            if( d == 0.0 ){

                cov = covparms(4);
                dcovmat(i1,i2,4) += 1.0;

            } else {

                cov = normcon2*pow( d, smooth2 )*boost::math::cyl_bessel_k( smooth2, d );
                // variance parameter
                dcovmat(i1,i2,4) += cov/covparms(4);

	        // derivatives of (d)^nu K_nu(d) with respect to covparms
	        // deriv of (d)^nu K_nu(d) with respect to d is negative!
	        // note that normcon contains covparms(4) 
                double cov_nu_m1 = normcon2*pow( d, smooth2 - 1.0 )*
                    boost::math::cyl_bessel_k(smooth2 - 1.0, d);  
	        double dM_dd = -d*cov_nu_m1;
	        
                double dd_dc5 = -1/d*d1/covparms(5);
                double dd_dc6 = -1/d*d2/covparms(6);

	        dcovmat(i1,i2,5) = dM_dd * dd_dc5;
	        dcovmat(i1,i2,6) = dM_dd * dd_dc6;

                // smoothness parameter (finite differencing)
                dcovmat(i1,i2,7) += 
                    ( normconeps2*pow(d,smooth2+eps)*boost::math::cyl_bessel_k(smooth2 + eps, d) -
                      cov )/eps;
            }
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
#endif
