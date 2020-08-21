#ifndef COVMATRIX_FUNS_01_H
#define COVMATRIX_FUNS_01_H



// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "basis.h"
#include <boost/math/special_functions.hpp>

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]


//' Isotropic exponential covariance function
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
//' \deqn{ M(x,y) = \sigma^2 exp( - || x - y ||/ \alpha )}
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat exponential_isotropic(arma::vec covparms, arma::mat locs ){

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
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        // calculate distance
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = std::sqrt( d );
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

//' @describeIn exponential_isotropic Derivatives of isotropic exponential covariance
// [[Rcpp::export]]
arma::cube d_exponential_isotropic(arma::vec covparms, arma::mat locs ){

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
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = std::sqrt( d );
        
        dcovmat(i1,i2,0) += std::exp(-d);
        dcovmat(i1,i2,1) += covparms(0)*std::exp(-d)*d/covparms(1);
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


// [[Rcpp::export]]
arma::mat exponential_isotropic_fast(arma::vec covparms, arma::mat locs ){

    std::chrono::steady_clock::time_point t1;
    std::chrono::steady_clock::time_point t2;
    std::chrono::steady_clock::time_point t3;
    std::chrono::steady_clock::time_point t4;
    std::chrono::steady_clock::time_point t5;
    std::chrono::steady_clock::time_point t6;
    std::chrono::steady_clock::time_point t7;
    std::chrono::steady_clock::time_point t8;

	t1 = std::chrono::steady_clock::now();

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

	/*
    // create scaled locations
	double ls[n][dim];
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            ls[i][j] = locs(i,j)/covparms(1);
        }
    }
	*/

	t2 = std::chrono::steady_clock::now();

    // calculate covariances
    arma::mat covmat(n,n);
	arma::mat distmat(n,n, fill::zeros);
    //double dmat[n][n];
    //double cmat[n][n];

	/*
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<n; i2++){ 
        dmat[i1][i2] = 0.0;
    }}
	*/

	t3 = std::chrono::steady_clock::now();

	/*
    for(int j=0; j<dim; j++){
		for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
				double dd = ls[i1][j]-ls[i2][j];
				dmat[i1][i2] += dd*dd;
	    }}	
    }
	*/

	t4 = std::chrono::steady_clock::now();
		
    for(int j=0; j<dim; j++){
	    //for(arma::uword i1=0; i1<n; i1++){ for(arma::uword i2=0; i2<=i1; i2++){
		for(int i1=0; i1<n; i1++){ for(int i2=0; i2<n; i2++){
				double dd = locs_scaled(i1,j)-locs_scaled(i2,j);
				distmat(i2,i1) += dd*dd;
	    }}	
    }

	/*
	for(arma::uword j=0; j<dim; j++){
        for(arma::uword i2=0; i2<n; i2++){
            double ls2 = locs_scaled(i2,j);
            arma::mat::iterator it_end = covmat.end_col(i2);
            for(arma::mat::iterator it = covmat.begin_col(i2); it != it_end; ++it ){
                (*it) += 
    }
	*/


	t5 = std::chrono::steady_clock::now();

	//covmat = distmat;

	t6 = std::chrono::steady_clock::now();

    for(arma::uword i1=0; i1<n; i1++){ for(arma::uword i2=0; i2<=i1; i2++){
		distmat(i2,i1) = std::sqrt( distmat(i2,i1) );
		covmat(i2,i1) = covparms(0)*std::exp( -distmat(i2,i1) );
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
		covmat(i1,i2) = covmat(i2,i1);
    }}

	/*
	arma::mat::iterator it_end = covmat.end();
	for(arma::mat::iterator it = covmat.begin(); it != it_end; ++it){
		(*it) = covparms(0)*std::exp( -std::sqrt( (*it) ) );
	}

	t7 = std::chrono::steady_clock::now();

	for(int i=0; i<n; i++){
        covmat(i,i) += nugget; 
    } 
	*/
	
	t8 = std::chrono::steady_clock::now();

	/*
    for(arma::uword i1=0; i1<n; i1++){ for(arma::uword i2=0; i2<n; i2++){
		distmat(i2,i1) = std::sqrt( distmat(i2,i1) );
		covmat(i2,i1) = covparms(0)*std::exp( -distmat(i2,i1) );
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
		covmat(i1,i2) = covmat(i2,i1);
    }}

	arma::mat::iterator it_end = covmat.end();
	for(arma::mat::iterator it = covmat.begin(); it != it_end; ++it){
		distmat(it) = std::sqrt( distmat(it) );
		covmat(it) = covparms(0)*std::exp( -distmat(it) );
	}
	for(int i=0; i<n; i++){
        covmat(i,i) += nugget; 
    } 
	*/

	/*
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
		dmat[i1][i2] = std::sqrt( dmat[i1][i2] );
		cmat[i1][i2] = covparms(0)*std::exp( -dmat[i1][i2] );
        if( i1 == i2 ){ cmat[i2][i2] += nugget; } 
    }}

	t7 = std::chrono::steady_clock::now();

    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
		covmat(i2,i1) = cmat[i1][i2];
        covmat(i1,i2) = cmat[i1][i2];
    }}
	*/
		
	/*
	Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() << endl;
	Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t3-t2).count() << endl;
	Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count() << endl;
	Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t5-t4).count() << endl;
	Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t6-t5).count() << endl;
	Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t7-t6).count() << endl;
	Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t8-t7).count() << endl;
	*/

    return covmat;
}

// [[Rcpp::export]]
arma::cube d_exponential_isotropic_fast(arma::vec covparms, arma::mat locs ){

    std::chrono::steady_clock::time_point t1;
    std::chrono::steady_clock::time_point t2;
    std::chrono::steady_clock::time_point t3;
    std::chrono::steady_clock::time_point t4;
    std::chrono::steady_clock::time_point t5;
    std::chrono::steady_clock::time_point t6;
	double cc = 1.0/covparms(1);
	
    int dim = locs.n_cols;
    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 2 );
    // create scaled locations

	t1 = std::chrono::steady_clock::now();
	
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }

	t2 = std::chrono::steady_clock::now();

    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
	arma::mat distmat(n,n, fill::zeros);
    for(int j=0; j<dim; j++){
	    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
			double dd = locs_scaled(i1,j)-locs_scaled(i2,j);
            distmat(i2,i1) += dd*dd;
	    }}	
    }

	t3 = std::chrono::steady_clock::now();

    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){

        distmat(i2,i1) = std::sqrt( distmat(i2,i1) );
        dcovmat(i2,i1,0) += std::exp( -distmat(i2,i1) );
        dcovmat(i2,i1,1) += covparms(0)*dcovmat(i2,i1,0)*distmat(i2,i1)*cc;
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i2,i1,0) += covparms(2);
            dcovmat(i2,i1,2) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i1,i2,j) = dcovmat(i2,i1,j);
            }
		    dcovmat(i1,i2,0) = dcovmat(i2,i1,0);
        }
    }}
	
	t4 = std::chrono::steady_clock::now();

	//Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() << endl;
	//Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t3-t2).count() << endl;
	//Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count() << endl;
	//	Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t5-t4).count() << endl;
	//	Rcout << std::chrono::duration_cast<std::chrono::microseconds>(t6-t5).count() << endl;
	
    return dcovmat;
}

//' Isotropic Matern covariance function
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs gives a point in R^d.
//' @param covparms A vector with covariance parameters
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
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern_isotropic(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*boost::math::tgamma(covparms(2) ));
    
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
                    pow( d, covparms(2) )*boost::math::cyl_bessel_k(covparms(2), d);
            }
            // add nugget
            if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
            // fill in opposite entry
            else { covmat(i1,i2) = covmat(i2,i1); }
        }    
    }
    return covmat;
}



//' @describeIn matern_isotropic Derivatives of isotropic Matern covariance
// [[Rcpp::export]]
arma::cube d_matern_isotropic(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;
    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 3 );
    double normcon = covparms(0)/(pow(2.0,covparms(2)-1.0)*boost::math::tgamma(covparms(2) ));
    double eps = 1e-8;
    double normconeps = 
        covparms(0)/(pow(2.0,covparms(2)+eps-1.0)*boost::math::tgamma(covparms(2)+eps ));
    
    // create scaled locations
    mat locs_scaled(n,dim);
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







//' Isotropic Matern covariance function, smoothness = 1.5
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
//' \deqn{ M(x,y) = \sigma^2 (1 + || x - y || ) exp( - || x - y ||/ \alpha )}
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern15_isotropic(arma::vec covparms, arma::mat locs ){

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
            covmat(i2,i1) = covparms(0)*(1 + d)*std::exp( -d );
        }
        // add nugget
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
        // fill in opposite entry
        else { covmat(i1,i2) = covmat(i2,i1); }
    }}

    return covmat;
}

//' @describeIn exponential_isotropic Derivatives of isotropic 
//' matern covariance with smoothness 1.5
// [[Rcpp::export]]
arma::cube d_matern15_isotropic(arma::vec covparms, arma::mat locs ){

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
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        dcovmat(i1,i2,0) += (1 + d)*std::exp(-d);
        dcovmat(i1,i2,1) += covparms(0)*std::exp(-d)*d*d/covparms(1);
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









#endif
