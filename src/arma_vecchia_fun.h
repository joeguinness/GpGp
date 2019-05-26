#ifndef ARMAVECCHIAFUN_H
#define ARMAVECCHIAFUN_H

#include <assert.h>
#include <RcppArmadillo.h>
#include "arma_covmatrix_funs.h"
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* 
// need to fix this, probably need to pass a pointer to the
// covariance function pointer. Too complicated for my brain

inline void update_covfun(std::string covfun_name_string, 
    arma::mat (*p_covfun)(Rcpp::NumericVector, Rcpp::NumericMatrix) ){
    
    if( covfun_name_string.compare("arma_matern_isotropic") == 0 )
    { p_covfun = &arma_matern_isotropic; } 
    else if( covfun_name_string.compare("arma_exponential_isotropic") == 0 )
    { p_covfun = &arma_exponential_isotropic; }
    else   // stop the program
    {
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
        assert(0);
    }
}
*/


void arma_vecchia(NumericVector covparms, StringVector covfun_name,
                                  const NumericMatrix locs, IntegerMatrix NNarray,
                                  NumericVector& y, NumericMatrix* Linv,
                                  NumericVector* ll, int whichreturn){

    bool trip_chol = false;

    int n = NNarray.nrow();               // length of response
    int m = NNarray.ncol();           // number of neighbors + 1
    int dim = locs.ncol();            

    // the user inputs a covariance function name as a string (covfun_name)
    // then based on the string, we assign a pointer to an allowable
    // c++ covariance function
    double cparms[3];       // non-nugget cov parameters (defined for matern)
    double nugget;
    
    mat (*p_covfun)(NumericVector, NumericMatrix);
    cube (*p_d_covfun)(NumericVector, NumericMatrix);

    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];
    
    // would prefer to have this in a separate function
    // that function is not working right now
    if( covfun_name_string.compare("arma_matern_isotropic") == 0 )
    { p_covfun = &arma_matern_isotropic; } 
    else if( covfun_name_string.compare("arma_exponential_isotropic") == 0 )
    { 
        p_covfun = &arma_exponential_isotropic; 
        p_d_covfun = &d_arma_exponential_isotropic;
    }
    else if( covfun_name_string.compare("arma_matern_scaledim") == 0 )
    { p_covfun = &arma_matern_scaledim; }
    else { // stop the program
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
        assert(0);
    }

    // big loop over observations starts here
    for(int i=0; i<n; i++){

        int bsize = std::min(i+1,m);
        vec ysub(bsize);                   // subset of response values
        NumericMatrix locsub(bsize, dim);
        
        Rcpp::checkUserInterrupt();
        
        // first, fill in ysub and locsub in reverse order
        for(int j=bsize-1; j>=0; j--){
            for(int k=0;k<dim;k++){ locsub(bsize-1-j,k) = locs( NNarray(i,j)-1, k ); }
            ysub(bsize-1-j) = y( NNarray(i,j)-1 );
        }
        
        // compute covariance matrix
        mat covmat = (*p_covfun)( covparms, locsub );

        // take cholesky
        mat cholmat = eye( size(covmat) );
        chol( cholmat, covmat, "lower");

        // do solve
        vec z = solve( trimatl(cholmat), ysub );

        // calculate contribution to likelihood
        (*ll)(0) += -z(bsize-1)*z(bsize-1)*0.5 - log( cholmat(bsize-1,bsize-1) );

    }

    // add the constant
    (*ll)(0) += -n*log(2*M_PI)/2;
    if( trip_chol ){
        Rcpp::Rcout << "*" << std::endl;
    }

}



#endif

