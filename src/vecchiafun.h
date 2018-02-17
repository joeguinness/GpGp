#ifndef VECCHIAFUN_H
#define VECCHIAFUN_H

#include "covfuns.h"
#include <Rcpp.h>
#include <assert.h>

using namespace Rcpp;

void vecchia(NumericVector covparms, StringVector covfun_name,
                                  const NumericMatrix locs, IntegerMatrix NNarray,
                                  NumericVector& y, NumericMatrix* Linv,
                                  NumericVector* ll, int whichreturn){

    // utility integers
    int i, j, k, el;
    bool trip_chol = false;

    int n = NNarray.nrow();               // length of response
    int m = NNarray.ncol();           // number of neighbors + 1
    double d;                         // utility double
//    double ysub[m];                   // subset of response values
    std::vector<double> ysub (m,0);                   // subset of response values

    // cholesky factor
    std::vector<std::vector<double> > L(m, std::vector<double>(m, 0));
    std::vector<double> g (m,0);     // cholesky row
    std::vector<double> sig (m,0);   // covariance vector

    // the user inputs a covariance function name as a string (covfun_name)
    // then based on the string, we assign a pointer to an allowable
    // c++ covariance function
    double cparms[3];       // non-nugget cov parameters
    double nugget;
    double (*p_covfun)(const std::vector<double>* loc1, const std::vector<double>* loc2, double* cparms);

    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];
    
    // update everything based on the covariance functin name
    NumericMatrix locs_scaled( locs.nrow(), locs.ncol() );
    for(i=0; i<locs.nrow(); i++){
        for(j=0; j<locs.ncol(); j++){
            locs_scaled(i,j) = locs(i,j);
        }
    }
    update_vars_based_on_covfun(covfun_name_string, cparms, &nugget, &locs_scaled, covparms); 
    p_covfun = &matern_isotropic_internal; // pointer to covariance fun
    // only matern functions implemented so far


    int dim = locs_scaled.ncol();            // dimension of newly defined locations
    // subvector of locations
    std::vector<std::vector<double> > locsub(m, std::vector<double>(dim, 0));

    // big loop over observations starts here
    for(i=0; i<n; i++){

        int bsize = std::min(i+1,m);

        // first, fill in ysub and locsub in reverse order
        for(j=bsize-1; j>=0; j--){
            for(k=0;k<dim;k++){ locsub[bsize-1-j][k] = locs_scaled( NNarray(i,j)-1, k ); }
            ysub[bsize-1-j] = y( NNarray(i,j)-1 );
        }

        // make sure all elements of Cholesky matrix are zero
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ L[k][j] = 0.0; }}

        // get 0,0 entry
        L[0][0] = pow( (*p_covfun)(&locsub[0], &locsub[0], cparms) + nugget, 0.5 );

        for(j=2; j<bsize+1; j++){  // j = row of Li

            // initialize g
            for(k=0; k<m; k++){ g[k] = 0.0; }

            // get first j-1 entries of jth row of L
            for(k=0; k<j-1; k++){

                // compute covariance, would prefer commented version below it
                sig[k] = (*p_covfun)(&locsub[k], &locsub[j-1], cparms);

                // solve lower triangular system to get L[j-1][k]
                // might be faster for work on L[j-1][k] directly
                // instead of through g
                g[k] = sig[k];
                if(k>0){
                    for(el=0; el<k; el++){
                        g[k] -= L[k][el]*g[el];
                    }
                }
                g[k] = g[k]/L[k][k];
                L[j-1][k] = g[k];
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            d = (*p_covfun)(&locsub[j-1],&locsub[j-1],cparms) + nugget - d;
            if( d <= 0 ){ 
                d = cparms[0]*1000000; 
                trip_chol = true;
            }
            L[j-1][j-1] = pow( d, 0.5 );
        }

        // get g = L^{-1}y
        if( whichreturn == 1 ){
            g[0] = ysub[0]/L[0][0];
            for(j=1; j<bsize; j++){
                g[j] = ysub[j];
                for(k=0; k<j; k++){
                    g[j] -= L[j][k]*g[k];
                }
                g[j] = g[j]/L[j][j];
            }
        } else if( whichreturn == 2 ) {

            g[bsize-1] = 1.0/L[bsize-1][bsize-1];
            (*Linv)(i, (bsize-1) - (bsize-1) ) = g[bsize-1];
            for(j=bsize-2; j >-1; j--){
                g[j] = 0.0;
                for(k = j+1; k < bsize; k++){
                    g[j] += L[k][j]*g[k];
                }
                g[j] = -g[j]/L[j][j];
               (*Linv)(i, (bsize-1) - j ) = g[j];
            }

        }


        // get contribution to likelihood
        (*ll)(0) += -g[bsize-1]*g[bsize-1]/2 - log( L[bsize-1][bsize-1] );
    }

    // add the constant
    (*ll)(0) += -n*log(2*M_PI)/2;
    if( trip_chol ){
        Rcpp::Rcout << "*" << std::endl;
    }

}



#endif

