#ifndef VECCHIAFUN_H
#define VECCHIAFUN_H

#include "covfuns.h"
#include <Rcpp.h>
#include <assert.h>

using namespace Rcpp;

void vecchia(NumericVector covparms, StringVector covfun_name,
                                  NumericMatrix locs, IntegerMatrix NNarray,
                                  NumericVector& y, NumericMatrix* Linv,
                                  NumericVector* ll, int whichreturn){

    // utility integers
    int i, j, k, el;

    int n = NNarray.nrow();               // length of response
    int m = NNarray.ncol();           // number of neighbors + 1
    double d;                         // utility double
    double ysub[m];                   // subset of response values

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

    // set p_covfun, cparms, and locations based on covfun_name_string
    if( covfun_name_string.compare("matern_isotropic") == 0 )
    {
        for(k=0; k<3; k++){ cparms[k] = covparms[k]; }  // re-assign non-nugget parameters
        nugget = covparms[0]*covparms[3];               // separate variable for nugget
        p_covfun = &matern_isotropic_internal;                    // pointer to covariance fun
    }
    else if( covfun_name_string.compare("matern_sphere") == 0 )
    {
        for(k=0; k<3; k++){ cparms[k] = covparms[k]; }  // re-assign non-nugget parameters
        nugget = covparms[0]*covparms[3];               // separate variable for nugget
        double lonrad;                                  // longitude
        double latrad;                                  // latitude
        Rcpp::NumericMatrix xyz(n, 3);
        for(i = 0; i < n; i++){
            lonrad = 2*M_PI*locs(i,0)/360;
            latrad = 2*M_PI*(locs(i,1)+90)/360;
            xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
            xyz(i,1) = sin(latrad)*sin(lonrad);
            xyz(i,2) = cos(latrad);
        }
        locs = xyz;
        p_covfun = &matern_isotropic_internal; 
    }
    else if( covfun_name_string.compare("matern_sphere_time") == 0 )
    {
        cparms[0] = covparms[0];                    // variance
        cparms[1] = 1;                              // locations scaled below, so set range = 1
        cparms[2] = covparms[3];                    // smoothness
        nugget = covparms[0]*covparms[4];           // nugget
        double lonrad;
        double latrad;
        Rcpp::NumericMatrix xyzt(n, 4);
        for(i = 0; i < n; i++){
            lonrad = 2*M_PI*locs(i,0)/360;
            latrad = 2*M_PI*(locs(i,1)+90)/360;
            xyzt(i,0) = sin(latrad)*cos(lonrad)/covparms[1];   // convert lon,lat,time to
            xyzt(i,1) = sin(latrad)*sin(lonrad)/covparms[1];   // scaled x,y,z, and scaled time
            xyzt(i,2) = cos(latrad)/covparms[1];
            xyzt(i,3) = locs(i,2)/covparms[2];
        }
        locs = xyzt;
        p_covfun = &matern_isotropic_internal; 
    }
    else   // stop the program
    {
        std::cout << "Unrecognized Covariance Function Name \n";
        assert(0);
    }

    int dim = locs.ncol();            // dimension of newly defined locations
    // subvector of locations
    std::vector<std::vector<double> > locsub(m, std::vector<double>(dim, 0));

    // big loop over observations starts here
    for(i=0; i<n; i++){

        int bsize = std::min(i+1,m);

        // first, fill in ysub and locsub in reverse order
        for(j=bsize-1; j>=0; j--){
            for(k=0;k<dim;k++){ locsub[bsize-1-j][k] = locs( NNarray(i,j)-1, k ); }
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
            L[j-1][j-1] = pow( (*p_covfun)(&locsub[j-1],&locsub[j-1],cparms) + nugget - d, 0.5 );
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

}



#endif

