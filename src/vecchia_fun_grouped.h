#ifndef VECCHIAFUNGROUPED_H
#define VECCHIAFUNGROUPED_H

#include "covfuns.h"
#include <Rcpp.h>
#include <assert.h>

using namespace Rcpp;

void vecchia_grouped(NumericVector covparms, StringVector covfun_name,
                                  NumericMatrix locs, List NNlist,
                                  NumericVector& y, NumericVector* Linv,
                                  NumericVector* ll, int whichreturn){

    // utility integers
    int i, j, k, el, first_ind, last_ind, first_resp, last_resp;
    int first_L_ind = 0;
    
    // number of obs
    int n = y.length();
    
    // vector of all indices
    std::vector<int> all_inds = Rcpp::as<std::vector<int> >(NNlist["all_inds"]);
    // vector of local response indices
    std::vector<int> local_resp_inds = Rcpp::as<std::vector<int> >(NNlist["local_resp_inds"]);
    // vector of global response indices
    std::vector<int> global_resp_inds = Rcpp::as<std::vector<int> >(NNlist["global_resp_inds"]);
    // last index of each block in all_inds
    std::vector<int> last_ind_of_block = Rcpp::as<std::vector<int> >(NNlist["last_ind_of_block"]);
    // last response index of each block in local_resp_inds and global_resp_inds
    std::vector<int> last_resp_of_block = Rcpp::as<std::vector<int> >(NNlist["last_resp_of_block"]);

    int nb = last_ind_of_block.size();  // number of blocks
    double d;                           // utility double
    
    //if( whichreturn == 2 ){
    //    int nentries = 0;
    //    for(j=0; j<n; j++){ nentries += local_resp_inds[j]; }
    //}    
    
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
    
    // figure out size of largest block
    // for preallocating L, g, ysub, and sig
    int maxbsize = last_ind_of_block[0];
    for(j=1; j<nb; j++){
        if( last_ind_of_block[j]-last_ind_of_block[j-1] > maxbsize ){
            maxbsize = last_ind_of_block[j]-last_ind_of_block[j-1];
        }
    }
    
    // subvector of locations
    std::vector<std::vector<double> > locsub(maxbsize, std::vector<double>(dim, 0));
    std::vector<double> ysub (maxbsize,0);   // covariance vector
    // cholesky factor
    std::vector<std::vector<double> > L(maxbsize, std::vector<double>(maxbsize, 0));
    std::vector<double> g (maxbsize,0);     // cholesky row
    std::vector<double> sig (maxbsize,0);   // covariance vector
    

    // big loop over observations starts here
    for(i=0; i<nb; i++){
        
        //Rcpp::IntegerMatrix NNi = Rcpp::as<Rcpp::IntegerMatrix>(NNlist(i));
  
        if(i==0){ first_ind = 0; } else {first_ind = last_ind_of_block[i-1]; }
        last_ind = last_ind_of_block[i];

        if(i==0){ first_resp = 0; } else {first_resp = last_resp_of_block[i-1]; }
        last_resp = last_resp_of_block[i];
        
                
        int bsize = last_ind - first_ind;
        // these are full vector indices
        std::vector<int> inds(all_inds.begin() + first_ind, 
                              all_inds.begin() + last_ind);  
        // these are subvector indices for response
        std::vector<int> whichinds(local_resp_inds.begin() + first_resp,
                                   local_resp_inds.begin() + last_resp);
        for(j=0; j < whichinds.size(); j++){ whichinds[j] -= 1; }
        
        // first, fill in ysub and locsub in forward order
        for(j=0; j<bsize; j++){
            for(k=0;k<dim;k++){ locsub[j][k] = locs( inds[j] - 1, k ); }
            ysub[j] = y( inds[j] - 1 );
        }

        // make sure all elements of Cholesky matrix are zero
        for(k=0;k<bsize;k++){ for(j=0;j<bsize;j++){ L[k][j] = 0.0; }}

        // get 0,0 entry
        L[0][0] = pow( (*p_covfun)(&locsub[0], &locsub[0], cparms) + nugget, 0.5 );

        for(j=1; j<bsize; j++){  // j = row of Li

            // initialize g
            for(k=0; k<bsize; k++){ g[k] = 0.0; }

            // get first j-1 entries of jth row of L
            for(k=0; k<j; k++){

                // compute covariance, would prefer commented version below it
                sig[k] = (*p_covfun)(&locsub[k], &locsub[j], cparms);

                // solve lower triangular system to get L[j][k]
                // might be faster for work on L[j][k] directly
                // instead of through g
                g[k] = sig[k];
                //L[j][k] = sig[k];
                if(k>0){
                    for(el=0; el<k; el++){
                        g[k] -= L[k][el]*g[el];
                        //L[j][k] -= L[k][el]*L[j][el];
                    }
                }
                g[k] = g[k]/L[k][k];
                L[j][k] = g[k];
                //L[j][k] /= L[k][k];
            }

            // get diagonal entry
            d = 0.0;
            for( k=0;k<j;k++ ){ d += g[k]*g[k]; }
            //for(k=0; k<j; k++ ){ d += L[j][k]*L[j][k]; }
            L[j][j] = pow( (*p_covfun)(&locsub[j],&locsub[j],cparms) + nugget - d, 0.5 );

        }

        if( whichreturn == 1 ){
            // get g = L^{-1}y
            g[0] = ysub[0]/L[0][0];
            for(j=1; j<bsize; j++){
                g[j] = ysub[j];
                for(k=0; k<j; k++){
                    g[j] -= L[j][k]*g[k];
                }
                g[j] = g[j]/L[j][j];
            }
            // add contribution to likelihood
            for(el=0; el<whichinds.size(); el++){
                j = whichinds[el];
                (*ll)(0) += -g[j]*g[j]/2 - log( L[j][j] );
            }
            
        } else if( whichreturn == 2 ) {
            // get L^{-1}
            //Rcpp::NumericMatrix Linvi( whichinds.size(), bsize );
            for(el=0; el<whichinds.size(); el++){
                int bsize_cur = whichinds[el] + 1;
                g[bsize_cur-1] = 1.0/L[bsize_cur-1][bsize_cur-1];
                (*Linv)(first_L_ind + bsize_cur - 1) = g[bsize_cur-1];
                for(j=bsize_cur-2; j >-1; j--){
                    g[j] = 0.0;
                    for(k = j+1; k < bsize_cur; k++){
                        g[j] += L[k][j]*g[k];
                    }
                    g[j] = -g[j]/L[j][j];
                    (*Linv)(first_L_ind + j) = g[j];
                }
                first_L_ind += bsize_cur;
            }
        }
        

     }
    // add the constant
    
    (*ll)(0) += -n*log(2*M_PI)/2;

}



#endif

