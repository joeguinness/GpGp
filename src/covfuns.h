
#ifndef COVFUNS_H
#define COVFUNS_H

#include <assert.h>

inline double matern_function(double d, double *cparms){

    // has special cases for 1/2 and 3/2
    if( d == 0.0 ){
        d = cparms[0];
    } else {
        if( cparms[2] == 0.5 ){
            d = cparms[0]*exp(-d/cparms[1]);
        } else if( cparms[2] == 1.5 ){
            d = cparms[0]*(1+d/cparms[1])*exp(-d/cparms[1]);
        } else {
            double normcon = cparms[0]/(pow(2.0,cparms[2]-1)*Rf_gammafn(cparms[2]));
            d = normcon*pow( d/cparms[1], cparms[2] )*
            Rf_bessel_k(d/cparms[1],cparms[2],1.0);
        }
    }
    return d;
}

inline double matern_isotropic_internal( const std::vector<double>* loc1, const std::vector<double>* loc2, double* cparms){

    // inputs two vectors loc1 and loc2 of length dim and returns
    // isotropic matern covariance using euclidean distance
    // between loc1 and loc2
    int dim = (*loc1).size();
    double d = 0.0;
    int i;
    for(i=0; i<dim; i++){
        d += pow( (*loc1)[i] - (*loc2)[i], 2);
    }
    d = pow(d, 0.5);

    d = matern_function(d,cparms);
    return d;

}

inline void update_vars_based_on_covfun(std::string covfun_name_string, 
    double cparms[], double* nugget, Rcpp::NumericMatrix* locs_scaled, 
    const Rcpp::NumericVector covparms){
    
    int n = (*locs_scaled).nrow();
    
        // set p_covfun, cparms, and locations based on covfun_name_string
    if( covfun_name_string.compare("matern_isotropic") == 0 )
    {
        for(int k=0; k<3; k++){ cparms[k] = covparms[k]; }  // re-assign non-nugget parameters
        *nugget = covparms[0]*covparms[3];               // separate variable for nugget
    }
    else if( covfun_name_string.compare("matern_sphere") == 0 )
    {
        for(int k=0; k<3; k++){ cparms[k] = covparms[k]; }  // re-assign non-nugget parameters
        *nugget = covparms[0]*covparms[3];               // separate variable for nugget
        double lonrad;                                  // longitude
        double latrad;                                  // latitude
        Rcpp::NumericMatrix xyz(n, 3);
        for(int i = 0; i < n; i++){
            lonrad = 2*M_PI*(*locs_scaled)(i,0)/360;
            latrad = 2*M_PI*((*locs_scaled)(i,1)+90)/360;
            xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
            xyz(i,1) = sin(latrad)*sin(lonrad);
            xyz(i,2) = cos(latrad);
        }
        *locs_scaled = xyz;
    }
    else if( covfun_name_string.compare("matern_sphere_time") == 0 )
    {
        cparms[0] = covparms[0];                    // variance
        cparms[1] = 1;                              // locations scaled below, so set range = 1
        cparms[2] = covparms[3];                    // smoothness
        *nugget = covparms[0]*covparms[4];           // nugget
        double lonrad;
        double latrad;
        Rcpp::NumericMatrix xyzt(n, 4);
        for(int i = 0; i < n; i++){
            lonrad = 2*M_PI*(*locs_scaled)(i,0)/360;
            latrad = 2*M_PI*((*locs_scaled)(i,1)+90)/360;
            xyzt(i,0) = sin(latrad)*cos(lonrad)/covparms[1];   // convert lon,lat,time to
            xyzt(i,1) = sin(latrad)*sin(lonrad)/covparms[1];   // scaled x,y,z, and scaled time
            xyzt(i,2) = cos(latrad)/covparms[1];
            xyzt(i,3) = (*locs_scaled)(i,2)/covparms[2];
        }
        *locs_scaled = xyzt;
    }
    else if( covfun_name_string.compare("matern_space_time") == 0 )
    {
        int d = (*locs_scaled).ncol() - 1;
        cparms[0] = covparms[0];                    // variance
        cparms[1] = 1;                              // locations scaled below, so set range = 1
        cparms[2] = covparms[3];                    // smoothness
        *nugget = covparms[0]*covparms[4];          // nugget
        for(int i = 0; i < n; i++){
            for(int j=0; j<d; j++){
                (*locs_scaled)(i,j) = (*locs_scaled)(i,j)/covparms[1];
            }
            (*locs_scaled)(i,d) = (*locs_scaled)(i,d)/covparms[2];
        }

    }
    else   // stop the program
    {
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
        assert(0);
    }
    
    
}

#endif

