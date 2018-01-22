
#ifndef COVFUNS_H
#define COVFUNS_H


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


#endif

