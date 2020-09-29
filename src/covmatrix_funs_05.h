#ifndef COVMATRIX_FUNS_spheretime_H
#define COVMATRIX_FUNS_spheretime_H

// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include "basis.h"
#include "covmatrix_funs_01.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]


//' Matern covariance function on sphere x time
//'
//' From a matrix of longitudes, latitudes, and times, and a vector covariance parameters of the form
//' (variance, range_1, range_2, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param lonlattime A matrix with \code{n} rows and three columns: longitudes in (-180,180),
//' latitudes in (-90,90), and times.
//' Each row of lonlattime describes a point on the sphere x time.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range_1, range_2, smoothness, nugget), where range_1 is a 
//' spatial range (assuming sphere of radius 1), and range_2 is a temporal range.
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlattime[i,]} and
//' \code{lonlattime[j,]}.
//' @section Covariances on spheres:
//' The function first calculates the (x,y,z) 3D coordinates, and then inputs
//' the resulting locations into \code{matern_spacetime}. This means that we construct
//' covariances on the sphere by embedding the sphere in a 3D space. There has been some
//' concern expressed in the literature that such embeddings may produce distortions.
//' The source and nature of such distortions has never been articulated,
//' and to date, no such distortions have been documented. Guinness and
//' Fuentes (2016) argue that 3D embeddings produce reasonable models for data on spheres.
// [[Rcpp::export]]
arma::mat matern_spheretime(arma::vec covparms, arma::mat lonlattime ){
    
    int n = lonlattime.n_rows;
    // matrix to hold (x,y,z,t)
    arma::mat locs(n, 4);
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

//' @describeIn matern_spheretime Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_matern_spheretime(arma::vec covparms, arma::mat lonlattime){
    
    int n = lonlattime.n_rows;
    // matrix to hold (x,y,z,t)
    arma::mat locs(n, 4);
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

//' Exponential covariance function on sphere x time
//'
//' From a matrix of longitudes, latitudes, and times, and a vector covariance parameters of the form
//' (variance, range_1, range_2, nugget), return the square matrix of
//' all pairwise covariances.
//' @param lonlattime A matrix with \code{n} rows and three columns: longitudes in (-180,180),
//' latitudes in (-90,90), and times.
//' Each row of lonlattime describes a point on the sphere x time.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range_1, range_2, nugget), where range_1 is a 
//' spatial range (assuming sphere of radius 1), and range_2 is a temporal range.
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlattime[i,]} and
//' \code{lonlattime[j,]}.
//' @section Covariances on spheres:
//' The function first calculates the (x,y,z) 3D coordinates, and then inputs
//' the resulting locations into \code{exponential_spacetime}. This means that we construct
//' covariances on the sphere by embedding the sphere in a 3D space. There has been some
//' concern expressed in the literature that such embeddings may produce distortions.
//' The source and nature of such distortions has never been articulated,
//' and to date, no such distortions have been documented. Guinness and
//' Fuentes (2016) argue that 3D embeddings produce reasonable models for data on spheres.
// [[Rcpp::export]]
arma::mat exponential_spheretime(arma::vec covparms, arma::mat lonlattime ){
    
    int n = lonlattime.n_rows;
    // matrix to hold (x,y,z,t)
    arma::mat locs(n, 4);
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
    arma::mat covmat = exponential_spacetime( covparms, locs );
    return covmat;
}

//' @describeIn exponential_spheretime Derivatives with respect to parameters.
// [[Rcpp::export]]
arma::cube d_exponential_spheretime(arma::vec covparms, arma::mat lonlattime){
    
    int n = lonlattime.n_rows;
    // matrix to hold (x,y,z,t)
    arma::mat locs(n, 4);
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
    arma::cube dcovmat = d_exponential_spacetime( covparms, locs );
    return dcovmat;
}



#endif
