#ifndef COVMATRIX_FUNS_sphere_H
#define COVMATRIX_FUNS_sphere_H

// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include "covmatrix_funs_01.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]


//' Isotropic Matern covariance function on sphere
//'
//' From a matrix of longitudes and latitudes and a vector covariance parameters of the form
//' (variance, range, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param lonlat A matrix with \code{n} rows and one column with longitudes in (-180,180)
//' and one column of latitudes in (-90,90).
//' Each row of lonlat describes a point on the sphere.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range, smoothness, nugget). Range parameter assumes that
//' the sphere has radius 1 (units are radians).
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlat[i,]} and
//' \code{lonlat[j,]}.
//' @section Matern on Sphere Domain:
//' The function first calculates the (x,y,z) 3D coordinates, and then inputs
//' the resulting locations into \code{matern_isotropic}. This means that we construct
//' covariances on the sphere by embedding the sphere in a 3D space. There has been some
//' concern expressed in the literature that such embeddings may produce distortions.
//' The source and nature of such distortions has never been articulated,
//' and to date, no such distortions have been documented. Guinness and
//' Fuentes (2016) argue that 3D embeddings produce reasonable models for data on spheres.
// [[Rcpp::export]]
arma::mat matern_sphere(arma::vec covparms, arma::mat lonlat ){

    int n = lonlat.n_rows;
    double lonrad;                                  // longitude
    double latrad;                                  // latitude
    arma::mat xyz(n, 3);
    for(int i = 0; i < n; i++){
        lonrad = 2*M_PI*lonlat(i,0)/360;
        latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }

    arma::mat covmat = matern_isotropic( covparms, xyz );
    return covmat;
}

//' @describeIn matern_sphere Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_matern_sphere(arma::vec covparms, arma::mat lonlat ){

    int n = lonlat.n_rows;
    arma::mat xyz(n, 3);
    for(int i = 0; i < n; i++){
        double lonrad = 2*M_PI*lonlat(i,0)/360;
        double latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }

    arma::cube dcovmat = d_matern_isotropic( covparms, xyz );
    return dcovmat;
}


//' Isotropic exponential covariance function on sphere
//'
//' From a matrix of longitudes and latitudes and a vector covariance parameters of the form
//' (variance, range, nugget), return the square matrix of
//' all pairwise covariances.
//' @param lonlat A matrix with \code{n} rows and one column with longitudes in (-180,180)
//' and one column of latitudes in (-90,90).
//' Each row of lonlat describes a point on the sphere.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range, nugget). Range parameter assumes that
//' the sphere has radius 1 (units are radians).
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlat[i,]} and
//' \code{lonlat[j,]}.
//' @section Covariances on spheres:
//' The function first calculates the (x,y,z) 3D coordinates, and then inputs
//' the resulting locations into \code{exponential_isotropic}. This means that we construct
//' covariances on the sphere by embedding the sphere in a 3D space. There has been some
//' concern expressed in the literature that such embeddings may produce distortions.
//' The source and nature of such distortions has never been articulated,
//' and to date, no such distortions have been documented. Guinness and
//' Fuentes (2016) argue that 3D embeddings produce reasonable models for data on spheres.
// [[Rcpp::export]]
arma::mat exponential_sphere(arma::vec covparms, arma::mat lonlat ){

    int n = lonlat.n_rows;
    double lonrad;                                  // longitude
    double latrad;                                  // latitude
    arma::mat xyz(n, 3);
    for(int i = 0; i < n; i++){
        lonrad = 2*M_PI*lonlat(i,0)/360;
        latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }

    arma::mat covmat = exponential_isotropic( covparms, xyz );
    return covmat;
}

//' @describeIn exponential_sphere Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_exponential_sphere(arma::vec covparms, arma::mat lonlat ){

    int n = lonlat.n_rows;
    arma::mat xyz(n, 3);
    for(int i = 0; i < n; i++){
        double lonrad = 2*M_PI*lonlat(i,0)/360;
        double latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }

    arma::cube dcovmat = d_exponential_isotropic( covparms, xyz );
    return dcovmat;
}


#endif
