#ifndef COVMATRIX_FUNS_spherewarp_H
#define COVMATRIX_FUNS_spherewarp_H

// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include "basis.h"
#include "covmatrix_funs_01.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]



//' Deformed Matern covariance function on sphere
//'
//' From a matrix of longitudes and latitudes and a vector covariance parameters of the form
//' (variance, range, smoothness, nugget, <5 warping parameters>), return the square matrix of
//' all pairwise covariances.
//' @param lonlat A matrix with \code{n} rows and one column with longitudes in (-180,180)
//' and one column of latitudes in (-90,90).
//' Each row of lonlat describes a point on the sphere.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range, smoothness, nugget, <5 warping parameters>). 
//' Range parameter assumes that the sphere has radius 1 (units are radians).
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlat[i,]} and
//' \code{lonlat[j,]}.
//' @section Warpings:
//' The function first calculates the (x,y,z) 3D coordinates, and then "warps"
//' the locations to \eqn{(x,y,z) + \Phi(x,y,z)}, where \eqn{\Phi} is a warping
//' function composed of gradients of spherical harmonic functions of degree 2.
//' See Guinness (2019, "Gaussian Process Learning via Fisher Scoring of 
//' Vecchia's Approximation") for details.
//' The warped locations are input into \code{matern_isotropic}. 
// [[Rcpp::export]]
arma::mat matern_sphere_warp(arma::vec covparms, arma::mat lonlat ){

    int n = lonlat.n_rows;
    int nisoparms = 4;
    arma::vec isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.n_elem - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    arma::mat xyz(n, 3);
    for(int i = 0; i < n; i++){
        double lonrad = 2*M_PI*lonlat(i,0)/360;
        double latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }
    
    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyz, Lmax );

    // warp locations
    for(int i=0; i<n; i++){
        for(int k=0; k<3; k++){
            for(int j=0; j<nbasis; j++){
                xyz(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
            }
        }
    }
    
    // compute covariances
    arma::mat covmat = matern_isotropic( isoparms, xyz );
    return covmat;
}

//' @describeIn matern_sphere_warp Derivatives with respect to parameters.
// [[Rcpp::export]]
arma::cube d_matern_sphere_warp(arma::vec covparms, arma::mat lonlat ){

    int n = lonlat.n_rows;
    int nisoparms = 4;
    arma::vec isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.n_elem - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    arma::mat xyz(n, 3);
    for(int i = 0; i < n; i++){
        double lonrad = 2*M_PI*lonlat(i,0)/360;
        double latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }
    
    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyz, Lmax );

    // warp locations
    for(int i=0; i<n; i++){
        for(int k=0; k<3; k++){
            for(int j=0; j<nbasis; j++){
                xyz(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
            }
        }
    }

    arma::cube ddcov = d_matern_isotropic( isoparms, xyz );
    arma::cube dcovmat(n,n,covparms.n_elem);
    for(int i=0; i<nisoparms; i++){ dcovmat.slice(i) = ddcov.slice(i); }
    for(int j=0; j<nbasis; j++){
        for(int i1=0; i1<n; i1++){ for(int i2=i1; i2<n; i2++){
            // compute distance
            double dd = 0.0;
            for(int k=0; k<3; k++){ dd += pow( xyz(i2,k) - xyz(i1,k), 2.0); }
            dd = pow(dd,0.5);
            if(dd==0.0){ 
                dcovmat(i2,i1,nisoparms+j) = 0.0; 
            } else {
                // make use of already computed derivatives wrt covparms(1)
                dcovmat(i2,i1,nisoparms+j) = -dcovmat(i2,i1,1)*covparms(1)/dd;
                // pd wrt to location component * pd wrt to basis
                double pd = 0.0;
                for(int k=0; k<3; k++){ 
                    pd += (xyz(i2,k)-xyz(i1,k))/dd * grad_basis(i2,j,k);
                    pd += (xyz(i1,k)-xyz(i2,k))/dd * grad_basis(i1,j,k);
                }
                dcovmat(i2,i1,nisoparms+j) *= pd;
            }
            // fill in opposite side
            dcovmat(i1,i2,nisoparms+j) = dcovmat(i2,i1,nisoparms+j);
        }}
    }
    return dcovmat;
}





//' Deformed exponential covariance function on sphere
//'
//' From a matrix of longitudes and latitudes and a vector covariance parameters of the form
//' (variance, range, nugget, <5 warping parameters>), return the square matrix of
//' all pairwise covariances.
//' @param lonlat A matrix with \code{n} rows and one column with longitudes in (-180,180)
//' and one column of latitudes in (-90,90).
//' Each row of lonlat describes a point on the sphere.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range, nugget, <5 warping parameters>). 
//' Range parameter assumes that the sphere has radius 1 (units are radians).
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlat[i,]} and
//' \code{lonlat[j,]}.
//' @section Warpings:
//' The function first calculates the (x,y,z) 3D coordinates, and then "warps"
//' the locations to \eqn{(x,y,z) + \Phi(x,y,z)}, where \eqn{\Phi} is a warping
//' function composed of gradients of spherical harmonic functions of degree 2.
//' See Guinness (2019, "Gaussian Process Learning via Fisher Scoring of 
//' Vecchia's Approximation") for details.
//' The warped locations are input into \code{exponential_isotropic}. 
// [[Rcpp::export]]
arma::mat exponential_sphere_warp(arma::vec covparms, arma::mat lonlat ){

    int n = lonlat.n_rows;
    int nisoparms = 3;
    arma::vec isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.n_elem - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    arma::mat xyz(n, 3);
    for(int i = 0; i < n; i++){
        double lonrad = 2*M_PI*lonlat(i,0)/360;
        double latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }
    
    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyz, Lmax );

    // warp locations
    for(int i=0; i<n; i++){
        for(int k=0; k<3; k++){
            for(int j=0; j<nbasis; j++){
                xyz(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
            }
        }
    }
    
    // compute covariances
    arma::mat covmat = exponential_isotropic( isoparms, xyz );
    return covmat;
}

//' @describeIn exponential_sphere_warp Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_exponential_sphere_warp(arma::vec covparms, arma::mat lonlat ){

    int n = lonlat.n_rows;
    int nisoparms = 3;
    arma::vec isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.n_elem - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    arma::mat xyz(n, 3);
    for(int i = 0; i < n; i++){
        double lonrad = 2*M_PI*lonlat(i,0)/360;
        double latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }
    
    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyz, Lmax );

    // warp locations
    for(int i=0; i<n; i++){
        for(int k=0; k<3; k++){
            for(int j=0; j<nbasis; j++){
                xyz(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
            }
        }
    }

    arma::cube ddcov = d_exponential_isotropic( isoparms, xyz );
    arma::cube dcovmat(n,n,covparms.n_elem);
    for(int i=0; i<nisoparms; i++){ dcovmat.slice(i) = ddcov.slice(i); }
    for(int j=0; j<nbasis; j++){
        for(int i1=0; i1<n; i1++){ for(int i2=i1; i2<n; i2++){
            // compute distance
            double dd = 0.0;
            for(int k=0; k<3; k++){ dd += pow( xyz(i2,k) - xyz(i1,k), 2.0); }
            dd = pow(dd,0.5);
            if(dd==0.0){ 
                dcovmat(i2,i1,nisoparms+j) = 0.0; 
            } else {
                // make use of already computed derivatives wrt covparms(1)
                dcovmat(i2,i1,nisoparms+j) = -dcovmat(i2,i1,1)*covparms(1)/dd;
                // pd wrt to location component * pd wrt to basis
                double pd = 0.0;
                for(int k=0; k<3; k++){ 
                    pd += (xyz(i2,k)-xyz(i1,k))/dd * grad_basis(i2,j,k);
                    pd += (xyz(i1,k)-xyz(i2,k))/dd * grad_basis(i1,j,k);
                }
                dcovmat(i2,i1,nisoparms+j) *= pd;
            }
            // fill in opposite side
            dcovmat(i1,i2,nisoparms+j) = dcovmat(i2,i1,nisoparms+j);
        }}
    }
    return dcovmat;
}

#endif
