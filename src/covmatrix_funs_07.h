

#ifndef COVMATRIX_FUNS_spheretimewarp_H
#define COVMATRIX_FUNS_spheretimewarp_H

// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include "basis.h"
#include "covmatrix_funs_03.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]





//' Deformed Matern covariance function on sphere
//'
//' From a matrix of longitudes, latitudes, times, and a vector covariance parameters of the form
//' (variance, range_1, range_2, smoothness, nugget, <5 warping parameters>), return the square matrix of
//' all pairwise covariances.
//' @param lonlattime A matrix with \code{n} rows and three columns: longitudes in (-180,180),
//' latitudes in (-90,90), and times.
//' Each row of lonlattime describes a point on the sphere x time.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range_1, range_2, smoothness, nugget, <5 warping parameters>). 
//' range_1 is a spatial range parameter that assumes that the sphere 
//' has radius 1 (units are radians). range_2 is a temporal range parameter.
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlat[i,]} and
//' \code{lonlat[j,]}.
//' @section Warpings:
//' The function first calculates the (x,y,z) 3D coordinates, and then "warps"
//' the locations to \eqn{(x,y,z) + \Phi(x,y,z)}, where \eqn{\Phi} is a warping
//' function composed of gradients of spherical harmonic functions of degree 2.
//' See Guinness (2019, "Gaussian Process Learning via Fisher Scoring of 
//' Vecchia's Approximation") for details.
//' The warped locations are input into \code{matern_spacetime}. The function
//' does not do temporal warping.
// [[Rcpp::export]]
arma::mat matern_spheretime_warp(arma::vec covparms, arma::mat lonlattime ){

    int n = lonlattime.n_rows;
    int nisoparms = 5;
    arma::vec isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.n_elem - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    arma::mat xyzt(n, 4);
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattime(i,0)/360;
        double latrad = 2*M_PI*(lonlattime(i,1)+90)/360;
        xyzt(i,0) = sin(latrad)*cos(lonrad); // convert lon,lat to x,y,z
        xyzt(i,1) = sin(latrad)*sin(lonrad);
        xyzt(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){ xyzt(i,3) = lonlattime(i,2); }

    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyzt , Lmax );

    // warp locations
    for(int i=0; i<n; i++){ for(int k=0; k<3; k++){ for(int j=0; j<nbasis; j++){
        xyzt(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
    }}}

    // compute covariances
    arma::mat covmat = matern_spacetime( isoparms, xyzt );
    return covmat;
}

//' @describeIn matern_spheretime_warp Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_matern_spheretime_warp(arma::vec covparms, arma::mat lonlattime ){

    int n = lonlattime.n_rows;
    int nisoparms = 5;
    arma::vec isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.n_elem - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    arma::mat xyzt(n, 4);
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattime(i,0)/360;
        double latrad = 2*M_PI*(lonlattime(i,1)+90)/360;
        xyzt(i,0) = sin(latrad)*cos(lonrad); // convert lon,lat to x,y,z
        xyzt(i,1) = sin(latrad)*sin(lonrad);
        xyzt(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){ xyzt(i,3) = lonlattime(i,2); }

    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyzt , Lmax );
    
    // warp locations
    for(int i=0; i<n; i++){ for(int k=0; k<3; k++){ for(int j=0; j<nbasis; j++){
        xyzt(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
    }}}
    
    arma::cube ddcov = d_matern_spacetime( isoparms, xyzt );
    arma::cube dcovmat(n,n,covparms.n_elem);
    for(int i=0; i<nisoparms; i++){ dcovmat.slice(i) = ddcov.slice(i); }

    for(int j=0; j<nbasis; j++){
        for(int i1=0; i1<n; i1++){ for(int i2=i1; i2<n; i2++){
            // compute distance
            double dd = 0.0;
            for(int k=0; k<3; k++){ 
                dd += pow( (xyzt(i2,k) - xyzt(i1,k))/covparms(1), 2.0); 
            }
            dd += pow( (xyzt(i2,3)-xyzt(i1,3))/covparms(2), 2.0 );
            dd = pow(dd,0.5);
            if(dd==0.0){ 
                dcovmat(i2,i1,nisoparms+j) = 0.0; 
            } else {
                // make use of already computed derivatives wrt covparms(1)
                double pth = 0.0;
                // derivative of distance wrt to temporal range parameter
                pth += -(xyzt(i1,3)-xyzt(i2,3))/pow(covparms(2),3)/dd*xyzt(i1,3);
                pth += -(xyzt(i2,3)-xyzt(i1,3))/pow(covparms(2),3)/dd*xyzt(i2,3);
                // derivative of distance wrt to spatial range parameter
                for(int k=0; k<3; k++){ 
                    pth += -(xyzt(i1,k)-xyzt(i2,k))/pow(covparms(1),3)/dd * xyzt(i1,k);
                    pth += -(xyzt(i2,k)-xyzt(i1,k))/pow(covparms(1),3)/dd * xyzt(i2,k);
                }
                dcovmat(i2,i1,nisoparms+j) = (dcovmat(i2,i1,1)+dcovmat(i2,i1,2))/pth;
                // pd wrt to location component * pd wrt to basis
                double pd = 0.0;
                for(int k=0; k<3; k++){ 
                    pd += (xyzt(i2,k)-xyzt(i1,k))/pow(covparms(1),2)/dd * grad_basis(i2,j,k);
                    pd += (xyzt(i1,k)-xyzt(i2,k))/pow(covparms(1),2)/dd * grad_basis(i1,j,k);
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
//' From a matrix of longitudes, latitudes, times, and a vector covariance parameters of the form
//' (variance, range_1, range_2, nugget, <5 warping parameters>), return the square matrix of
//' all pairwise covariances.
//' @param lonlattime A matrix with \code{n} rows and three columns: longitudes in (-180,180),
//' latitudes in (-90,90), and times.
//' Each row of lonlattime describes a point on the sphere x time.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range_1, range_2, nugget, <5 warping parameters>). 
//' range_1 is a spatial range parameter that assumes that the sphere 
//' has radius 1 (units are radians). range_2 is a temporal range parameter.
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlat[i,]} and
//' \code{lonlat[j,]}.
//' @section Warpings:
//' The function first calculates the (x,y,z) 3D coordinates, and then "warps"
//' the locations to \eqn{(x,y,z) + \Phi(x,y,z)}, where \eqn{\Phi} is a warping
//' function composed of gradients of spherical harmonic functions of degree 2.
//' See Guinness (2019, "Gaussian Process Learning via Fisher Scoring of 
//' Vecchia's Approximation") for details.
//' The warped locations are input into \code{exponential_spacetime}. The function
//' does not do temporal warping.
// [[Rcpp::export]]
arma::mat exponential_spheretime_warp(arma::vec covparms, arma::mat lonlattime ){

    int n = lonlattime.n_rows;
    int nisoparms = 4;
    arma::vec isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.n_elem - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    arma::mat xyzt(n, 4);
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattime(i,0)/360;
        double latrad = 2*M_PI*(lonlattime(i,1)+90)/360;
        xyzt(i,0) = sin(latrad)*cos(lonrad); // convert lon,lat to x,y,z
        xyzt(i,1) = sin(latrad)*sin(lonrad);
        xyzt(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){ xyzt(i,3) = lonlattime(i,2); }

    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyzt , Lmax );

    // warp locations
    for(int i=0; i<n; i++){ for(int k=0; k<3; k++){ for(int j=0; j<nbasis; j++){
        xyzt(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
    }}}

    // compute covariances
    arma::mat covmat = exponential_spacetime( isoparms, xyzt );
    return covmat;
}

//' @describeIn exponential_spheretime_warp Derivatives with respect to parameters
// [[Rcpp::export]]
arma::cube d_exponential_spheretime_warp(arma::vec covparms, arma::mat lonlattime ){

    int n = lonlattime.n_rows;
    int nisoparms = 4;
    arma::vec isoparms(nisoparms);
    for(int i=0; i<nisoparms; i++){ isoparms(i) = covparms(i); }
    int nbasis = covparms.n_elem - nisoparms;
    int Lmax = pow( nbasis + 4, 0.5 ) - 1;
    
    arma::mat xyzt(n, 4);
    for(int i=0; i<n; i++){
        double lonrad = 2*M_PI*lonlattime(i,0)/360;
        double latrad = 2*M_PI*(lonlattime(i,1)+90)/360;
        xyzt(i,0) = sin(latrad)*cos(lonrad); // convert lon,lat to x,y,z
        xyzt(i,1) = sin(latrad)*sin(lonrad);
        xyzt(i,2) = cos(latrad);
    }
    // plug in time
    for(int i=0; i<n; i++){ xyzt(i,3) = lonlattime(i,2); }

    // generate warping basis
    arma::cube grad_basis = sph_grad_xyz( xyzt , Lmax );
    
    // warp locations
    for(int i=0; i<n; i++){ for(int k=0; k<3; k++){ for(int j=0; j<nbasis; j++){
        xyzt(i,k) += covparms(nisoparms+j)*grad_basis(i,j,k);
    }}}
    
    arma::cube ddcov = d_exponential_spacetime( isoparms, xyzt );
    arma::cube dcovmat(n,n,covparms.n_elem);
    for(int i=0; i<nisoparms; i++){ dcovmat.slice(i) = ddcov.slice(i); }

    for(int j=0; j<nbasis; j++){
        for(int i1=0; i1<n; i1++){ for(int i2=i1; i2<n; i2++){
            // compute distance
            double dd = 0.0;
            for(int k=0; k<3; k++){ 
                dd += pow( (xyzt(i2,k) - xyzt(i1,k))/covparms(1), 2.0); 
            }
            dd += pow( (xyzt(i2,3)-xyzt(i1,3))/covparms(2), 2.0 );
            dd = pow(dd,0.5);
            if(dd==0.0){ 
                dcovmat(i2,i1,nisoparms+j) = 0.0; 
            } else {
                // make use of already computed derivatives wrt covparms(1)
                double pth = 0.0;
                // derivative of distance wrt to temporal range parameter
                pth += -(xyzt(i1,3)-xyzt(i2,3))/pow(covparms(2),3)/dd*xyzt(i1,3);
                pth += -(xyzt(i2,3)-xyzt(i1,3))/pow(covparms(2),3)/dd*xyzt(i2,3);
                // derivative of distance wrt to spatial range parameter
                for(int k=0; k<3; k++){ 
                    pth += -(xyzt(i1,k)-xyzt(i2,k))/pow(covparms(1),3)/dd * xyzt(i1,k);
                    pth += -(xyzt(i2,k)-xyzt(i1,k))/pow(covparms(1),3)/dd * xyzt(i2,k);
                }
                dcovmat(i2,i1,nisoparms+j) = (dcovmat(i2,i1,1)+dcovmat(i2,i1,2))/pth;
                // pd wrt to location component * pd wrt to basis
                double pd = 0.0;
                for(int k=0; k<3; k++){ 
                    pd += (xyzt(i2,k)-xyzt(i1,k))/pow(covparms(1),2)/dd * grad_basis(i2,j,k);
                    pd += (xyzt(i1,k)-xyzt(i2,k))/pow(covparms(1),2)/dd * grad_basis(i1,j,k);
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
