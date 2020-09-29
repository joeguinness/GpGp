#ifndef ARMA_BASIS_H
#define ARMA_BASIS_H

// basis functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]

//' compute gradient of spherical harmonics functions
//' 
//' @param xyz xyz coordinates of locations on sphere
//' @param Lmax largest degree of spherical harmonics. 
//' Current only Lmax=2 supported
// [[Rcpp::export]]
arma::cube sph_grad_xyz( arma::mat xyz, int Lmax ){
    
    int nbasis = (Lmax+1)*(Lmax+1) - 4;
    int n = xyz.n_rows;
    cube grad_basis = cube( n, nbasis, 3, fill::zeros );
    
    if(Lmax > 1){
        double c2 = pow( 15.0/ M_PI, 0.5 );
        double c3 = pow(  5.0/ M_PI, 0.5 );
        for(int i=0; i<n; i++){
            // first basis
            grad_basis(i,0,0) = 0.5*c2*xyz( i , 1 );
            grad_basis(i,0,1) = 0.5*c2*xyz( i , 0 );
            // second basis
            grad_basis(i,1,1) = 0.5*c2*xyz( i , 2 );
            grad_basis(i,1,2) = 0.5*c2*xyz( i , 1 );
            // third basis
            grad_basis(i,2,0) = 0.25*c3*( -2.0*xyz( i , 0 ) );
            grad_basis(i,2,1) = 0.25*c3*( -2.0*xyz( i , 1 ) );
            grad_basis(i,2,2) = 0.25*c3*( +4.0*xyz( i , 2 ) );
            // fourth basis
            grad_basis(i,3,0) = 0.5*c2*xyz( i , 2 );
            grad_basis(i,3,2) = 0.5*c2*xyz( i , 0 );
            // fifth basis
            grad_basis(i,4,0) = 0.25*c2*( +2.0*xyz( i , 0 ) );
            grad_basis(i,4,1) = 0.25*c2*( -2.0*xyz( i , 1 ) );
        }
    }

    return grad_basis;
    
}

#endif
