#ifndef ARMAONEPASS_H
#define ARMAONEPASS_H

#include <assert.h>
#include <RcppArmadillo.h>
#include "covmatrix_funs.h"
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

 
void arma_onepass_compute_pieces(
    NumericVector covparms, 
    StringVector covfun_name,
    const NumericMatrix locs, 
    IntegerMatrix NNarray,
    NumericVector& y, 
    NumericMatrix X,
    mat* XSX,
    vec* ySX,
    double* ySy,
    double* logdet,
    cube* dXSX,
    mat* dySX,
    vec* dySy,
    vec* dlogdet,
    mat* ainfo
){

    // data dimensions
    int n = y.length();
    int m = NNarray.ncol();
    int p = X.ncol();
    int nparms = covparms.length();
    int dim = locs.ncol();
    
    // declare covariance fun and derivative
    mat (*p_covfun)(NumericVector, NumericMatrix);
    cube (*p_d_covfun)(NumericVector, NumericMatrix);
    
    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];

    // would prefer to have this in a separate function
    // that function is not working right now
    if( covfun_name_string.compare("arma_matern_isotropic") == 0 )
    { 
        p_covfun = &arma_matern_isotropic; 
        p_d_covfun = &d_arma_matern_isotropic;
    } 
    else if( covfun_name_string.compare("arma_matern_anisotropic2D") == 0 )
    { 
        p_covfun = &arma_matern_anisotropic2D; 
        p_d_covfun = &d_arma_matern_anisotropic2D;
    } 
    else if( covfun_name_string.compare("arma_exponential_anisotropic3D") == 0 )
    { 
        p_covfun = &arma_exponential_anisotropic3D; 
        p_d_covfun = &d_arma_exponential_anisotropic3D;
    } 
    else if( covfun_name_string.compare("arma_matern_anisotropic3D") == 0 )
    { 
        p_covfun = &arma_matern_anisotropic3D; 
        p_d_covfun = &d_arma_matern_anisotropic3D;
    } 
    else if( covfun_name_string.compare("arma_matern_nonstat_var") == 0 )
    { 
        p_covfun = &arma_matern_nonstat_var; 
        p_d_covfun = &d_arma_matern_nonstat_var;
    } 
    else if( covfun_name_string.compare("arma_exponential_isotropic") == 0 )
    { 
        p_covfun = &arma_exponential_isotropic; 
        p_d_covfun = &d_arma_exponential_isotropic;
    }
    else if( covfun_name_string.compare("arma_matern_sphere") == 0 )
    { 
        p_covfun = &arma_matern_sphere; 
        p_d_covfun = &d_arma_matern_sphere;
    }
    else if( covfun_name_string.compare("arma_matern_sphere_warp") == 0 )
    { 
        p_covfun = &arma_matern_sphere_warp; 
        p_d_covfun = &d_arma_matern_sphere_warp;
    }
    else if( covfun_name_string.compare("arma_matern_spheretime_warp") == 0 )
    { 
        p_covfun = &arma_matern_spheretime_warp; 
        p_d_covfun = &d_arma_matern_spheretime_warp;
    }
    else if( covfun_name_string.compare("arma_matern_spheretime") == 0 )
    { 
        p_covfun = &arma_matern_spheretime; 
        p_d_covfun = &d_arma_matern_spheretime;
    }
    else if( covfun_name_string.compare("arma_matern_spheretime_nonstatvar") == 0 )
    { 
        p_covfun = &arma_matern_spheretime_nonstatvar; 
        p_d_covfun = &d_arma_matern_spheretime_nonstatvar;
    }
    else if( covfun_name_string.compare("arma_matern_spacetime") == 0 )
    { 
        p_covfun = &arma_matern_spacetime; 
        p_d_covfun = &d_arma_matern_spacetime;
    }
    else if( covfun_name_string.compare("arma_matern_scaledim") == 0 )
    { 
        p_covfun = &arma_matern_scaledim;
        p_d_covfun = &d_arma_matern_scaledim;
    }
    else { // stop the program
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
        assert(0);
    }

    // loop over every observation    
    for(int i=0; i<n; i++){
    
        Rcpp::checkUserInterrupt();
        int bsize = std::min(i+1,m);

        // first, fill in ysub, locsub, and X0 in reverse order
        NumericMatrix locsub(bsize, dim);
        arma::vec ysub(bsize);
        arma::mat X0( bsize, p );
        for(int j=bsize-1; j>=0; j--){
            ysub(bsize-1-j) = y( NNarray(i,j)-1 );
            for(int k=0;k<dim;k++){ locsub(bsize-1-j,k) = locs( NNarray(i,j)-1, k ); }
            for(int k=0;k<p;k++){ X0(bsize-1-j,k) = X( NNarray(i,j)-1, k ); }
        }
        
        // compute covariance matrix and derivatives and take cholesky
        arma::mat covmat = (*p_covfun)( covparms, locsub );
        arma::cube dcovmat = (*p_d_covfun)( covparms, locsub );
        arma::mat cholmat = eye( size(covmat) );
        chol( cholmat, covmat, "lower" );
        
        // i1 is conditioning set, i2 is response        
        arma::span i1 = span(0,bsize-2);
        arma::span i2 = span(bsize-1,bsize-1);
        
        // get last row of cholmat
        arma::vec onevec = zeros(bsize);
        onevec(bsize-1) = 1.0;
        arma::vec choli2 = solve( trimatu(cholmat.t()), onevec );
        
        bool cond = bsize > 1;
        double fac = 1.0;
        
        // do solves with X and y
        arma::mat LiX0 = solve( trimatl(cholmat), X0 );
        arma::vec Liy0 = solve( trimatl(cholmat), ysub );
        
        // loglik objects
        *logdet += 2.0*std::log( as_scalar(cholmat(i2,i2)) ); 
        *XSX +=    LiX0.rows(i2).t() * LiX0.rows(i2);
        *ySy +=    pow( as_scalar(Liy0(i2)), 2 );
        *ySX += ( Liy0(i2) * LiX0.rows(i2) ).t();
        
        // gradient objects
        // LidSLi3 is last column of Li * (dS_j) * Lit for 1 parameter i
        // LidSLi2 stores these columns in a matrix for all parameters
        arma::mat LidSLi2(bsize,nparms);
        
        if(cond){ // if we condition on anything
            
            for(int j=0; j<nparms; j++){
                // compute last column of Li * (dS_j) * Lit
                arma::vec LidSLi3 = solve( trimatl(cholmat), dcovmat.slice(j) * choli2 );
                // store LiX0.t() * LidSLi3 and Liy0.t() * LidSLi3
                arma::vec v1 = LiX0.t() * LidSLi3;
                double s1 = as_scalar( Liy0.t() * LidSLi3 ); 
                // update all quantities
                // bottom-right corner gets double counted, so need to subtract it off
                (*dXSX).slice(j) += v1 * LiX0.rows(i2) + ( v1 * LiX0.rows(i2) ).t() - 
                    as_scalar(LidSLi3(i2)) * ( LiX0.rows(i2).t() * LiX0.rows(i2) );
                (*dySy)(j) += as_scalar( 2.0 * s1 * Liy0(i2)  - 
                    LidSLi3(i2) * Liy0(i2) * Liy0(i2) );
                (*dySX).col(j) += (  s1 * LiX0.rows(i2) + ( v1 * Liy0(i2) ).t() -  
                    as_scalar( LidSLi3(i2) ) * LiX0.rows(i2) * as_scalar( Liy0(i2))).t();
                (*dlogdet)(j) += as_scalar( LidSLi3(i2) );
                // store last column of Li * (dS_j) * Lit
                LidSLi2.col(j) = LidSLi3;
            }

            // fisher information object
            // bottom right corner gets double counted, so subtract it off
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (*ainfo)(i,j) += 
                    1.0*accu( LidSLi2.col(i) % LidSLi2.col(j) ) - 
                    0.5*accu( LidSLi2.rows(i2).col(j) %
                              LidSLi2.rows(i2).col(i) );
            }}
            
        } else { // similar calculations, but for when there is no conditioning set
            for(int j=0; j<nparms; j++){
                arma::mat LidSLi = solve( trimatl(cholmat), dcovmat.slice(j) );
                LidSLi = solve( trimatl(cholmat), LidSLi.t() );
                (*dXSX).slice(j) += LiX0.t() *  LidSLi * LiX0; 
                (*dySy)(j) += as_scalar( Liy0.t() * LidSLi * Liy0 );
                (*dySX).col(j) += ( ( Liy0.t() * LidSLi ) * LiX0 ).t();
                (*dlogdet)(j) += trace( LidSLi );
                LidSLi2.col(j) = LidSLi;
            }
            
            // fisher information object
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (*ainfo)(i,j) += 0.5*accu( LidSLi2.col(i) % LidSLi2.col(j) ); 
            }}

        }        

    }
}    

    

    
    
    
    
void arma_onepass_compute_pieces_grouped(
    NumericVector covparms, 
    StringVector covfun_name,
    const NumericMatrix locs, 
    List NNlist,
    NumericVector& y, 
    NumericMatrix X,
    mat* XSX,
    vec* ySX,
    double* ySy,
    double* logdet,
    cube* dXSX,
    mat* dySX,
    vec* dySy,
    vec* dlogdet,
    mat* ainfo
){
    

    // data dimensions
    int n = y.length();
    //int m = NNarray.ncol();
    int p = X.ncol();
    int nparms = covparms.length();
    int dim = locs.ncol();
    
    // declare covariance fun and derivative
    mat (*p_covfun)(NumericVector, NumericMatrix);
    cube (*p_d_covfun)(NumericVector, NumericMatrix);

    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];

    // would prefer to have this in a separate function
    // that function is not working right now
    if( covfun_name_string.compare("arma_matern_isotropic") == 0 )
    { 
        p_covfun = &arma_matern_isotropic; 
        p_d_covfun = &d_arma_matern_isotropic;
    } 
    else if( covfun_name_string.compare("arma_matern_anisotropic2D") == 0 )
    { 
        p_covfun = &arma_matern_anisotropic2D; 
        p_d_covfun = &d_arma_matern_anisotropic2D;
    } 
    else if( covfun_name_string.compare("arma_exponential_anisotropic3D") == 0 )
    { 
        p_covfun = &arma_exponential_anisotropic3D; 
        p_d_covfun = &d_arma_exponential_anisotropic3D;
    } 
    else if( covfun_name_string.compare("arma_matern_anisotropic3D") == 0 )
    { 
        p_covfun = &arma_matern_anisotropic3D; 
        p_d_covfun = &d_arma_matern_anisotropic3D;
    } 
    else if( covfun_name_string.compare("arma_matern_nonstat_var") == 0 )
    { 
        p_covfun = &arma_matern_nonstat_var; 
        p_d_covfun = &d_arma_matern_nonstat_var;
    } 
    else if( covfun_name_string.compare("arma_exponential_isotropic") == 0 )
    { 
        p_covfun = &arma_exponential_isotropic; 
        p_d_covfun = &d_arma_exponential_isotropic;
    }
    else if( covfun_name_string.compare("arma_matern_sphere") == 0 )
    { 
        p_covfun = &arma_matern_sphere; 
        p_d_covfun = &d_arma_matern_sphere;
    }
    else if( covfun_name_string.compare("arma_matern_sphere_warp") == 0 )
    { 
        p_covfun = &arma_matern_sphere_warp; 
        p_d_covfun = &d_arma_matern_sphere_warp;
    }
    else if( covfun_name_string.compare("arma_matern_spheretime_warp") == 0 )
    { 
        p_covfun = &arma_matern_spheretime_warp; 
        p_d_covfun = &d_arma_matern_spheretime_warp;
    }
    else if( covfun_name_string.compare("arma_matern_spheretime") == 0 )
    { 
        p_covfun = &arma_matern_spheretime; 
        p_d_covfun = &d_arma_matern_spheretime;
    }
    else if( covfun_name_string.compare("arma_matern_spheretime_nonstatvar") == 0 )
    { 
        p_covfun = &arma_matern_spheretime_nonstatvar; 
        p_d_covfun = &d_arma_matern_spheretime_nonstatvar;
    }
    else if( covfun_name_string.compare("arma_matern_spacetime") == 0 )
    { 
        p_covfun = &arma_matern_spacetime; 
        p_d_covfun = &d_arma_matern_spacetime;
    }
    else if( covfun_name_string.compare("arma_matern_scaledim") == 0 )
    { 
        p_covfun = &arma_matern_scaledim;
        p_d_covfun = &d_arma_matern_scaledim;
    }
    else { // stop the program
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
        assert(0);
    }

    // vector of all indices
    IntegerVector all_inds = NNlist["all_inds"];
    // vector of local response indices
    IntegerVector local_resp_inds = as<IntegerVector>(NNlist["local_resp_inds"]);
    // vector of global response indices
    IntegerVector global_resp_inds = as<IntegerVector>(NNlist["global_resp_inds"]);
    // last index of each block in all_inds
    IntegerVector last_ind_of_block = as<IntegerVector>(NNlist["last_ind_of_block"]);
    // last response index of each block in local_resp_inds and global_resp_inds
    IntegerVector last_resp_of_block = as<IntegerVector>(NNlist["last_resp_of_block"]);

    int nb = last_ind_of_block.size();  // number of blocks


    // loop over every block
    for(int i=0; i<nb; i++){

        Rcpp::checkUserInterrupt();
        
        // first ind and last ind are the positions in all_inds
        // of the observations for block i.
        // these come in 1-indexing and are converted to 0-indexing here
        int first_ind;
        int last_ind;
        if(i==0){ first_ind = 0; } else {first_ind = last_ind_of_block[i-1]+1-1; }
        last_ind = last_ind_of_block[i]-1;
        int bsize = last_ind - first_ind + 1;
        
        int first_resp;
        int last_resp;
        if(i==0){ first_resp = 0; } else {first_resp = last_resp_of_block[i-1]+1-1; }
        last_resp = last_resp_of_block[i]-1;
        int rsize = last_resp - first_resp + 1;
        
        arma::uvec whichresp(rsize);
        for(int j=0; j<rsize; j++){
            whichresp(j) = local_resp_inds(first_resp+j) - 1;
        }

        // fill in ysub, locsub, and X0 in forward order
        NumericMatrix locsub(bsize, dim);
        arma::vec ysub(bsize);
        arma::mat X0( bsize, p );
        for(int j=0; j<bsize; j++){
            int jglobal = all_inds[first_ind + j] - 1;
            ysub(j) = y( jglobal );
            for(int k=0;k<dim;k++){ locsub(j,k) = locs( jglobal, k ); }
            for(int k=0;k<p;k++){ X0(j,k) = X( jglobal, k ); }
        }

        // compute covariance matrix and derivatives and take cholesky
        arma::mat covmat = (*p_covfun)( covparms, locsub );
        arma::cube dcovmat = (*p_d_covfun)( covparms, locsub );
        arma::mat cholmat = eye( size(covmat) );
        chol( cholmat, covmat, "lower" );

        // get response rows of inverse cholmat, put in column vectors
        arma::mat onemat = zeros(bsize,rsize);
        for(int j=0; j<rsize; j++){ 
            onemat(whichresp(j),j) = 1.0;
        }
        arma::mat choli2 = solve( trimatu(cholmat.t()), onemat );
        bool cond = bsize > rsize;
        double fac = 1.0;
        
        // do solves with X and y
        arma::mat LiX0 = solve( trimatl(cholmat), X0 );
        arma::vec Liy0 = solve( trimatl(cholmat), ysub );
        // loglik objects
        for(int j=0; j<rsize; j++){
            int ii = whichresp(j);
            *logdet += 2.0*std::log( as_scalar(cholmat(ii,ii)) ); 
            *XSX +=    LiX0.row(ii).t() * LiX0.row(ii);
            *ySy +=    pow( as_scalar(Liy0(ii)), 2 );
            *ySX += ( Liy0(ii) * LiX0.row(ii) ).t();
        }    
        
        // gradient objects
        // LidSLi3 is last column of Li * (dS_j) * Lit for 1 parameter i
        // LidSLi2 stores these columns in a matrix for all parameters
        if(cond){ // if we condition on anything
            
            arma::cube LidSLi2 = arma::cube(bsize,rsize,nparms,fill::zeros);
            for(int j=0; j<nparms; j++){
                // compute last column of Li * (dS_j) * Lit
                arma::mat LidSLi4 = 
                    solve( trimatl(cholmat), 
                           dcovmat.slice(j) * choli2 );
                
                for(int k=0; k<rsize; k++){
                    int i2 = whichresp(k);
                    arma::span i1 = span(0,i2);
                    arma::vec LidSLi3 = LidSLi4(i1,k); 
                    // store LiX0.t() * LidSLi3 and Liy0.t() * LidSLi3
                    arma::vec v1 = LiX0.rows(i1).t() * LidSLi3;
                    double s1 = dot( Liy0(i1), LidSLi3 ); 
                    // update all quantities
                    // bottom-right corner gets double counted, so need to subtract it off
                    (*dXSX).slice(j) += v1 * LiX0.row(i2) + ( v1 * LiX0.row(i2) ).t() - 
                        LidSLi3(i2) * LiX0.row(i2).t() * LiX0.row(i2);
                    (*dySy)(j) += 2.0 * s1 * Liy0(i2)  - 
                        LidSLi3(i2) * Liy0(i2) * Liy0(i2);
                    (*dySX).col(j) += (  s1 * LiX0.row(i2) + ( v1 * Liy0(i2) ).t() -  
                        LidSLi3(i2) * LiX0.row(i2) * Liy0(i2) ).t();
                    (*dlogdet)(j) += LidSLi3(i2);
                    // store last column of Li * (dS_j) * Lit
                    LidSLi2.subcube(i1, span(k,k), span(j,j)) = LidSLi3;
                }
            }
            // fisher information object
            // bottom right corner gets double counted, so subtract it off
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (*ainfo)(i,j) += accu( LidSLi2.slice(i) % LidSLi2.slice(j) );
                for(int k=0; k<rsize; k++){
                    int i2 = whichresp(k);
                    (*ainfo)(i,j) -= 0.5*LidSLi2(i2,k,j) * LidSLi2(i2,k,i);
                }
            }}
        } else { // similar calculations, but for when there is no conditioning set
            arma::cube LidSLi2(bsize,bsize,nparms);
            for(int j=0; j<nparms; j++){
                arma::mat LidSLi = solve( trimatl(cholmat), dcovmat.slice(j) );
                LidSLi = solve( trimatl(cholmat), LidSLi.t() );
                (*dXSX).slice(j) += LiX0.t() *  LidSLi * LiX0; 
                (*dySy)(j) += as_scalar( Liy0.t() * LidSLi * Liy0 );
                (*dySX).col(j) += ( ( Liy0.t() * LidSLi ) * LiX0 ).t();
                (*dlogdet)(j) += trace( LidSLi );
                LidSLi2.slice(j) = LidSLi;
            }
            
            // fisher information object
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (*ainfo)(i,j) += 0.5*accu( LidSLi2.slice(i) % LidSLi2.slice(j) ); 
            }}
        }
    }
}    

    
    

void arma_onepass_synthesize(
    NumericVector covparms, 
    StringVector covfun_name,
    const NumericMatrix locs, 
    IntegerMatrix NNarray,
    NumericVector& y, 
    NumericMatrix X,
    NumericVector* ll, 
    NumericVector* betahat,
    NumericVector* grad,
    NumericMatrix* info,
    NumericMatrix* betainfo ){

    // data dimensions
    int n = y.length();
    int m = NNarray.ncol();
    int p = X.ncol();
    int nparms = covparms.length();
    int dim = locs.ncol();
    
    // likelihood objects
    arma::mat XSX = arma::mat(p, p, fill::zeros);
    arma::vec ySX = arma::vec(p, fill::zeros);
    double ySy = 0.0;
    double logdet = 0.0;
    
    // gradient objects    
    arma::cube dXSX = arma::cube(p,p,nparms,fill::zeros);
    arma::mat dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec dySy = arma::vec(nparms, fill::zeros);
    arma::vec dlogdet = arma::vec(nparms, fill::zeros);
    // fisher information
    arma::mat ainfo = arma::mat(nparms, nparms, fill::zeros);

    // this is where the big computation happens
    arma_onepass_compute_pieces(
        covparms, covfun_name, locs, NNarray, y, X,
        &XSX, &ySX, &ySy, &logdet, &dXSX, &dySX, &dySy, &dlogdet, &ainfo
    );
        
    // synthesize everything and update loglik, grad, beta, betainfo, info
    
    // betahat and dbeta
    arma::vec abeta = solve( XSX, ySX );
    for(int j=0; j<p; j++){ (*betahat)(j) = abeta(j); };

    arma::mat dbeta(p,nparms);
    for(int j=0; j<nparms; j++){
        dbeta.col(j) = solve( XSX, dySX.col(j) - dXSX.slice(j) * abeta );
    }
    // get sigmahatsq
    double sig2 = ( ySy - 2.0*as_scalar( ySX.t() * abeta ) + 
        as_scalar( abeta.t() * XSX * abeta ) )/n;
    // loglikelihood
    (*ll)(0) = -0.5*( n*std::log(2.0*M_PI) + logdet + n*sig2 ); 
    // gradient
    for(int j=0; j<nparms; j++){
        (*grad)(j) = 0.0;
        (*grad)(j) -= 0.5*dlogdet(j);
        (*grad)(j) += 0.5*dySy(j);
        (*grad)(j) -= 1.0*as_scalar( abeta.t() * dySX.col(j) );
        (*grad)(j) += 1.0*as_scalar( ySX.t() * dbeta.col(j) );
        (*grad)(j) += 0.5*as_scalar( abeta.t() * dXSX.slice(j) * abeta );
        (*grad)(j) -= 1.0*as_scalar( abeta.t() * XSX * dbeta.col(j) );
    }
    // betainfo
    for(int i=0; i<p; i++){ for(int j=0; j<i+1; j++){
        (*betainfo)(i,j) = XSX(i,j);
        (*betainfo)(j,i) = XSX(j,i);
    }}
    // fisher information
    for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
        (*info)(i,j) = ainfo(i,j);
        (*info)(j,i) = (*info)(i,j);
    }}

}


void arma_onepass_synthesize_grouped(
    NumericVector covparms, 
    StringVector covfun_name,
    const NumericMatrix locs, 
    List NNlist,
    NumericVector& y, 
    NumericMatrix X,
    NumericVector* ll, 
    NumericVector* betahat,
    NumericVector* grad,
    NumericMatrix* info,
    NumericMatrix* betainfo ){

    // data dimensions
    int n = y.length();
    //int m = NNarray.ncol();
    int p = X.ncol();
    int nparms = covparms.length();
    int dim = locs.ncol();
    
    // likelihood objects
    arma::mat XSX = arma::mat(p, p, fill::zeros);
    arma::vec ySX = arma::vec(p, fill::zeros);
    double ySy = 0.0;
    double logdet = 0.0;
    
    // gradient objects    
    arma::cube dXSX = arma::cube(p,p,nparms,fill::zeros);
    arma::mat dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec dySy = arma::vec(nparms, fill::zeros);
    arma::vec dlogdet = arma::vec(nparms, fill::zeros);
    // fisher information
    arma::mat ainfo = arma::mat(nparms, nparms, fill::zeros);
    // this is where the big computation happens
    arma_onepass_compute_pieces_grouped(
        covparms, covfun_name, locs, NNlist, y, X,
        &XSX, &ySX, &ySy, &logdet, &dXSX, &dySX, &dySy, &dlogdet, &ainfo
    );
    
        
        
    // synthesize everything and update loglik, grad, beta, betainfo, info
    
    // betahat and dbeta
    arma::vec abeta = solve( XSX, ySX );
    for(int j=0; j<p; j++){ (*betahat)(j) = abeta(j); };
    
    arma::mat dbeta(p,nparms);
    for(int j=0; j<nparms; j++){
        dbeta.col(j) = solve( XSX, dySX.col(j) - dXSX.slice(j) * abeta );
    }
    // get sigmahatsq
    double sig2 = ( ySy - 2.0*as_scalar( ySX.t() * abeta ) + 
        as_scalar( abeta.t() * XSX * abeta ) )/n;
    // loglikelihood
    (*ll)(0) = -0.5*( n*std::log(2.0*M_PI) + logdet + n*sig2 ); 
    // gradient
    for(int j=0; j<nparms; j++){
        (*grad)(j) = 0.0;
        (*grad)(j) -= 0.5*dlogdet(j);
        (*grad)(j) += 0.5*dySy(j);
        (*grad)(j) -= 1.0*as_scalar( abeta.t() * dySX.col(j) );
        (*grad)(j) += 1.0*as_scalar( ySX.t() * dbeta.col(j) );
        (*grad)(j) += 0.5*as_scalar( abeta.t() * dXSX.slice(j) * abeta );
        (*grad)(j) -= 1.0*as_scalar( abeta.t() * XSX * dbeta.col(j) );
    }
    // betainfo
    for(int i=0; i<p; i++){ for(int j=0; j<i+1; j++){
        (*betainfo)(i,j) = XSX(i,j);
        (*betainfo)(j,i) = XSX(j,i);
    }}
    // fisher information
    for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
        (*info)(i,j) = ainfo(i,j);
        (*info)(j,i) = (*info)(i,j);
    }}

}



/*

void arma_onepass_profile_synthesize(
    NumericVector subparms, 
    double* sigmasq,
    StringVector covfun_name,
    const NumericMatrix locs, 
    IntegerMatrix NNarray,
    NumericVector& y, 
    NumericMatrix X,
    NumericVector* ll, 
    NumericVector* betahat,
    NumericVector* grad,
    NumericMatrix* info,
    NumericMatrix* betainfo ){
    
    // define covparms, fill in 1 for sigmasq, rest with subparms
    NumericVector covparms(subparms.length()+1);
    covparms(0) = 1.0;
    for(int i=0; i<subparms.length(); i++){
        covparms(i+1) = subparms(i);
    }
    
    // data dimensions
    int n = y.length();
    int m = NNarray.ncol();
    int p = X.ncol();
    int nparms = covparms.length();
    int dim = locs.ncol();
    
    // likelihood objects
    arma::mat XSX = arma::mat(p, p, fill::zeros);
    arma::vec ySX = arma::vec(p, fill::zeros);
    double ySy = 0.0;
    double logdet = 0.0;
    
    // gradient objects    
    arma::cube dXSX = arma::cube(p,p,nparms,fill::zeros);
    arma::mat dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec dySy = arma::vec(nparms, fill::zeros);
    arma::vec dlogdet = arma::vec(nparms, fill::zeros);
    // fisher information
    arma::mat ainfo = arma::mat(nparms, nparms, fill::zeros);
    
    // this is where the big computation happens
    arma_onepass_compute_pieces(
        covparms, covfun_name, locs, NNarray, y, X,
        &XSX, &ySX, &ySy, &logdet, &dXSX, &dySX, &dySy, &dlogdet, &ainfo
    );
        
        
    // synthesize everything and update loglik, grad, beta, betainfo, info
    
    // betahat and dbeta
    arma::vec abeta = solve( XSX, ySX );
    for(int j=0; j<p; j++){ (*betahat)(j) = abeta(j); };

    arma::mat dbeta(p,nparms);
    for(int j=1; j<nparms; j++){
        dbeta.col(j-1) = solve( XSX, dySX.col(j) - dXSX.slice(j) * abeta );
    }
    
    // get sigmasq
    *sigmasq = ( ySy - 2.0*as_scalar( ySX.t() * abeta ) + 
        as_scalar( abeta.t() * XSX * abeta ) )/n;
    
    // loglikelihood
    (*ll)(0) = -0.5*( n*std::log(2.0*M_PI) + logdet + n*log(*sigmasq) + n ); 
    
    // gradient
    for(int j=1; j<nparms; j++){
        (*grad)(j-1) = 0.0;
        (*grad)(j-1) -= 0.5*dlogdet(j);
        (*grad)(j-1) += 0.5*dySy(j)/(*sigmasq);
        (*grad)(j-1) -= 1.0*as_scalar( abeta.t() * dySX.col(j) )/(*sigmasq);
        (*grad)(j-1) += 1.0*as_scalar( ySX.t() * dbeta.col(j) )/(*sigmasq);
        (*grad)(j-1) += 0.5*as_scalar( abeta.t() * dXSX.slice(j) * abeta )/(*sigmasq);
        (*grad)(j-1) -= 1.0*as_scalar( abeta.t() * XSX * dbeta.col(j) )/(*sigmasq);
    }
    // betainfo
    for(int i=0; i<p; i++){ for(int j=0; j<i+1; j++){
        (*betainfo)(i,j) = XSX(i,j)/(*sigmasq);
        (*betainfo)(j,i) = XSX(j,i)/(*sigmasq);
    }}
    // fisher information
    for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
        (*info)(i,j) = ainfo(i,j);
        (*info)(j,i) = (*info)(i,j);
    }}

}

*/

// [[Rcpp::export]]
NumericMatrix vecchia_Linv(
    NumericVector covparms,
    StringVector covfun_name,
    const NumericMatrix locs,
    IntegerMatrix NNarray ){
    
    // data dimensions
    int n = locs.nrow();
    int m = NNarray.ncol();
    int nparms = covparms.length();
    int dim = locs.ncol();
    NumericMatrix Linv(n,m);
    
    // declare covariance fun and derivative
    mat (*p_covfun)(NumericVector, NumericMatrix);
    cube (*p_d_covfun)(NumericVector, NumericMatrix);
    
    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];

    // would prefer to have this in a separate function
    // that function is not working right now
    if( covfun_name_string.compare("arma_matern_isotropic") == 0 )
    { 
        p_covfun = &arma_matern_isotropic; 
        p_d_covfun = &d_arma_matern_isotropic;
    } 
    else if( covfun_name_string.compare("arma_matern_anisotropic2D") == 0 )
    { 
        p_covfun = &arma_matern_anisotropic2D; 
        p_d_covfun = &d_arma_matern_anisotropic2D;
    } 
    else if( covfun_name_string.compare("arma_exponential_anisotropic3D") == 0 )
    { 
        p_covfun = &arma_exponential_anisotropic3D; 
        p_d_covfun = &d_arma_exponential_anisotropic3D;
    } 
    else if( covfun_name_string.compare("arma_matern_anisotropic3D") == 0 )
    { 
        p_covfun = &arma_matern_anisotropic3D; 
        p_d_covfun = &d_arma_matern_anisotropic3D;
    } 
    else if( covfun_name_string.compare("arma_matern_nonstat_var") == 0 )
    { 
        p_covfun = &arma_matern_nonstat_var; 
        p_d_covfun = &d_arma_matern_nonstat_var;
    } 
    else if( covfun_name_string.compare("arma_exponential_isotropic") == 0 )
    { 
        p_covfun = &arma_exponential_isotropic; 
        p_d_covfun = &d_arma_exponential_isotropic;
    }
    else if( covfun_name_string.compare("arma_matern_sphere") == 0 )
    { 
        p_covfun = &arma_matern_sphere; 
        p_d_covfun = &d_arma_matern_sphere;
    }
    else if( covfun_name_string.compare("arma_matern_sphere_warp") == 0 )
    { 
        p_covfun = &arma_matern_sphere_warp; 
        p_d_covfun = &d_arma_matern_sphere_warp;
    }
    else if( covfun_name_string.compare("arma_matern_spheretime_warp") == 0 )
    { 
        p_covfun = &arma_matern_spheretime_warp; 
        p_d_covfun = &d_arma_matern_spheretime_warp;
    }
    else if( covfun_name_string.compare("arma_matern_spheretime") == 0 )
    { 
        p_covfun = &arma_matern_spheretime; 
        p_d_covfun = &d_arma_matern_spheretime;
    }
    else if( covfun_name_string.compare("arma_matern_spheretime_nonstatvar") == 0 )
    { 
        p_covfun = &arma_matern_spheretime_nonstatvar; 
        p_d_covfun = &d_arma_matern_spheretime_nonstatvar;
    }
    else if( covfun_name_string.compare("arma_matern_spacetime") == 0 )
    { 
        p_covfun = &arma_matern_spacetime; 
        p_d_covfun = &d_arma_matern_spacetime;
    }
    else if( covfun_name_string.compare("arma_matern_scaledim") == 0 )
    { 
        p_covfun = &arma_matern_scaledim;
        p_d_covfun = &d_arma_matern_scaledim;
    }
    else { // stop the program
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
        assert(0);
    }

    // loop over every observation    
    for(int i=0; i<n; i++){
    
        Rcpp::checkUserInterrupt();
        int bsize = std::min(i+1,m);

        // first, fill in ysub, locsub, and X0 in reverse order
        NumericMatrix locsub(bsize, dim);
        for(int j=bsize-1; j>=0; j--){
            for(int k=0;k<dim;k++){ locsub(bsize-1-j,k) = locs( NNarray(i,j)-1, k ); }
        }
        
        // compute covariance matrix and derivatives and take cholesky
        arma::mat covmat = (*p_covfun)( covparms, locsub );
        arma::mat cholmat = eye( size(covmat) );
        bool checkchol = chol( cholmat, covmat, "lower" );
        //if(i==n-1){ Rcout << checkchol << endl; }
        // get last row of cholmat
        arma::vec onevec = zeros(bsize);
        onevec(bsize-1) = 1.0;
        arma::vec choli2;
        if( checkchol == 0 ){
            choli2 = onevec;
            //Rcout << "failed chol" << endl;
            //Rcout << checkchol << endl;
        } else {
            choli2 = solve( trimatu(cholmat.t()), onevec );
        }
        for(int j=bsize-1; j>=0; j--){
            Linv(i,bsize-1-j) = choli2(j);
        }
    }    
    return Linv;    
}


#endif
    
