#ifndef ONEPASS_H
#define ONEPASS_H

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <RcppArmadillo.h>
#include "covmatrix_funs.h"

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

 
void compute_pieces(
    arma::vec covparms, 
    StringVector covfun_name,
    arma::mat locs, 
    arma::Mat<uword> NNarray,
    arma::vec y, 
    arma::mat X,
    arma::mat* XSX,
    arma::vec* ySX,
    double* ySy,
    double* logdet,
    arma::cube* dXSX,
    arma::mat* dySX,
    arma::vec* dySy,
    arma::vec* dlogdet,
    arma::mat* ainfo,
    int profbeta,
    int grad_info
){

    // data dimensions
    uword n = y.n_elem;
    uword m = NNarray.n_cols;
    uword p = X.n_cols;
    uword nparms = covparms.n_elem;
    uword dim = locs.n_cols;
    

    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];
    
    // assign covariance fun and derivative based on covfun_name_string
    covfun_t covstruct = get_covfun(covfun_name_string);
    mat (*p_covfun)(arma::vec, arma::mat) = covstruct.p_covfun;
    cube (*p_d_covfun)(arma::vec, arma::mat) = covstruct.p_d_covfun;


   /* arma::mat l_XSX = arma::mat(p, p, fill::zeros);
    arma::vec l_ySX = arma::vec(p, fill::zeros);
    double l_ySy = 0.0;
    double l_logdet = 0.0;
    arma::cube l_dXSX = arma::cube(p, p, nparms,fill::zeros);
    arma::mat l_dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec l_dySy = arma::vec(nparms, fill::zeros);
    arma::vec l_dlogdet = arma::vec(nparms, fill::zeros);
    arma::mat l_ainfo = arma::mat(nparms, nparms, fill::zeros);*/

#pragma omp parallel shared(covparms,NNarray,locs,y,X,n,m,p,nparms,dim,p_covfun,p_d_covfun)
{
    arma::mat l_XSX = arma::mat(p, p, fill::zeros);
    arma::vec l_ySX = arma::vec(p, fill::zeros);
    double l_ySy = 0.0;
    double l_logdet = 0.0;
    arma::cube l_dXSX = arma::cube(p, p, nparms,fill::zeros);
    arma::mat l_dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec l_dySy = arma::vec(nparms, fill::zeros);
    arma::vec l_dlogdet = arma::vec(nparms, fill::zeros);
    arma::mat l_ainfo = arma::mat(nparms, nparms, fill::zeros);
    
#pragma omp for
    for(uword i=0;i<n;i++){
     
	uword bsize = std::min(i+1,m);

        // first, fill in ysub, locsub, and X0 in reverse order
        arma::mat locsub(bsize, dim);
        arma::vec ysub(bsize);
        arma::mat X0(bsize, p);
        
	for(uword j=0;j<bsize;j++){
            ysub(j) = y(NNarray(i,bsize-1-j)-1);
            for(uword k=0;k<dim;k++){ locsub(j,k) = locs(NNarray(i,bsize-1-j)-1, k); }
            if(profbeta){ 
                for(uword k=0;k<p;k++){ X0(j,k) = X(NNarray(i,bsize-1-j)-1, k); } 
            }
        }
        
        // compute covariance matrix and derivatives and take cholesky
        arma::mat covmat = (*p_covfun)(covparms, locsub);
        arma::cube dcovmat;
        if(grad_info){ 
            dcovmat = (*p_d_covfun)( covparms, locsub ); 
        }
		
        arma::mat cholmat = eye( size(covmat) );
        chol( cholmat, covmat, "lower" );
        
        // i1 is conditioning set, i2 is response        
        //arma::span i1 = span(0,bsize-2);
        arma::span i2 = span(bsize-1,bsize-1);
        
        // get last row of cholmat
        arma::vec onevec = zeros(bsize);
        onevec(bsize-1) = 1.0;
        arma::vec choli2;
        if(grad_info){
            choli2 = solve( trimatu(cholmat.t()), onevec );
        }
        
        bool cond = bsize > 1;

        //double fac = 1.0;
        
        // do solves with X and y
        arma::mat LiX0;
        if(profbeta){
            LiX0 = solve( trimatl(cholmat), X0 );
        }
        arma::vec Liy0 = solve( trimatl(cholmat), ysub );
        
        // loglik objects
        l_logdet += 2.0*std::log( as_scalar(cholmat(i2,i2)) ); 
        l_ySy +=    pow( as_scalar(Liy0(i2)), 2 );
        if(profbeta){
            l_XSX +=   LiX0.rows(i2).t() * LiX0.rows(i2);
            l_ySX += ( Liy0(i2) * LiX0.rows(i2) ).t();
        }
        
        if( grad_info ){
        // gradient objects
        // LidSLi3 is last column of Li * (dS_j) * Lit for 1 parameter i
        // LidSLi2 stores these columns in a matrix for all parameters
        arma::mat LidSLi2(bsize,nparms);
        
        if(cond){ // if we condition on anything
            
            for(uword j=0; j<nparms; j++){
                // compute last column of Li * (dS_j) * Lit
                arma::vec LidSLi3 = solve( trimatl(cholmat), dcovmat.slice(j) * choli2 );
                // store LiX0.t() * LidSLi3 and Liy0.t() * LidSLi3
                arma::vec v1 = LiX0.t() * LidSLi3;
                double s1 = as_scalar( Liy0.t() * LidSLi3 ); 
                // update all quantities
                // bottom-right corner gets double counted, so need to subtract it off
                l_dXSX.slice(j) += v1 * LiX0.rows(i2) + ( v1 * LiX0.rows(i2) ).t() - 
                    as_scalar(LidSLi3(i2)) * ( LiX0.rows(i2).t() * LiX0.rows(i2) );
                l_dySy(j) += as_scalar( 2.0 * s1 * Liy0(i2)  - 
                    LidSLi3(i2) * Liy0(i2) * Liy0(i2) );
                l_dySX.col(j) += (  s1 * LiX0.rows(i2) + ( v1 * Liy0(i2) ).t() -  
                    as_scalar( LidSLi3(i2) ) * LiX0.rows(i2) * as_scalar( Liy0(i2))).t();
                l_dlogdet(j) += as_scalar( LidSLi3(i2) );
                // store last column of Li * (dS_j) * Lit
                LidSLi2.col(j) = LidSLi3;
            }

            // fisher information object
            // bottom right corner gets double counted, so subtract it off
            for(uword i=0; i<nparms; i++){ 
	      for(uword j=0; j<i+1; j++){
                 l_ainfo(i,j) += 1.0*accu( LidSLi2.col(i) % LidSLi2.col(j) ) - 0.5*accu( LidSLi2.rows(i2).col(j) % LidSLi2.rows(i2).col(i) );
            }}
            
        } else { // similar calculations, but for when there is no conditioning set
            for(uword j=0; j<nparms; j++){
                arma::mat LidSLi = solve( trimatl(cholmat), dcovmat.slice(j) );
                LidSLi = solve( trimatl(cholmat), LidSLi.t() );
                l_dXSX.slice(j) += LiX0.t() *  LidSLi * LiX0; 
                l_dySy(j) += as_scalar( Liy0.t() * LidSLi * Liy0 );
                l_dySX.col(j) += ( ( Liy0.t() * LidSLi ) * LiX0 ).t();
                l_dlogdet(j) += trace( LidSLi );
                LidSLi2.col(j) = LidSLi;
            }
            
            // fisher information object
            for(uword i=0; i<nparms; i++){ 
	      for(uword j=0; j<i+1; j++){
                l_ainfo(i,j) += 0.5*accu( LidSLi2.col(i) % LidSLi2.col(j) ); 
            }}

        }
        
        }

    }
#pragma omp critical
{  
    *XSX += l_XSX;
    *ySX += l_ySX;
    *ySy += l_ySy;
    *logdet += l_logdet;
    *dXSX += l_dXSX;
    *dySX += l_dySX;
    *dySy += l_dySy;
    *dlogdet += l_dlogdet;
    *ainfo += l_ainfo;
}
}

}
    
    
    
    
void compute_pieces_grouped(
    arma::vec covparms, 
    StringVector covfun_name,
    arma::mat locs, 
    List NNlist,
    arma::vec y, 
    arma::mat X,
    mat* XSX,
    vec* ySX,
    double* ySy,
    double* logdet,
    arma::cube* dXSX,
    arma::mat* dySX,
    arma::vec* dySy,
    arma::vec* dlogdet,
    arma::mat* ainfo,
    bool profbeta,
    bool grad_info
){

    // data dimensions
    uword p = X.n_cols;
    uword nparms = covparms.n_elem;
    uword dim = locs.n_cols;
    

    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];
    
    // assign covariance fun and derivative based on covfun_name_string
    covfun_t covstruct = get_covfun(covfun_name_string);
    mat (*p_covfun)(arma::vec, arma::mat) = covstruct.p_covfun;
    cube (*p_d_covfun)(arma::vec, arma::mat) = covstruct.p_d_covfun;

    // vector of all indices
    arma::uvec all_inds = NNlist["all_inds"];
    // vector of local response indices
    arma::uvec local_resp_inds = as<arma::uvec>(NNlist["local_resp_inds"]);
    // vector of global response indices
    arma::uvec global_resp_inds = as<arma::uvec>(NNlist["global_resp_inds"]);
    // last index of each block in all_inds
    arma::uvec last_ind_of_block = as<arma::uvec>(NNlist["last_ind_of_block"]);
    // last response index of each block in local_resp_inds and global_resp_inds
    arma::uvec last_resp_of_block = as<arma::uvec>(NNlist["last_resp_of_block"]);

    uword nb = last_ind_of_block.n_elem;  // number of blocks

    arma::mat l_XSX = arma::mat(p, p, fill::zeros);
    arma::vec l_ySX = arma::vec(p, fill::zeros);
    double l_ySy = 0.0;
    double l_logdet = 0.0;
    arma::cube l_dXSX = arma::cube(p,p,nparms,fill::zeros);
    arma::mat l_dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec l_dySy = arma::vec(nparms, fill::zeros);
    arma::vec l_dlogdet = arma::vec(nparms, fill::zeros);
    arma::mat l_ainfo = arma::mat(nparms, nparms, fill::zeros);
    

    // loop over every block
    for(uword i=0; i<nb; i++){

        // first ind and last ind are the positions in all_inds
        // of the observations for block i.
        // these come in 1-indexing and are converted to 0-indexing here
        uword first_ind;
        uword last_ind;
        if(i==0){ first_ind = 0; } else {first_ind = last_ind_of_block[i-1]+1-1; }
        last_ind = last_ind_of_block[i]-1;
        uword bsize = last_ind - first_ind + 1;
        
        uword first_resp;
        uword last_resp;
        if(i==0){ first_resp = 0; } else {first_resp = last_resp_of_block[i-1]+1-1; }
        last_resp = last_resp_of_block[i] - 1;
        uword rsize = last_resp - first_resp + 1;
        
        arma::uvec whichresp(rsize);
        for(uword j=0; j<rsize; j++){
            whichresp(j) = local_resp_inds(first_resp+j) - 1;
        }

        // fill in ysub, locsub, and X0 in forward order
        arma::mat locsub(bsize, dim);
        arma::vec ysub(bsize);
        arma::mat X0(bsize,p);
        for(uword j=0; j<bsize; j++){
            uword jglobal = all_inds[first_ind + j] - 1;
            ysub(j) = y( jglobal );
            for(uword k=0;k<dim;k++){ locsub(j,k) = locs( jglobal, k ); }
            if(profbeta){
                for(uword k=0;k<p;k++){ X0(j,k) = X( jglobal, k ); }
            }
        }

        // compute covariance matrix and derivatives and take cholesky
        arma::mat covmat(bsize,bsize);
        covmat = (*p_covfun)( covparms, locsub );
        arma::cube dcovmat(bsize,bsize,nparms);
        if(grad_info){
            dcovmat = (*p_d_covfun)( covparms, locsub );
        }

        arma::mat cholmat = eye( size(covmat) );
        chol( cholmat, covmat, "lower" );

        // get response rows of inverse cholmat, put in column vectors
        arma::mat onemat = zeros(bsize,rsize);
        for(uword j=0; j<rsize; j++){ 
            onemat(whichresp(j),j) = 1.0;
        }

        arma::mat choli2(bsize,bsize);
        if(grad_info){
            choli2 = solve( trimatu(cholmat.t()), onemat );
        }

        bool cond = bsize > rsize;
        //double fac = 1.0;
        
        // do solves with X and y
        arma::mat LiX0( bsize, p );
        if(profbeta){
            LiX0 = solve( trimatl(cholmat), X0 );
        }

        arma::vec Liy0 = solve( trimatl(cholmat), ysub );
        // loglik objects
        for(uword j=0; j<rsize; j++){
            uword ii = whichresp(j);
            l_logdet += 2.0*std::log( as_scalar(cholmat(ii,ii)) ); 
            l_ySy +=    pow( as_scalar(Liy0(ii)), 2 );
            if(profbeta){
                l_XSX +=    LiX0.row(ii).t() * LiX0.row(ii);
                l_ySX += ( Liy0(ii) * LiX0.row(ii) ).t();
            }
        }    
        
        if(grad_info){
        // gradient objects
        // LidSLi3 is last column of Li * (dS_j) * Lit for 1 parameter i
        // LidSLi2 stores these columns in a matrix for all parameters
        if(cond){ // if we condition on anything
            
            arma::cube LidSLi2 = arma::cube(bsize,rsize,nparms,fill::zeros);
            for(uword j=0; j<nparms; j++){
                // compute last column of Li * (dS_j) * Lit
                arma::mat LidSLi4 = 
                    solve( trimatl(cholmat), 
                           dcovmat.slice(j) * choli2 );
                
                for(uword k=0; k<rsize; k++){
                    uword i2 = whichresp(k);
                    arma::span i1 = span(0,i2);
                    arma::vec LidSLi3 = LidSLi4(i1,k); 
                    // store LiX0.t() * LidSLi3 and Liy0.t() * LidSLi3
                    arma::vec v1 = LiX0.rows(i1).t() * LidSLi3;
                    double s1 = dot( Liy0(i1), LidSLi3 ); 
                    // update all quantities
                    // bottom-right corner gets double counted, so need to subtract it off
                    (l_dXSX).slice(j) += v1 * LiX0.row(i2) + ( v1 * LiX0.row(i2) ).t() - 
                        LidSLi3(i2) * LiX0.row(i2).t() * LiX0.row(i2);
                    (l_dySy)(j) += 2.0 * s1 * Liy0(i2)  - 
                        LidSLi3(i2) * Liy0(i2) * Liy0(i2);
                    (l_dySX).col(j) += (  s1 * LiX0.row(i2) + ( v1 * Liy0(i2) ).t() -  
                        LidSLi3(i2) * LiX0.row(i2) * Liy0(i2) ).t();
                    (l_dlogdet)(j) += LidSLi3(i2);
                    // store last column of Li * (dS_j) * Lit
                    LidSLi2.subcube(i1, span(k,k), span(j,j)) = LidSLi3;
                }
            }

            // fisher information object
            // bottom right corner gets double counted, so subtract it off
            for(uword i=0; i<nparms; i++){ for(uword j=0; j<i+1; j++){
                (l_ainfo)(i,j) += accu( LidSLi2.slice(i) % LidSLi2.slice(j) );
                for(uword k=0; k<rsize; k++){
                    uword i2 = whichresp(k);
                    (l_ainfo)(i,j) -= 0.5*LidSLi2(i2,k,j) * LidSLi2(i2,k,i);
                }
            }}
        } else { // similar calculations, but for when there is no conditioning set
            arma::cube LidSLi2(bsize,bsize,nparms);
            for(uword j=0; j<nparms; j++){
                arma::mat LidSLi = solve( trimatl(cholmat), dcovmat.slice(j) );
                LidSLi = solve( trimatl(cholmat), LidSLi.t() );
                (l_dXSX).slice(j) += LiX0.t() *  LidSLi * LiX0; 
                (l_dySy)(j) += as_scalar( Liy0.t() * LidSLi * Liy0 );
                (l_dySX).col(j) += ( ( Liy0.t() * LidSLi ) * LiX0 ).t();
                (l_dlogdet)(j) += trace( LidSLi );
                LidSLi2.slice(j) = LidSLi;
            }
            
            // fisher information object
            for(uword i=0; i<nparms; i++){ for(uword j=0; j<i+1; j++){
                (l_ainfo)(i,j) += 0.5*accu( LidSLi2.slice(i) % LidSLi2.slice(j) ); 
            }}
        }
        }
    }

    *XSX += l_XSX;
    *ySX += l_ySX;
    *ySy += l_ySy;
    *logdet += l_logdet;
    *dXSX += l_dXSX;
    *dySX += l_dySX;
    *dySy += l_dySy;
    *dlogdet += l_dlogdet;
    *ainfo += l_ainfo;
}    

    
    

void synthesize(
    NumericVector covparms, 
    StringVector covfun_name,
    NumericMatrix locs, 
    NumericMatrix NNarray,
    NumericVector y, 
    NumericMatrix X,
    NumericVector* ll, 
    NumericVector* betahat,
    NumericVector* grad,
    NumericMatrix* info,
    NumericMatrix* betainfo,
    bool profbeta,
    bool grad_info ){

    // data dimensions
    uword n = static_cast<uword>(y.length());
    uword p = static_cast<uword>(X.ncol());
    uword dim = static_cast<uword>(locs.ncol());
    uword m = static_cast<uword>(NNarray.ncol());
    uword nparms = static_cast<uword>(covparms.length());
    
    
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
    // first convert Numeric- to arma
    arma::vec covparms_c = arma::vec(covparms.begin(), nparms);
    arma::mat locs_c = arma::mat(locs.begin(), n, dim);
    arma::Mat<uword> NNarray_c = conv_to<Mat<uword>>::from(arma::mat(NNarray.begin(), n, m));
    arma::vec y_c = arma::vec(y.begin(), n);
    arma::mat X_c = arma::mat(X.begin(), n, p);
    
    compute_pieces(
        covparms_c, covfun_name, locs_c, NNarray_c, y_c, X_c,
        &XSX, &ySX, &ySy, &logdet, &dXSX, &dySX, &dySy, &dlogdet, &ainfo,
        profbeta, grad_info
    );
        
    // synthesize everything and update loglik, grad, beta, betainfo, info
    
    // betahat and dbeta
    arma::vec abeta = arma::vec( p, fill::zeros );
    if(profbeta){ abeta = solve( XSX, ySX ); }
    for(uword j=0; j<p; j++){ (*betahat)(j) = abeta(j); };

    arma::mat dbeta = arma::mat(p,nparms, fill::zeros);
    if( profbeta && grad_info){
    for(uword j=0; j<nparms; j++){
        dbeta.col(j) = solve( XSX, dySX.col(j) - dXSX.slice(j) * abeta );
    }
    }
    // get sigmahatsq
    double sig2 = ( ySy - 2.0*as_scalar( ySX.t() * abeta ) + 
        as_scalar( abeta.t() * XSX * abeta ) )/n;
    // loglikelihood
    (*ll)(0) = -0.5*( n*std::log(2.0*M_PI) + logdet + n*sig2 ); 
    
    if(profbeta){
    // betainfo
    for(uword i=0; i<p; i++){ for(uword j=0; j<i+1; j++){
        (*betainfo)(i,j) = XSX(i,j);
        (*betainfo)(j,i) = XSX(j,i);
    }}
    }

    if(grad_info){
    // gradient
    for(uword j=0; j<nparms; j++){
        (*grad)(j) = 0.0;
        (*grad)(j) -= 0.5*dlogdet(j);
        (*grad)(j) += 0.5*dySy(j);
        (*grad)(j) -= 1.0*as_scalar( abeta.t() * dySX.col(j) );
        (*grad)(j) += 1.0*as_scalar( ySX.t() * dbeta.col(j) );
        (*grad)(j) += 0.5*as_scalar( abeta.t() * dXSX.slice(j) * abeta );
        (*grad)(j) -= 1.0*as_scalar( abeta.t() * XSX * dbeta.col(j) );
    }
    // fisher information
    for(uword i=0; i<nparms; i++){ for(uword j=0; j<i+1; j++){
        (*info)(i,j) = ainfo(i,j);
        (*info)(j,i) = (*info)(i,j);
    }}
    }
}


void synthesize_grouped(
    NumericVector covparms, 
    StringVector covfun_name,
    NumericMatrix locs, 
    List NNlist,
    NumericVector y, 
    NumericMatrix X,
    NumericVector* ll, 
    NumericVector* betahat,
    NumericVector* grad,
    NumericMatrix* info,
    NumericMatrix* betainfo,
    bool profbeta,
    bool grad_info){

    // data dimensions
    uword n = static_cast<uword>(y.length());
    uword p = static_cast<uword>(X.ncol());
    uword nparms = static_cast<uword>(covparms.length());
    uword dim = static_cast<uword>(locs.ncol());
    
    // likelihood objects
    arma::mat XSX = arma::mat(p, p, fill::zeros);
    arma::vec ySX = arma::vec(p, fill::zeros);
    double ySy = 0.0;
    double logdet = 0.0;
    
    // gradient objects    
    arma::cube dXSX = arma::cube(p, p,nparms,fill::zeros);
    arma::mat dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec dySy = arma::vec(nparms, fill::zeros);
    arma::vec dlogdet = arma::vec(nparms, fill::zeros);
    
    // fisher information
    arma::mat ainfo = arma::mat(nparms, nparms, fill::zeros);
    
    // this is where the big computation happens
    
    // convert Numeric- to arma
    
    arma::vec covparms_c = arma::vec(covparms.begin(), nparms);
    arma::mat locs_c = arma::mat(locs.begin(), n, dim);
    arma::vec y_c = arma::vec(y.begin(), n);
    arma::mat X_c = arma::mat(X.begin(), n, p);
    
    
    compute_pieces_grouped(
        covparms_c, covfun_name, locs_c, NNlist, y_c, X_c,
        &XSX, &ySX, &ySy, &logdet, &dXSX, &dySX, &dySy, &dlogdet, &ainfo,
        profbeta, grad_info
    );
        
    // synthesize everything and update loglik, grad, beta, betainfo, info
    
    // betahat and dbeta
    arma::vec abeta = arma::vec(p,fill::zeros);
    if(profbeta){ abeta = solve( XSX, ySX ); }
    for(uword j=0; j<p; j++){ (*betahat)(j) = abeta(j); };
    
    arma::mat dbeta(p,nparms);
    if(profbeta && grad_info){
    for(uword j=0; j<nparms; j++){
        dbeta.col(j) = solve( XSX, dySX.col(j) - dXSX.slice(j) * abeta );
    }
    }

    // get sigmahatsq
    double sig2 = ( ySy - 2.0*as_scalar( ySX.t() * abeta ) + 
        as_scalar( abeta.t() * XSX * abeta ) )/n;
    // loglikelihood
    (*ll)(0) = -0.5*( n*std::log(2.0*M_PI) + logdet + n*sig2 ); 
    
    if(profbeta){    
    // betainfo
    for(uword i=0; i<p; i++){ for(uword j=0; j<i+1; j++){
        (*betainfo)(i,j) = XSX(i,j);
        (*betainfo)(j,i) = XSX(j,i);
    }}
    }
    
    if(grad_info){
    // gradient
    for(uword j=0; j<nparms; j++){
        (*grad)(j) = 0.0;
        (*grad)(j) -= 0.5*dlogdet(j);
        (*grad)(j) += 0.5*dySy(j);
        (*grad)(j) -= 1.0*as_scalar( abeta.t() * dySX.col(j) );
        (*grad)(j) += 1.0*as_scalar( ySX.t() * dbeta.col(j) );
        (*grad)(j) += 0.5*as_scalar( abeta.t() * dXSX.slice(j) * abeta );
        (*grad)(j) -= 1.0*as_scalar( abeta.t() * XSX * dbeta.col(j) );
    }
    // fisher information
    for(uword i=0; i<nparms; i++){ for(uword j=0; j<i+1; j++){
        (*info)(i,j) = ainfo(i,j);
        (*info)(j,i) = (*info)(i,j);
    }}
    }
}


//' Entries of inverse Cholesky approximation
//' 
//' This function returns the entries of the inverse Cholesky
//' factor of the covariance matrix implied by Vecchia's approximation.
//' For return matrix \code{Linv}, \code{Linv[i,]} contains 
//' the non-zero entries of row \code{i} of
//' the inverse Cholesky matrix. The columns of the non-zero entries
//' are specified in \code{NNarray[i,]}.
//' @inheritParams vecchia_meanzero_loglik
//' @param start_ind Compute entries of Linv only for rows \code{start_ind}
//' until the last row. 
//' @return matrix containing entries of inverse Cholesky
//' @examples
//' n1 <- 40
//' n2 <- 40
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' covparms <- c(2, 0.2, 0.75, 0)
//' NNarray <- find_ordered_nn(locs,20)
//' Linv <- vecchia_Linv(covparms, "matern_isotropic", locs, NNarray)
//' @export
// [[Rcpp::export]]
NumericMatrix vecchia_Linv(
    NumericVector covparms,
    StringVector covfun_name,
    NumericMatrix locs,
    IntegerMatrix NNarray, 
    int start_ind = 1){
    
    int n = locs.nrow();
    int m = NNarray.ncol();
    int dim = locs.ncol();;
    NumericMatrix Linv(n,m);
    
    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];
    
    // assign covariance fun and derivative based on covfun_name_string
    covfun_t covstruct = get_covfun(covfun_name_string);
    mat (*p_covfun)(arma::vec, arma::mat) = covstruct.p_covfun;
    //cube (*p_d_covfun)(NumericVector, NumericMatrix) = covstruct.p_d_covfun;


    // loop over every observation    
    for(int i=start_ind-1; i<n; i++){
    
        Rcpp::checkUserInterrupt();
        int bsize = std::min(i+1,m);

        // first, fill in ysub, locsub, and X0 in reverse order
        arma::mat locsub(bsize, dim);
        for(int j=0; j<bsize; j++){
            for(int k=0;k<dim;k++){ locsub(j,k) = locs( NNarray(i,bsize-1-j) - 1, k ); }
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
        for(int j=0; j<bsize; j++){
            Linv(i,j) = choli2(bsize-1-j);
        }
    }    
    return Linv;    
}


#endif
    
