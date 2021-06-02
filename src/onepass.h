#ifndef ONEPASS_H
#define ONEPASS_H

#include <RcppArmadillo.h>
#include <string>
//[[Rcpp::depends(RcppArmadillo)]]

#include "covmatrix_funs.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

arma::vec forward_solve( arma::mat cholmat, arma::vec b ){

    int n = cholmat.n_rows;
    arma::vec x(n);
    x(0) = b(0)/cholmat(0,0);

    for(int i=1; i<n; i++){
        double dd = 0.0;
        for(int j=0; j<i; j++){
            dd += cholmat(i,j)*x(j);
        }
        x(i) = (b(i)-dd)/cholmat(i,i);
    }    
    return x;

} 

arma::mat forward_solve_mat( arma::mat cholmat, arma::mat b ){

    int n = cholmat.n_rows;
    int p = b.n_cols;
    arma::mat x(n,p);
    for(int k=0; k<p; k++){ x(0,k) = b(0,k)/cholmat(0,0); }

    for(int i=1; i<n; i++){
	for(int k=0; k<p; k++){
            double dd = 0.0;
            for(int j=0; j<i; j++){
                dd += cholmat(i,j)*x(j,k);
            }
            x(i,k) = (b(i,k)-dd)/cholmat(i,i);
       	}
    }    
    return x;
} 

arma::vec backward_solve( arma::mat lower, arma::vec b ){

    int n = lower.n_rows;
    arma::vec x(n);
    x(n-1) = b(n-1)/lower(n-1,n-1);

    for(int i=n-2; i>=0; i--){
        double dd = 0.0;
        for(int j=n-1; j>i; j--){
            dd += lower(j,i)*x(j);
        }
        x(i) = (b(i)-dd)/lower(i,i);
    }    
    return x;
} 

arma::mat backward_solve_mat( arma::mat cholmat, arma::mat b ){

    int n = cholmat.n_rows;
    int p = b.n_cols;
    arma::mat x(n,p);
    for(int k=0; k<p; k++){ x(n-1,k) = b(n-1,k)/cholmat(n-1,n-1); }

    for(int i=n-2; i>=0; i--){
	for(int k=0; k<p; k++){
            double dd = 0.0;
            for(int j=n-1; j>i; j--){
                dd += cholmat(j,i)*x(j,k);
            }
            x(i,k) = (b(i,k)-dd)/cholmat(i,i);
       	}
    }    
    return x;
} 


arma::mat mychol( arma::mat A ){

    arma::uword n = A.n_rows;
    arma::mat L(n,n);
    bool pd = true;
    
    // upper-left entry
    if( A(0,0) < 0 ){
	pd = false;
	L(0,0) = 1.0;
    } else {
        L(0,0) = std::sqrt(A(0,0));
    }
    if( n > 1 ){
	// second row
	L(1,0) = A(1,0)/L(0,0);
	double f = A(1,1) - L(1,0)*L(1,0);
	if( f < 0 ){
	    pd = false;
	    L(1,1) = 1.0;
	} else {
	    L(1,1) = std::sqrt( f );
	}
	// rest of the rows
	if( n > 2 ){
            for(uword i=2; i<n; i++){
    	        // leftmost entry in row i
    	        L(i,0) = A(i,0)/L(0,0);
    	        // middle entries in row i 
    	        for(uword j=1; j<i; j++){
    	            double d = A(i,j);
    	            for(uword k=0; k<j; k++){
    	        	d -= L(i,k)*L(j,k);
    	            }
    	            L(i,j) = d/L(j,j);
    	        }
		// diagonal entry in row i
    	        double e = A(i,i);
    	        for(uword k=0; k<i; k++){
    	            e -= L(i,k)*L(i,k);
    	        }
		if( e < 0 ){
		    pd = false;
		    L(i,i) = 1.0;
		} else {
    	            L(i,i) = std::sqrt(e);
		}
	    }
	}
    }
    return L;	
}



void compute_pieces(
    arma::vec covparms, 
    StringVector covfun_name,
    arma::mat locs, 
    arma::mat NNarray,
    arma::vec y, 
    arma::mat X,
    mat* XSX,
    vec* ySX,
    double* ySy,
    double* logdet,
    cube* dXSX,
    mat* dySX,
    vec* dySy,
    vec* dlogdet,
    mat* ainfo,
    int profbeta,
    int grad_info
){

    // data dimensions
    int n = y.n_elem;
    int m = NNarray.n_cols;
    int p = X.n_cols;
    int nparms = covparms.n_elem;
    int dim = locs.n_cols;
    
    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];
    
    // assign covariance fun and derivative based on covfun_name_string

    /* p_covfun is an array of length 1. Its entry is a pointer to a function which takes
     in arma::vec and arma::mat and returns mat. p_d_covfun is analogous. This was a workaround for the solaris bug*/

    mat (*p_covfun[1])(arma::vec, arma::mat);
    cube (*p_d_covfun[1])(arma::vec, arma::mat);
    get_covfun(covfun_name_string, p_covfun, p_d_covfun);
    

#pragma omp parallel 
{   
    arma::mat l_XSX = arma::mat(p, p, fill::zeros);
    arma::vec l_ySX = arma::vec(p, fill::zeros);
    double l_ySy = 0.0;
    double l_logdet = 0.0;
    arma::cube l_dXSX = arma::cube(p,p, nparms, fill::zeros);
    arma::mat l_dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec l_dySy = arma::vec(nparms, fill::zeros);
    arma::vec l_dlogdet = arma::vec(nparms, fill::zeros);
    arma::mat l_ainfo = arma::mat(nparms, nparms, fill::zeros);

    #pragma omp for	    
    for(int i=0; i<n; i++){
    
        int bsize = std::min(i+1,m);

	//std::vector<std::chrono::steady_clock::time_point> tt;

	//tt.push_back( std::chrono::steady_clock::now() );

        // first, fill in ysub, locsub, and X0 in reverse order
        arma::mat locsub(bsize, dim);
        arma::vec ysub(bsize);
        arma::mat X0( bsize, p );
        for(int j=bsize-1; j>=0; j--){
            ysub(bsize-1-j) = y( NNarray(i,j)-1 );
            for(int k=0;k<dim;k++){ locsub(bsize-1-j,k) = locs( NNarray(i,j)-1, k ); }
            if(profbeta){ 
                for(int k=0;k<p;k++){ X0(bsize-1-j,k) = X( NNarray(i,j)-1, k ); } 
            }
        }
        
        // compute covariance matrix and derivatives and take cholesky
        arma::mat covmat = p_covfun[0]( covparms, locsub );	
	
        arma::cube dcovmat;
        if(grad_info){ 
            dcovmat = p_d_covfun[0]( covparms, locsub ); 
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
            //choli2 = solve( trimatu(cholmat.t()), onevec );
            choli2 = backward_solve( cholmat, onevec );
        }
        
        bool cond = bsize > 1;
        //double fac = 1.0;
        
        // do solves with X and y
        arma::mat LiX0;
        if(profbeta){
            //LiX0 = solve( trimatl(cholmat), X0 );
            LiX0 = forward_solve_mat( cholmat, X0 );
        }

        //arma::vec Liy0 = solve( trimatl(cholmat), ysub );
        arma::vec Liy0 = forward_solve( cholmat, ysub );
        
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
            
            for(int j=0; j<nparms; j++){
                // compute last column of Li * (dS_j) * Lit
                //arma::vec LidSLi3 = solve( trimatl(cholmat), dcovmat.slice(j) * choli2 );
                arma::vec LidSLi3 = forward_solve( cholmat, dcovmat.slice(j) * choli2 );
                // store LiX0.t() * LidSLi3 and Liy0.t() * LidSLi3
                arma::vec v1 = LiX0.t() * LidSLi3;
                double s1 = as_scalar( Liy0.t() * LidSLi3 ); 
                // update all quantities
                // bottom-right corner gets double counted, so need to subtract it off
                (l_dXSX).slice(j) += v1 * LiX0.rows(i2) + ( v1 * LiX0.rows(i2) ).t() - 
                    as_scalar(LidSLi3(i2)) * ( LiX0.rows(i2).t() * LiX0.rows(i2) );
                (l_dySy)(j) += as_scalar( 2.0 * s1 * Liy0(i2)  - 
                    LidSLi3(i2) * Liy0(i2) * Liy0(i2) );
                (l_dySX).col(j) += (  s1 * LiX0.rows(i2) + ( v1 * Liy0(i2) ).t() -  
                    as_scalar( LidSLi3(i2) ) * LiX0.rows(i2) * as_scalar( Liy0(i2))).t();
                (l_dlogdet)(j) += as_scalar( LidSLi3(i2) );
                // store last column of Li * (dS_j) * Lit
                LidSLi2.col(j) = LidSLi3;
            }

            // fisher information object
            // bottom right corner gets double counted, so subtract it off
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (l_ainfo)(i,j) += 
                    1.0*accu( LidSLi2.col(i) % LidSLi2.col(j) ) - 
                    0.5*accu( LidSLi2.rows(i2).col(j) %
                              LidSLi2.rows(i2).col(i) );
            }}
            
        } else { // similar calculations, but for when there is no conditioning set
            for(int j=0; j<nparms; j++){
                //arma::mat LidSLi = solve( trimatl(cholmat), dcovmat.slice(j) );
                arma::mat LidSLi = forward_solve_mat( cholmat, dcovmat.slice(j) );
                //LidSLi = solve( trimatl(cholmat), LidSLi.t() );
                LidSLi = forward_solve_mat( cholmat, LidSLi.t() );
                (l_dXSX).slice(j) += LiX0.t() *  LidSLi * LiX0; 
                (l_dySy)(j) += as_scalar( Liy0.t() * LidSLi * Liy0 );
                (l_dySX).col(j) += ( ( Liy0.t() * LidSLi ) * LiX0 ).t();
                (l_dlogdet)(j) += trace( LidSLi );
                LidSLi2.col(j) = LidSLi;
            }
            
            // fisher information object
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (l_ainfo)(i,j) += 0.5*accu( LidSLi2.col(i) % LidSLi2.col(j) ); 
            }}

        }
        
        }

	// if(i % 100 == 0 ){
	//     for(int k=1; k<tt.size(); k++ ){
	//         cout << std::chrono::duration_cast<std::chrono::microseconds>(tt[k]-tt[k-1]).count();
	//         cout << "  ";
	//     }
	//     cout << endl;
	// }


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
    cube* dXSX,
    mat* dySX,
    vec* dySy,
    vec* dlogdet,
    mat* ainfo,
    bool profbeta,
    bool grad_info
){

    // data dimensions
    //int n = y.length();
    //int m = NNarray.ncol();
    int p = X.n_cols;
    int nparms = covparms.n_elem;
    int dim = locs.n_cols;
    

    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];
    
    // assign covariance fun and derivative based on covfun_name_string
    mat (*p_covfun[1])(arma::vec, arma::mat);
    cube (*p_d_covfun[1])(arma::vec, arma::mat);
    get_covfun(covfun_name_string, p_covfun, p_d_covfun);

    // vector of all indices
    arma::vec all_inds = NNlist["all_inds"];
    // vector of local response indices
    arma::vec local_resp_inds = as<arma::vec>(NNlist["local_resp_inds"]);
    // vector of global response indices
    arma::vec global_resp_inds = as<arma::vec>(NNlist["global_resp_inds"]);
    // last index of each block in all_inds
    arma::vec last_ind_of_block = as<arma::vec>(NNlist["last_ind_of_block"]);
    // last response index of each block in local_resp_inds and global_resp_inds
    arma::vec last_resp_of_block = as<arma::vec>(NNlist["last_resp_of_block"]);

    int nb = last_ind_of_block.n_elem;  // number of blocks
    
#pragma omp parallel 
{   
    arma::mat l_XSX = arma::mat(p, p, fill::zeros);
    arma::vec l_ySX = arma::vec(p, fill::zeros);
    double l_ySy = 0.0;
    double l_logdet = 0.0;
    arma::cube l_dXSX = arma::cube(p,p, nparms, fill::zeros);
    arma::mat l_dySX = arma::mat(p, nparms, fill::zeros);
    arma::vec l_dySy = arma::vec(nparms, fill::zeros);
    arma::vec l_dlogdet = arma::vec(nparms, fill::zeros);
    arma::mat l_ainfo = arma::mat(nparms, nparms, fill::zeros);

#pragma omp for	    
    for(int i=0; i<nb; i++){

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
        arma::mat locsub(bsize, dim);
        arma::vec ysub(bsize);
        arma::mat X0(bsize,p);
        for(int j=0; j<bsize; j++){
            int jglobal = all_inds[first_ind + j] - 1;
            ysub(j) = y( jglobal );
            for(int k=0;k<dim;k++){ locsub(j,k) = locs( jglobal, k ); }
            if(profbeta){
                for(int k=0;k<p;k++){ X0(j,k) = X( jglobal, k ); }
            }
        }

        // compute covariance matrix and derivatives and take cholesky
	arma::mat covmat = p_covfun[0]( covparms, locsub );

        arma::cube dcovmat(bsize,bsize,nparms);
        if(grad_info){
            dcovmat = p_d_covfun[0]( covparms, locsub );
        }

        //arma::mat cholmat = eye( size(covmat) );
        //chol( cholmat, covmat, "lower" );
	arma::mat cholmat = mychol(covmat);
	
        // get response rows of inverse cholmat, put in column vectors
        arma::mat onemat = zeros(bsize,rsize);
        for(int j=0; j<rsize; j++){ 
            onemat(whichresp(j),j) = 1.0;
        }

        arma::mat choli2(bsize,bsize);
        if(grad_info){
            //choli2 = solve( trimatu(cholmat.t()), onemat );
            choli2 = backward_solve_mat( cholmat, onemat );
        }

        bool cond = bsize > rsize;
        //double fac = 1.0;
        
        // do solves with X and y
        arma::mat LiX0( bsize, p );
        if(profbeta){
            //LiX0 = solve( trimatl(cholmat), X0 );
            LiX0 = forward_solve_mat( cholmat, X0 );
        }

        //arma::vec Liy0 = solve( trimatl(cholmat), ysub );
        arma::vec Liy0 = forward_solve( cholmat, ysub );

        // loglik objects
        for(int j=0; j<rsize; j++){
            int ii = whichresp(j);
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
            for(int j=0; j<nparms; j++){
                // compute last column of Li * (dS_j) * Lit
                //arma::mat LidSLi4 = solve( trimatl(cholmat), dcovmat.slice(j) * choli2 );
                arma::mat LidSLi4 = forward_solve_mat( cholmat, dcovmat.slice(j) * choli2 );
                
                for(int k=0; k<rsize; k++){
                    int i2 = whichresp(k);
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
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (l_ainfo)(i,j) += accu( LidSLi2.slice(i) % LidSLi2.slice(j) );
                for(int k=0; k<rsize; k++){
                    int i2 = whichresp(k);
                    (l_ainfo)(i,j) -= 0.5*LidSLi2(i2,k,j) * LidSLi2(i2,k,i);
                }
            }}
        } else { // similar calculations, but for when there is no conditioning set
            arma::cube LidSLi2(bsize,bsize,nparms);
            for(int j=0; j<nparms; j++){
                //arma::mat LidSLi = solve( trimatl(cholmat), dcovmat.slice(j) );
                arma::mat LidSLi = forward_solve_mat( cholmat, dcovmat.slice(j) );
                //LidSLi = solve( trimatl(cholmat), LidSLi.t() );
                LidSLi = forward_solve_mat( cholmat, LidSLi.t() );
                (l_dXSX).slice(j) += LiX0.t() *  LidSLi * LiX0; 
                (l_dySy)(j) += as_scalar( Liy0.t() * LidSLi * Liy0 );
                (l_dySX).col(j) += ( ( Liy0.t() * LidSLi ) * LiX0 ).t();
                (l_dlogdet)(j) += trace( LidSLi );
                LidSLi2.slice(j) = LidSLi;
            }
            
            // fisher information object
            for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
                (l_ainfo)(i,j) += 0.5*accu( LidSLi2.slice(i) % LidSLi2.slice(j) ); 
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

    
    

void synthesize(
    NumericVector covparms, 
    StringVector covfun_name,
    const NumericMatrix locs, 
    NumericMatrix NNarray,
    NumericVector& y, 
    NumericMatrix X,
    NumericVector* ll, 
    NumericVector* betahat,
    NumericVector* grad,
    NumericMatrix* info,
    NumericMatrix* betainfo,
    bool profbeta,
    bool grad_info ){

    // data dimensions
    int n = y.length();
    //int m = NNarray.ncol();
    int p = X.ncol();
    int nparms = covparms.length();
    //int dim = locs.ncol();
    
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
    arma::vec covparms_c = arma::vec(covparms.begin(),covparms.length());
    arma::mat locs_c = arma::mat(locs.begin(),locs.nrow(),locs.ncol());
    arma::mat NNarray_c = arma::mat(NNarray.begin(),NNarray.nrow(),NNarray.ncol());
    arma::vec y_c = arma::vec(y.begin(),y.length());
    arma::mat X_c = arma::mat(X.begin(),X.nrow(),X.ncol());
    
    compute_pieces(
        covparms_c, covfun_name, locs_c, NNarray_c, y_c, X_c,
        &XSX, &ySX, &ySy, &logdet, &dXSX, &dySX, &dySy, &dlogdet, &ainfo,
        profbeta, grad_info
    );
        
    // synthesize everything and update loglik, grad, beta, betainfo, info
    
    // betahat and dbeta
    arma::vec abeta = arma::vec( p, fill::zeros );
    if(profbeta){ abeta = solve( XSX, ySX ); }
    for(int j=0; j<p; j++){ (*betahat)(j) = abeta(j); };

    arma::mat dbeta = arma::mat(p,nparms, fill::zeros);
    if( profbeta && grad_info){
    for(int j=0; j<nparms; j++){
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
    for(int i=0; i<p; i++){ for(int j=0; j<i+1; j++){
        (*betainfo)(i,j) = XSX(i,j);
        (*betainfo)(j,i) = XSX(j,i);
    }}
    }

    if(grad_info){
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
    // fisher information
    for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
        (*info)(i,j) = ainfo(i,j);
        (*info)(j,i) = (*info)(i,j);
    }}
    }

}


void synthesize_grouped(
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
    NumericMatrix* betainfo,
    bool profbeta,
    bool grad_info){

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
    
    // convert Numeric- to arma
    arma::vec covparms_c = arma::vec(covparms.begin(),covparms.length());
    arma::mat locs_c = arma::mat(locs.begin(),locs.nrow(),locs.ncol());
    arma::vec y_c = arma::vec(y.begin(),y.length());
    arma::mat X_c = arma::mat(X.begin(),X.nrow(),X.ncol());
    
    
    compute_pieces_grouped(
        covparms_c, covfun_name, locs_c, NNlist, y_c, X_c,
        &XSX, &ySX, &ySy, &logdet, &dXSX, &dySX, &dySy, &dlogdet, &ainfo,
        profbeta, grad_info
    );
        
    // synthesize everything and update loglik, grad, beta, betainfo, info
    
    // betahat and dbeta
    arma::vec abeta = arma::vec(p,fill::zeros);
    if(profbeta){ abeta = solve( XSX, ySX ); }
    for(int j=0; j<p; j++){ (*betahat)(j) = abeta(j); };
    
    arma::mat dbeta(p,nparms);
    if(profbeta && grad_info){
    for(int j=0; j<nparms; j++){
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
    for(int i=0; i<p; i++){ for(int j=0; j<i+1; j++){
        (*betainfo)(i,j) = XSX(i,j);
        (*betainfo)(j,i) = XSX(j,i);
    }}
    }
    
    if(grad_info){
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
    // fisher information
    for(int i=0; i<nparms; i++){ for(int j=0; j<i+1; j++){
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
    arma::vec covparms,
    StringVector covfun_name,
    arma::mat locs,
    arma::mat NNarray, 
    int start_ind = 1){
    
    // data dimensions
    int n = locs.n_rows;
    int m = NNarray.n_cols;
    //int nparms = covparms.length();
    int dim = locs.n_cols;
    
    // convert StringVector to std::string to use .compare() below
    std::string covfun_name_string;
    covfun_name_string = covfun_name[0];
    
    // assign covariance fun and derivative based on covfun_name_string
    mat (*p_covfun[1])(arma::vec, arma::mat);
    cube (*p_d_covfun[1])(arma::vec, arma::mat);
    get_covfun(covfun_name_string, p_covfun, p_d_covfun);
    
    arma::mat Linv = arma::mat(n, m, fill::zeros);
    // loop over every observation    
#pragma omp parallel 
{
    arma::mat l_Linv = arma::mat(n, m, fill::zeros);
#pragma omp for
    for(int i=start_ind-1; i<n; i++){
    
        //Rcpp::checkUserInterrupt();
        int bsize = std::min(i+1,m);

        // first, fill in ysub, locsub, and X0 in reverse order
        arma::mat locsub(bsize, dim);
        for(int j=bsize-1; j>=0; j--){
            for(int k=0;k<dim;k++){ locsub(bsize-1-j,k) = locs( NNarray(i,j)-1, k ); }
        }
        
        // compute covariance matrix and derivatives and take cholesky
        arma::mat covmat = p_covfun[0]( covparms, locsub );
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
            l_Linv(i,bsize-1-j) = choli2(j);
        }
    }
#pragma omp critical
    Linv += l_Linv; 
}
    NumericMatrix Linv_r = wrap(Linv);
    return Linv_r;    
}


#endif
    
