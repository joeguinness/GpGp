
#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "arma_vecchia_fun.h"
#include "arma_onepass.h"
//[[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector arma_vecchia_loglik(
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    const NumericMatrix locs, 
    IntegerMatrix NNarray){

    NumericVector ll(1);        // loglikelihood to be returned
    NumericMatrix Linv(1,1);    // Linv not to be returned
    arma_vecchia(covparms, covfun_name, locs, NNarray, y, &Linv, &ll, 1);
    return ll;
}


// [[Rcpp::export]]
List arma_vecchia_grad_hess( 
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    IntegerMatrix NNarray ){
    
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    arma_onepass_synthesize(covparms, covfun_name, locs, NNarray, y, X,
        &ll, &betahat, &grad, &info, &betainfo 
    );
    
    List ret = List::create( Named("loglik") = ll, Named("betahat") = betahat,
        Named("grad") = grad, Named("info") = info, Named("betainfo") = betainfo );
    return ret;
        
}


// [[Rcpp::export]]
List arma_vecchia_grad_hess_grouped( 
    NumericVector covparms, 
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    List NNlist ){
    
    NumericVector ll(1);
    NumericVector grad( covparms.length() );
    NumericVector betahat( X.ncol() );
    NumericMatrix info( covparms.length(), covparms.length() );
    NumericMatrix betainfo( X.ncol(), X.ncol() );

    // this function calls arma_onepass_compute_pieces
    // then synthesizes the result into loglik, beta, grad, info, betainfo
    // maybe the synthesize functions should take in an argument that
    // says which compute_pieces function to use
    arma_onepass_synthesize_grouped(covparms, covfun_name, locs, NNlist, y, X,
        &ll, &betahat, &grad, &info, &betainfo 
    );
    
    
    List ret = List::create( Named("loglik") = ll, Named("betahat") = betahat,
        Named("grad") = grad, Named("info") = info, Named("betainfo") = betainfo );
    return ret;
        
}


// [[Rcpp::export]]
List arma_vecchia_profile_grad_hess( 
    NumericVector subparms, 
    StringVector covfun_name,
    NumericVector y,
    NumericMatrix X,
    const NumericMatrix locs,
    IntegerMatrix NNarray ){
    
    NumericVector ll(1);
    double sigmasq;
    NumericVector grad( subparms.length() );
    NumericVector betahat( X.ncol() );
    NumericMatrix info( subparms.length()+1, subparms.length()+1 );
    NumericMatrix betainfo( X.ncol(), X.ncol() );
    Rcout << "*" << endl;
    arma_onepass_profile_synthesize(subparms, &sigmasq, covfun_name, locs, NNarray, y, X,
        &ll, &betahat, &grad, &info, &betainfo 
    );
    
    int sublen = subparms.length();
    mat ainfo(sublen+1,sublen+1);
    for(int i=0; i<sublen+1; i++){ for(int j=0; j<sublen+1; j++){
        ainfo(i,j) = info(i,j);
    }}
    mat invinfo = inv(ainfo);
    mat invsubinfo(subparms.length(),subparms.length());
    for(int i=0; i<sublen; i++){ for(int j=0; j<sublen; j++){
            invsubinfo(i,j) = invinfo(i+1,j+1);
    }}
    mat subinfo = inv(invsubinfo);
    List ret = List::create( Named("loglik") = ll, Named("betahat") = betahat,
        Named("grad") = grad, Named("info") = subinfo, Named("betainfo") = betainfo,
        Named("sigmasq") = sigmasq);
    return ret;
        
}


