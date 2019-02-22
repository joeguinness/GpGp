#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "covfuns.h"
#include "vecchia_fun_grouped.h"

using namespace std;
using namespace Rcpp;


//' Grouped Vecchia's approximation to the Gaussian loglikelihood
//' 
//' This function returns the grouped version (Guinness, 2018) of 
//' Vecchia's (1988) approximation to the Gaussian
//' loglikelihood. The approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @param NNlist A list with grouped neighbor information. 
//' Usually the output from \code{group_obs(NNarray)}.
//' @inheritParams vecchia_loglik
//' @return grouped version of Vecchia's approximation to the Gaussian loglikelihood
//' @examples
//' n1 <- 40
//' n2 <- 40
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' covparms <- c(2, 0.2, 0.75, 0)
//' y <- fast_Gp_sim(covparms, "matern_isotropic", locs, 50 ) 
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' NNlist <- group_obs(NNarray)
//' loglik <- vecchia_loglik_grouped( covparms, "matern_isotropic", y, locs, NNlist )
//' @export
// [[Rcpp::export]]
NumericVector vecchia_loglik_grouped(NumericVector covparms, StringVector covfun_name,
                                  NumericVector y,
                                  NumericMatrix locs, List NNlist) {

    NumericVector ll(1);        // loglikelihood to be returned
    NumericVector Linv(1);    // Linv not to be returned
    vecchia_grouped(covparms, covfun_name, locs, NNlist, y, &Linv, &ll, 1);
    return ll;
}



//' Inverse Cholesky factor implied by Vecchia's approximation
//' 
//' This function returns the entries of the sparse approximation to 
//' the Cholesky factor implied by Vecchia's (1988) approximation to the Gaussian
//' loglikelihood. The approximation modifies the ordered conditional
//' specification of the joint density; rather than each term in the product
//' conditioning on all previous observations, each term conditions on
//' a small subset of previous observations.
//' @inheritParams vecchia_loglik
//' @inheritParams vecchia_loglik_grouped
//' @return vector containing entries of grouped approximation to inverse Cholesky
//' @examples
//' n1 <- 40
//' n2 <- 40
//' n <- n1*n2
//' locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
//' covparms <- c(2, 0.2, 0.75, 0)
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' Linv <- vecchia_Linv( covparms, "matern_isotropic", locs, NNarray )
//' @export
// [[Rcpp::export]]
NumericVector vecchia_Linv_grouped(NumericVector covparms, StringVector covfun_name,
                            NumericMatrix locs, List NNlist) {

    NumericVector y(locs.nrow());
    NumericVector ll(1);        // loglikelihood not to be returned
    NumericVector local_resp_inds = NNlist["local_resp_inds"];
    int nentries = 0;
    for(int j=0; j<local_resp_inds.length(); j++){
        nentries += local_resp_inds[j];
    }
    NumericVector Linv(nentries);    // Linv to be returned
    vecchia_grouped(covparms, covfun_name, locs, NNlist, y, &Linv, &ll, 2);
    return Linv;
}



//' Multiply approximate inverse Cholesky by a vector
//' 
//' Vecchia's approximation implies a sparse approximation to the 
//' inverse Cholesky factor of the covariance matrix. This function
//' returns the result of multiplying that matrix by a vector.
//' @param Linv Entries of the sparse inverse Cholesky factor,
//' usually the output from \code{vecchiaLinv}.
//' @param z the vector to be multiplied
//' @inheritParams vecchia_loglik_grouped
//' @return the product of the sparse inverse Cholesky factor with a vector
//' @examples
//' n <- 2000
//' locs <- matrix( runif(2*n), n, 2 )
//' covparms <- c(2, 0.2, 0.75, 0.1)
//' ord <- order_maxmin(locs)
//' NNarray <- find_ordered_nn(locs,20)
//' Linv <- vecchia_Linv( covparms, "matern_isotropic", locs, NNarray )
//' z1 <- rnorm(n)
//' y <- fast_Gp_sim_Linv(Linv,NNarray,z1)
//' z2 <- Linv_mult(Linv, y, NNarray)
//' print( sum( (z1-z2)^2 ) )
//' @export
// [[Rcpp::export]]
NumericVector Linv_mult_grouped(NumericVector Linv, NumericVector z,
                                  List NNlist) {

    // return x = Linv * z
    int i, j, k, nresp;

    int n = z.length();
    NumericVector x(n);

    // vector of all indices
    std::vector<int> all_inds = Rcpp::as<std::vector<int> >(NNlist["all_inds"]);
    // vector of local response indices
    std::vector<int> local_resp_inds = Rcpp::as<std::vector<int> >(NNlist["local_resp_inds"]);
    // vector of global response indices
    std::vector<int> global_resp_inds = Rcpp::as<std::vector<int> >(NNlist["global_resp_inds"]);
    // last index of each block in all_inds
    std::vector<int> last_ind_of_block = Rcpp::as<std::vector<int> >(NNlist["last_ind_of_block"]);
    // last response index of each block in local_resp_inds and global_resp_inds
    std::vector<int> last_resp_of_block = Rcpp::as<std::vector<int> >(NNlist["last_resp_of_block"]);

    int nb = last_ind_of_block.size();    
    int first_ind_block = 0;
    int first_L_ind = 0;
    int first_resp = 0;
    int cur_resp = 0;
    // rows 1 though n
    for(i=0; i<nb; i++){
        
        if(i==0){ first_resp = 0; } else { first_resp = last_resp_of_block[i-1]; }
        nresp = last_resp_of_block[i] - first_resp;
        
        for(j=0; j<nresp; j++){
            
            for(k=0; k<local_resp_inds[cur_resp]; k++){
                x( global_resp_inds[cur_resp] - 1 ) += 
                    Linv( first_L_ind + k )*z( all_inds[ first_ind_block + k ] - 1 );
            }

            // update first_L_ind
            first_L_ind += local_resp_inds[cur_resp];
            cur_resp++;
        }
        // update first ind block
        first_ind_block = last_ind_of_block[i];
    }
    return x;
}
