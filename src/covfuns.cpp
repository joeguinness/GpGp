// covariance functions
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "covfuns.h"

using namespace Rcpp;

//' Isotropic Matern covariance function
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs gives a point in R^d.
//' @param covparms A vector giving positive-valued covariance parameters
//' in the form (variance, range, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range, smoothness, nugget)
//' = \eqn{(\sigma^2,\alpha,\nu,\tau^2)}, and the covariance function is parameterized
//' as
//' \deqn{ M(x,y) = \sigma^2 2^{1-\nu}/\Gamma(\nu) ( || x - y || / \alpha )^\nu K_\nu( || x - y || / \alpha )  }
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. The reason for this choice
//' is for simpler profiling of \eqn{ \sigma^2 }.
// [[Rcpp::export]]
NumericMatrix matern_isotropic(NumericVector covparms, NumericMatrix locs ){


    int dim = locs.ncol();
    int n = locs.nrow();
    std::vector<double> cloc1(dim);
    std::vector<double> cloc2(dim);
    double cparms[3];
    for(int j = 0; j < 3; j++){
        cparms[j] = covparms(j);
    }
    double nugget = covparms( 0 )*covparms( 3 );

    NumericMatrix covmat(n,n);
    for(int i1 = 0; i1 < n; i1++){
        for(int i2 = 0; i2 <= i1; i2++){
            for(int j = 0; j < dim; j++){
                cloc1[j] = locs(i1,j);
                cloc2[j] = locs(i2,j);
            }
            covmat(i1,i2) = matern_isotropic_internal( &cloc1, &cloc2, cparms );
            if( i1 == i2 ){
                covmat(i2,i2) += nugget;
            } else {
                covmat(i2,i1) = covmat(i1,i2);
            }
        }
    }
    return covmat;
}




//' Isotropic Matern covariance function on sphere
//'
//' From a matrix of longitudes and latitudes and a vector covariance parameters of the form
//' (variance, range, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param lonlat A matrix with \code{n} rows and one column with longitudes in (-180,180)
//' and one column of latitudes in (-90,90).
//' Each row of locs describes a point on the sphere.
//' @param covparms A vector giving positive-valued covariance parameters
//' in the form (variance, range, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlat[i,]} and
//' \code{lonlat[j,]}.
//' @section Matern on Sphere Domain:
//' The function first calculates the (x,y,z) 3D coordinates, and then inputs
//' the resulting locations into \code{maternIsotropic}. This means that we construct
//' covariances on the sphere by embedding the sphere in a 3D space. There has been some
//' concern expressed in the literature that such embeddings may produce distortions.
//' The source and nature of such distortions has never been articulated,
//' and to date, no such distortions have been documented. Guinness and
//' Fuentes (2016) argue that 3D embeddings produce reasonable models for data on spheres.
// [[Rcpp::export]]
NumericMatrix matern_sphere(NumericVector covparms, NumericMatrix lonlat ){

    int n = lonlat.nrow();
    double lonrad;                                  // longitude
    double latrad;                                  // latitude
    Rcpp::NumericMatrix xyz(n, 3);
    for(int i = 0; i < n; i++){
        lonrad = 2*M_PI*lonlat(i,0)/360;
        latrad = 2*M_PI*(lonlat(i,1)+90)/360;
        xyz(i,0) = sin(latrad)*cos(lonrad);         // convert lon,lat to x,y,z
        xyz(i,1) = sin(latrad)*sin(lonrad);
        xyz(i,2) = cos(latrad);
    }

    NumericMatrix covmat = matern_isotropic( covparms, xyz );
    return covmat;
}


//' Space-time Matern-covariance function
//'
//' From a matrix of locations and times and a vector covariance parameters of the form
//' (variance, spatial range, temporal range, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locstime A matrix with \code{n} rows and d+1 columns. The first d columns
//' give a location in R^d, and the last column gives a time. Each row corresponds
//' to a space-time location.
//' @param covparms A vector giving positive-valued covariance parameters
//' in the form (variance, spatial range, temporal range, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the \code{i,j} entry
//' containing the covariance between observations at \code{locstime[i,]} and
//' \code{locstime[j,]}.
// [[Rcpp::export]]
NumericMatrix matern_space_time(NumericVector covparms, NumericMatrix locstime ){

    int n = locstime.nrow();
    int d = locstime.ncol() - 1;
    NumericVector covparms1(4);
    covparms1[0] = covparms[0];                    // variance
    covparms1[1] = 1;                              // locations scaled below, so set range = 1
    covparms1[2] = covparms[3];                 // smoothness
    covparms1[3] = covparms[4];
    for(int i = 0; i < n; i++){
        for(int j=0; j<d; j++){
            locstime(i,j) = locstime(i,j)/covparms[2];
        }
        locstime(i,d) = locstime(i,d)/covparms[3];
    }

    NumericMatrix covmat = matern_isotropic( covparms1, locstime );
    return covmat;
}


//' Isotropic Matern covariance function on sphere-time
//'
//' From a matrix of longitudes, latitudes, and times and a vector covariance parameters of the form
//' (variance, spatial range, temporal range, smoothness, nugget), return the square matrix of
//' all pairwise covariances.
//' @param lonlattime A matrix with \code{n} rows and one column with longitudes in (-180,180),
//' one column of latitudes in (-90,90), and one column of times.
//' Each row of locs describes a point on the sphere-time domain.
//' @param covparms A vector giving positive-valued covariance parameters
//' in the form (variance, spatial range, temporal range, smoothness, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{lonlattime[i,]} and
//' \code{lonlattime[j,]}.
//' @section Matern on Sphere-Time Domain:
//' The function first calculates the (x,y,z) 3D spatial coordinates, and scales the
//' spatial coordinates by the spatial range and the temporal coordinates by
//' the temporal range. Then the scaled coordinates are input into
//' \code{maternIsotropic}. This means that we construct
//' covariances on the sphere-time by embedding the sphere-time domain in a 4D space, with a
//' different range parameter for the three spatial dimensions versus the one
//' temporal dimension. There has been some
//' concern expressed in the literature that embedding points on the sphere into
//' a 3D domain may cause distortions.
//' The source and nature of such distortions has never been articulated,
//' and to date, no such distortions have been documented. Guinness and
//' Fuentes (2016) argue that 3D embeddings produce reasonable models for data on spheres.
// [[Rcpp::export]]
NumericMatrix matern_sphere_time(NumericVector covparms, NumericMatrix lonlattime ){

    int n = lonlattime.nrow();
    NumericVector covparms1(4);
    covparms1[0] = covparms[0];                    // variance
    covparms1[1] = 1;                              // locations scaled below, so set range = 1
    covparms1[2] = covparms[3];                 // smoothness
    covparms1[3] = covparms[4];
    double lonrad;
    double latrad;
    Rcpp::NumericMatrix xyzt(n, 4);
    for(int i = 0; i < n; i++){
        lonrad = 2*M_PI*lonlattime(i,0)/360;
        latrad = 2*M_PI*(lonlattime(i,1)+90)/360;
        xyzt(i,0) = sin(latrad)*cos(lonrad)/covparms[1];   // convert lon,lat,time to
        xyzt(i,1) = sin(latrad)*sin(lonrad)/covparms[1];   // scaled x,y,z, and scaled time
        xyzt(i,2) = cos(latrad)/covparms[1];
        xyzt(i,3) = lonlattime(i,2)/covparms[2];
    }

    NumericMatrix covmat = matern_isotropic( covparms1, xyzt );
    return covmat;
}
