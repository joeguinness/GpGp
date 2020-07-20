#ifndef ARMACOVMATRIX_FUNS_H
#define ARMACOVMATRIX_FUNS_H



#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include "basis.h"

// covariance functions
#include "covmatrix_funs_01.h"
#include "covmatrix_funs_02.h"
#include "covmatrix_funs_03.h"
#include "covmatrix_funs_04.h"
#include "covmatrix_funs_05.h"
#include "covmatrix_funs_06.h"
#include "covmatrix_funs_07.h"
#include "covmatrix_funs_08.h"
#include "covmatrix_funs_09.h"
#include "covmatrix_funs_10.h"
#include "covmatrix_funs_11.h"
#include "covmatrix_funs_12.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]


struct covfun_t {
    arma::mat    (*p_covfun)(arma::vec, arma::mat);
    arma::cube (*p_d_covfun)(arma::vec, arma::mat);
} ;


covfun_t get_covfun(std::string covfun_name_string)
{
    
    covfun_t covstruct;
    if( covfun_name_string.compare("matern_isotropic") == 0 )
    { 
        covstruct.p_covfun = &matern_isotropic; 
        covstruct.p_d_covfun = &d_matern_isotropic;
    } 
    else if( covfun_name_string.compare("matern15_isotropic") == 0 )
    { 
        covstruct.p_covfun = &matern15_isotropic; 
        covstruct.p_d_covfun = &d_matern15_isotropic;
    } 
    else if( covfun_name_string.compare("matern25_isotropic") == 0 )
    { 
        covstruct.p_covfun = &matern25_isotropic; 
        covstruct.p_d_covfun = &d_matern25_isotropic;
    } 
    else if( covfun_name_string.compare("matern35_isotropic") == 0 )
    { 
        covstruct.p_covfun = &matern35_isotropic; 
        covstruct.p_d_covfun = &d_matern35_isotropic;
    } 
    else if( covfun_name_string.compare("matern45_isotropic") == 0 )
    { 
        covstruct.p_covfun = &matern45_isotropic; 
        covstruct.p_d_covfun = &d_matern45_isotropic;
    } 
    else if( covfun_name_string.compare("matern_anisotropic2D") == 0 )
    { 
        covstruct.p_covfun = &matern_anisotropic2D; 
        covstruct.p_d_covfun = &d_matern_anisotropic2D;
    } 
    else if( covfun_name_string.compare("exponential_anisotropic2D") == 0 )
    { 
        covstruct.p_covfun = &exponential_anisotropic2D; 
        covstruct.p_d_covfun = &d_exponential_anisotropic2D;
    } 
    else if( covfun_name_string.compare("exponential_anisotropic3D") == 0 )
    { 
        covstruct.p_covfun = &exponential_anisotropic3D; 
        covstruct.p_d_covfun = &d_exponential_anisotropic3D;
    } 
    else if( covfun_name_string.compare("matern_anisotropic3D") == 0 )
    { 
        covstruct.p_covfun = &matern_anisotropic3D; 
        covstruct.p_d_covfun = &d_matern_anisotropic3D;
    } 
    else if( covfun_name_string.compare("matern_nonstat_var") == 0 )
    { 
        covstruct.p_covfun = &matern_nonstat_var; 
        covstruct.p_d_covfun = &d_matern_nonstat_var;
    } 
    else if( covfun_name_string.compare("exponential_nonstat_var") == 0 )
    { 
        covstruct.p_covfun = &exponential_nonstat_var; 
        covstruct.p_d_covfun = &d_exponential_nonstat_var;
    } 
    else if( covfun_name_string.compare("exponential_isotropic") == 0 )
    { 
        covstruct.p_covfun = &exponential_isotropic; 
        covstruct.p_d_covfun = &d_exponential_isotropic;
    }
    else if( covfun_name_string.compare("matern_sphere") == 0 )
    { 
        covstruct.p_covfun = &matern_sphere; 
        covstruct.p_d_covfun = &d_matern_sphere;
    }
    else if( covfun_name_string.compare("exponential_sphere") == 0 )
    { 
        covstruct.p_covfun = &exponential_sphere; 
        covstruct.p_d_covfun = &d_exponential_sphere;
    }
    else if( covfun_name_string.compare("matern_sphere_warp") == 0 )
    { 
        covstruct.p_covfun = &matern_sphere_warp; 
        covstruct.p_d_covfun = &d_matern_sphere_warp;
    }
    else if( covfun_name_string.compare("exponential_sphere_warp") == 0 )
    { 
        covstruct.p_covfun = &exponential_sphere_warp; 
        covstruct.p_d_covfun = &d_exponential_sphere_warp;
    }
    else if( covfun_name_string.compare("matern_spheretime_warp") == 0 )
    { 
        covstruct.p_covfun = &matern_spheretime_warp; 
        covstruct.p_d_covfun = &d_matern_spheretime_warp;
    }
    else if( covfun_name_string.compare("exponential_spheretime_warp") == 0 )
    { 
        covstruct.p_covfun = &exponential_spheretime_warp; 
        covstruct.p_d_covfun = &d_exponential_spheretime_warp;
    }
    else if( covfun_name_string.compare("matern_spheretime") == 0 )
    { 
        covstruct.p_covfun = &matern_spheretime; 
        covstruct.p_d_covfun = &d_matern_spheretime;
    }
    else if( covfun_name_string.compare("exponential_spheretime") == 0 )
    { 
        covstruct.p_covfun = &exponential_spheretime; 
        covstruct.p_d_covfun = &d_exponential_spheretime;
    }
    else if( covfun_name_string.compare("matern_spacetime") == 0 )
    { 
        covstruct.p_covfun = &matern_spacetime; 
        covstruct.p_d_covfun = &d_matern_spacetime;
    }
    else if( covfun_name_string.compare("exponential_spacetime") == 0 )
    { 
        covstruct.p_covfun = &exponential_spacetime; 
        covstruct.p_d_covfun = &d_exponential_spacetime;
    }
    else if( covfun_name_string.compare("matern_scaledim") == 0 )
    { 
        covstruct.p_covfun = &matern_scaledim;
        covstruct.p_d_covfun = &d_matern_scaledim;
    }
    else if( covfun_name_string.compare("matern15_scaledim") == 0 )
    { 
        covstruct.p_covfun = &matern15_scaledim;
        covstruct.p_d_covfun = &d_matern15_scaledim;
    }
    else if( covfun_name_string.compare("matern25_scaledim") == 0 )
    { 
        covstruct.p_covfun = &matern25_scaledim;
        covstruct.p_d_covfun = &d_matern25_scaledim;
    }
    else if( covfun_name_string.compare("matern35_scaledim") == 0 )
    { 
        covstruct.p_covfun = &matern35_scaledim;
        covstruct.p_d_covfun = &d_matern35_scaledim;
    }
    else if( covfun_name_string.compare("matern45_scaledim") == 0 )
    { 
        covstruct.p_covfun = &matern45_scaledim;
        covstruct.p_d_covfun = &d_matern45_scaledim;
    }
    else if( covfun_name_string.compare("exponential_scaledim") == 0 )
    { 
        covstruct.p_covfun = &exponential_scaledim;
        covstruct.p_d_covfun = &d_exponential_scaledim;
    }
    else { // stop the program
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
    }
    return covstruct;
}


#endif
