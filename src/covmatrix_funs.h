#ifndef ARMACOVMATRIX_FUNS_H
#define ARMACOVMATRIX_FUNS_H

#include <RcppArmadillo.h>
#include <iostream>
#include <vector>

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
#include "covmatrix_funs_13.h"

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]


void get_covfun(std::string covfun_name_string,  mat (*p_covfun[1])(arma::vec, arma::mat), cube (*p_d_covfun[1])(arma::vec, arma::mat)  )
{
    
    if( covfun_name_string.compare("matern_isotropic") == 0 )
    { 
        p_covfun[0] = matern_isotropic; 
        p_d_covfun[0] = d_matern_isotropic;
    } 
    else if( covfun_name_string.compare("exponential_isotropic") == 0 )
    { 
        p_covfun[0] = exponential_isotropic; 
        p_d_covfun[0] = d_exponential_isotropic;
    }
    else if( covfun_name_string.compare("matern_anisotropic2D") == 0 )
    { 
        p_covfun[0] = matern_anisotropic2D; 
        p_d_covfun[0] = d_matern_anisotropic2D;
    } 
    else if( covfun_name_string.compare("exponential_anisotropic2D") == 0 )
    { 
        p_covfun[0] = exponential_anisotropic2D; 
        p_d_covfun[0] = d_exponential_anisotropic2D;
    } 
    else if( covfun_name_string.compare("exponential_anisotropic3D") == 0 )
    { 
        p_covfun[0] = exponential_anisotropic3D; 
        p_d_covfun[0] = d_exponential_anisotropic3D;
    } 
    else if( covfun_name_string.compare("exponential_anisotropic3D_alt") == 0 )
    { 
        p_covfun[0] = exponential_anisotropic3D_alt; 
        p_d_covfun[0] = d_exponential_anisotropic3D_alt;
    } 
    else if( covfun_name_string.compare("matern_anisotropic3D") == 0 )
    { 
        p_covfun[0] = matern_anisotropic3D; 
        p_d_covfun[0] = d_matern_anisotropic3D;
    } 
    else if( covfun_name_string.compare("matern_anisotropic3D_alt") == 0 )
    { 
        p_covfun[0] = matern_anisotropic3D_alt; 
        p_d_covfun[0] = d_matern_anisotropic3D_alt;
    } 
    else if( covfun_name_string.compare("matern15_isotropic") == 0 )
    { 
        p_covfun[0] = matern15_isotropic; 
        p_d_covfun[0] = d_matern15_isotropic;
    } 
    else if( covfun_name_string.compare("matern_spheretime") == 0 )
    { 
        p_covfun[0] = matern_spheretime; 
        p_d_covfun[0] = d_matern_spheretime;
    }
    else if( covfun_name_string.compare("exponential_spheretime") == 0 )
    { 
        p_covfun[0] = exponential_spheretime; 
        p_d_covfun[0] = d_exponential_spheretime;
    }
    else if( covfun_name_string.compare("matern_spacetime") == 0 )
    { 
        p_covfun[0] = matern_spacetime; 
        p_d_covfun[0] = d_matern_spacetime;
    }
    else if( covfun_name_string.compare("exponential_spacetime") == 0 )
    { 
        p_covfun[0] = exponential_spacetime; 
        p_d_covfun[0] = d_exponential_spacetime;
    }
    else if( covfun_name_string.compare("matern_scaledim") == 0 )
    { 
        p_covfun[0] = matern_scaledim;
        p_d_covfun[0] = d_matern_scaledim;
    }
    else if( covfun_name_string.compare("exponential_scaledim") == 0 )
    { 
        p_covfun[0] = exponential_scaledim;
        p_d_covfun[0] = d_exponential_scaledim;
    }
    else if( covfun_name_string.compare("matern_sphere") == 0 )
    { 
        p_covfun[0] = matern_sphere; 
        p_d_covfun[0] = d_matern_sphere;
    }
    else if( covfun_name_string.compare("exponential_sphere") == 0 )
    { 
        p_covfun[0] = exponential_sphere; 
        p_d_covfun[0] = d_exponential_sphere;
    }
    else if( covfun_name_string.compare("matern_sphere_warp") == 0 )
    { 
        p_covfun[0] = matern_sphere_warp; 
        p_d_covfun[0] = d_matern_sphere_warp;
    }
    else if( covfun_name_string.compare("exponential_sphere_warp") == 0 )
    { 
        p_covfun[0] = exponential_sphere_warp; 
        p_d_covfun[0] = d_exponential_sphere_warp;
    }
    else if( covfun_name_string.compare("matern_spheretime_warp") == 0 )
    { 
        p_covfun[0] = matern_spheretime_warp; 
        p_d_covfun[0] = d_matern_spheretime_warp;
    }
    else if( covfun_name_string.compare("exponential_spheretime_warp") == 0 )
    { 
        p_covfun[0] = exponential_spheretime_warp; 
        p_d_covfun[0] = d_exponential_spheretime_warp;
    }
    else if( covfun_name_string.compare("matern_nonstat_var") == 0 )
    { 
        p_covfun[0] = matern_nonstat_var; 
        p_d_covfun[0] = d_matern_nonstat_var;
    } 
    else if( covfun_name_string.compare("exponential_nonstat_var") == 0 )
    { 
        p_covfun[0] = exponential_nonstat_var; 
        p_d_covfun[0] = d_exponential_nonstat_var;
    } 
    else if( covfun_name_string.compare("matern15_scaledim") == 0 )
    { 
        p_covfun[0] = matern15_scaledim;
        p_d_covfun[0] = d_matern15_scaledim;
    }
    else if( covfun_name_string.compare("matern25_isotropic") == 0 )
    { 
        p_covfun[0] = matern25_isotropic; 
        p_d_covfun[0] = d_matern25_isotropic;
    } 
    else if( covfun_name_string.compare("matern35_isotropic") == 0 )
    { 
        p_covfun[0] = matern35_isotropic; 
        p_d_covfun[0] = d_matern35_isotropic;
    } 
    else if( covfun_name_string.compare("matern45_isotropic") == 0 )
    { 
        p_covfun[0] = matern45_isotropic; 
        p_d_covfun[0] = d_matern45_isotropic;
    } 
    else if( covfun_name_string.compare("matern25_scaledim") == 0 )
    { 
        p_covfun[0] = matern25_scaledim;
        p_d_covfun[0] = d_matern25_scaledim;
    }
    else if( covfun_name_string.compare("matern35_scaledim") == 0 )
    { 
        p_covfun[0] = matern35_scaledim;
        p_d_covfun[0] = d_matern35_scaledim;
    }
    else if( covfun_name_string.compare("matern45_scaledim") == 0 )
    { 
        p_covfun[0] = matern45_scaledim;
        p_d_covfun[0] = d_matern45_scaledim;
    }
    else if( covfun_name_string.compare("matern_categorical") == 0 )
    { 
        p_covfun[0] = matern_categorical;
        p_d_covfun[0] = d_matern_categorical;
    }
    else if( covfun_name_string.compare("matern_spacetime_categorical") == 0 )
    { 
        p_covfun[0] = matern_spacetime_categorical;
        p_d_covfun[0] = d_matern_spacetime_categorical;
    }
    else if( covfun_name_string.compare("matern_spacetime_categorical_local") == 0 )
    { 
        p_covfun[0] = matern_spacetime_categorical_local;
        p_d_covfun[0] = d_matern_spacetime_categorical_local;
    }
    else { // stop the program
        Rcpp::Rcout << "Unrecognized Covariance Function Name \n";
    }

}


#endif
