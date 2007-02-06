/*** Mvtdist3.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  13/09/2006  as Mvtdist2.cpp
//    CHANGES:  11/12/2006: namespace Mvtdist3 created
//                          some parts of the Wishart sampling changed
//
//    rwishartEye3:  11/12/2006 from rwishartEye2 (21/09/2006)
//       rwishart3:  11/12/2006 from rwishart2    (21/09/2006)
//
//             rmvnorm2006:  06/12/2006
//         rmvnormZero2006:  30/01/2007
//            rmvnormQ2006:  06/12/2006
//        rmvnormQZero2006:  30/01/2007
//            rmvnormC2006:  06/12/2006
//            rmvnormR2006:  06/12/2006
//
// PURPOSE: Utilities for several multivariate distributions
//          
//     * partially taken from mvtdist.cpp and mvtdist.h of the bayesSurv package
//
//
/* ********************************************************************************* */

#include "Mvtdist3.h"

namespace Mvtdist3 {

/* ********************************************************************************* */
/* rwishartEye: Sample from the Wishart distribution W(nu, Eye),                     */
/*              where Eye stands for an identity matrix                              */
/*                                                                                   */
/* Bartlett's decomposition as described on p. 99 of                                 */
/* Ripley (1987). Stochastic simulation. New York: Wiley is used.                    */
/*                                                                                   */
/* See also p. 41 of my red notes.                                                   */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
/*                                                                                   */
/*    work:  working array of length LT(dim)                                         */
/*      nu:  degrees of freedom of the Wishart distribution                          */
/*     dim:  dimension of the Wishart distribution                                   */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
void
rwishartEye3(double *W, double *work, const double *nu, const int *dim)
{
  static int i, j, k;
  static double *V, *epsilon, *epsilon2, *epsilonBeg;
  static double shape, first_elem;
  static const double scale = 2.0;

    /*** V matrix (random sample from Wishart(nu, eye) ***/
    /* ================================================= */
  V = W;
  epsilonBeg = work;
  epsilon = epsilonBeg;

    /* 0th column */
  shape = *nu/2;                      /* REMEMBER: Indeces go from 0 in C++ */
  *V = rgamma(shape, scale);          /* V[0,0] = zeta[0]                   */
  *epsilon = sqrt(*V);
  first_elem = *epsilon;
  V++;   
  epsilon++;

  for (i = 1; i < *dim; i++){
    *epsilon = norm_rand();
    *V = *epsilon * first_elem;
    V++;
    epsilon++;
  }

    /* 1st, ..., (dim-1)th column */
  for (j = 1; j < *dim; j++){
    shape = (*nu - j)/2;                /* REMEMBER: Indeces go from 0 in C++ */
    *V = rgamma(shape, scale);          /* V[j,j] = zeta[j]                   */
    *epsilon = sqrt(*V);
    first_elem = *epsilon;
    V++;
    epsilon++;

    for (i = j+1; i < *dim; i++){
      *epsilon = rnorm(0, 1);
      *V = *epsilon * first_elem;
      V++;
      epsilon++;
    }

    epsilon2 = epsilonBeg +j;   /* go to epsilon[j,0]                     */
    for (k = 0; k < j; k++){
      V -= *dim - j;    /* go back to the beginning of the column */
      first_elem = *epsilon2;
      for (i = j; i < *dim; i++){
        *V += *epsilon2 * first_elem;
        V++;
        epsilon2++;
      }        
      epsilon2 += j - k - 1;   /* go to epsilon[j,k+1]  */
    }
  }
 
  return;
}


/* ********************************************************************************* */
/* rwishart2, rwishart2a: Sample from the Wishart distribution W(nu, S)              */
/*                        (parametrization as in Gelman, 2005)                       */
/*                                                                                   */
/* Algorithm from Ripley (1987), p. 99 is used.                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
/*                                                                                   */
/* Different versions are provided with respect to the specification of S^{-1}       */
/*            Let S^{-1} = L * t(L)                                                  */
/*                     S = t(L^{-1}) * L^{-1}                                        */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
/*                                                                                   */
/*    work:  working array of length 2*dim^2                                         */
/*      nu:  degrees of freedom of the Wishart distribution                          */
/*                                                                                   */
/* rwishart3:                                                                        */
/*      invS:  lower triangle of the matrix dim x dim with S^{-1}                    */
/*             * will be decomposed by this procedure                                */
/*             OR matrix L, where S^{-1} = L*t(L)                                    */
/*                                                                                   */
/* must_decomp:  1 if invS must be decomposed, 0 if invS is already decomposed       */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
void
rwishart3(double *W,  double *work,  const double *nu,  double *invS,  const int *dim,  const int &must_decomp)
{
  static int info[1];
  static double shape, scale;
  static int dim2;
  static double *tW1;

  if (*dim == 1){  
    /*** Univariate Wishart(nu, S) = gamma(shape=nu/2, rate=1/(2S)) = gamma(shape=nu/2, scale=2S) ***/
    shape = *nu/2;
    if (must_decomp){
      scale = 2/(*invS); 
    }
    else{
      scale = 2/((*invS)*(*invS)); 
    }
    *W = rgamma(shape, scale);
  }
  else{
    /*** Sample V ~ Wishart(nu, eye), store it in W ***/
    Mvtdist3::rwishartEye3(W, work, nu, dim);     

    /*** Decomposition of invS matrix: invS = L*t(L) ***/
    if (must_decomp){
      AK_BLAS_LAPACK::chol_dpptrf(invS, dim, info);
      if (*info) throw returnR("Mvtdist3.cpp: rwishart3(...) error. Scale matrix is not PD.", 1);
    }

    /*** Compute W = t(L)^{-1}*V*L^{-1}                                   ***/
    /*** a) copy V stored as LT to work, where in work, V is fully stored ***/
    AK_BLAS_LAPACK::LT2Rect(work, W, *dim);

    /** b) compute W1 = t(L)^{-1}*V,  that is solve t(L)*W1 = V  **/
    AK_BLAS_LAPACK::chol_solve_backward_system(work, invS, dim, dim);

    /** c) compute W = W1*L^{-1},  that is W*L = W1,  that is t(W1) = t(L)*t(W) **/
    /**    W is symmetric, so that t(W) = W                                     **/
    /**    so that, we have to solve t(L)*W = t(W1)                             **/
    dim2 = (*dim)*(*dim);
    tW1  = work + dim2;
    AK_BLAS_LAPACK::transposition(tW1, work, dim, dim);    
    AK_BLAS_LAPACK::chol_solve_backward_system(tW1, invS, dim, dim);

    /** Copy tW1 to W  **/
    AK_BLAS_LAPACK::Rect2LT(W, tW1, *dim);
  }

  return;
}


/***** R wrapper to rwishart3                                                           *****/
/**                                                                                        **/
/** ====================================================================================== **/
//
// x[LT(dim)*nrandom]:   OUTPUT: sampled values (lower triangles only)
// work[dim^2]:          working array
// nu[1]:                Wishart degrees of freedom
// invS[LT(dim)]:        inverse scale matrix of the Wishart distribution
// dim:                  dimension of the Wishart distribution
// nrandom:              number of sampled points
//
extern "C"{
  void
  rwishartR3(double *W,  double *work,  const double *nu,  double *invS,  const int *dim,  const int *nrandom)
  {
    try{
      GetRNGstate();

      double *WP = W;
      int lWP    = ((*dim)*(*dim+1))/2;
      int info[1];

      /*** Decomposition of invS matrix: invS = L*t(L) ***/
      AK_BLAS_LAPACK::chol_dpptrf(invS, dim, info);
      if (*info) throw returnR("Mvtdist3.cpp: rwishartR3(...) error. Scale matrix is not PD.", 1);

      for (int i = 0; i < *nrandom; i++){
	Mvtdist3::rwishart3(WP, work, nu, invS, dim, 0);
        WP += lWP;
      }

      PutRNGstate();
      return;
    }
    catch (returnR rr){
      PutRNGstate();
      return;
    }    
  }  
}


/***** rmvnorm2006:  Sample from N(mu, Sigma)                                      *****/
/**                                                                                   **/
/** Algorithm 2.3 on page 34 of Rue and Held (2005) is used                           **/
/**                                                                                   **/
/** ================================================================================= **/
//
// x[nx]:      OUTPUT: sampled value
// mu[nx]:     mean of the normal distribution
// L[LT(nx)]:  lower triangle of the Cholesky decomposition of matrix Sigma,
//             that is, Sigma = L*t(L)
// nx:         dimension of the normal distribution
//
void
rmvnorm2006(double *x,  const double *mu,  const double *L,  const int *nx)
{
  static int i;
  static double *xP;
  static const double *muP;

  /** Sample z ~ N(0, I) **/
  xP = x;
  for (i = 0; i < *nx; i++){
    *xP = norm_rand();
    xP++;
  }

  /** Compute v = L*z **/
  AK_BLAS_LAPACK::a_La(x, L, nx);

  /** Compute x = mu + v **/
  xP = x;
  muP = mu;
  for (i = 0; i < *nx; i++){
    *xP += *muP;
    xP++;
    muP++;
  }

  return;
}
 

/***** rmvnormZero2006:  Sample from N(0, Sigma)                                   *****/
/**                                                                                   **/
/** Algorithm 2.3 on page 34 of Rue and Held (2005) is used                           **/
/**                                                                                   **/
/** ================================================================================= **/
//
// x[nx]:      OUTPUT: sampled value
// L[LT(nx)]:  lower triangle of the Cholesky decomposition of matrix Sigma,
//             that is, Sigma = L*t(L)
// nx:         dimension of the normal distribution
//
void
rmvnormZero2006(double *x,  const double *L,  const int *nx)
{
  static int i;
  static double *xP;

  /** Sample z ~ N(0, I) **/
  xP = x;
  for (i = 0; i < *nx; i++){
    *xP = norm_rand();
    xP++;
  }

  /** Compute v = L*z **/
  AK_BLAS_LAPACK::a_La(x, L, nx);

  return;
}


/***** rmvnormQ2006:  Sample from N(mu, Q^{-1})                                    *****/
/**                                                                                   **/
/** Algorithm 2.4 on page 34 of Rue and Held (2005) is used                           **/
/**                                                                                   **/
/** ================================================================================= **/
//
// x[nx]:      OUTPUT: sampled value
// mu[nx]:     mean of the normal distribution
// L[LT(nx)]:  lower triangle of the Cholesky decomposition of matrix Q,
//             that is, Q = L*t(L)
// nx:         dimension of the normal distribution
//
void
rmvnormQ2006(double *x,  const double *mu,  const double *L,  const int *nx)
{
  static int i;
  static double *xP;
  static const double *muP;

  /** Sample z ~ N(0, I) **/
  xP = x;
  for (i = 0; i < *nx; i++){
    *xP = norm_rand();
    xP++;
  }

  /** Solve t(L)*v = z **/
  AK_BLAS_LAPACK::chol_solve_backward(x, L, nx);

  /** Compute x = mu + v **/
  xP = x;
  muP = mu;
  for (i = 0; i < *nx; i++){
    *xP += *muP;
    xP++;
    muP++;
  }

  return;
}


/***** rmvnormQZero2006:  Sample from N(0, Q^{-1})                                 *****/
/**                                                                                   **/
/** Algorithm 2.4 on page 34 of Rue and Held (2005) is used                           **/
/**                                                                                   **/
/** ================================================================================= **/
//
// x[nx]:      OUTPUT: sampled value
// L[LT(nx)]:  lower triangle of the Cholesky decomposition of matrix Q,
//             that is, Q = L*t(L)
// nx:         dimension of the normal distribution
//
void
rmvnormQZero2006(double *x,  const double *L,  const int *nx)
{
  static int i;
  static double *xP;

  /** Sample z ~ N(0, I) **/
  xP = x;
  for (i = 0; i < *nx; i++){
    *xP = norm_rand();
    xP++;
  }

  /** Solve t(L)*v = z **/
  AK_BLAS_LAPACK::chol_solve_backward(x, L, nx);

  return;
}


/***** rmvnormC2006:  Sample from N_C(b, Q) = N(Q^{-1}*b, Q^{-1})                  *****/
/**                                                                                   **/
/** Algorithm 2.5 on page 35 of Rue and Held (2005) is used                           **/
/**                                                                                   **/
/** ================================================================================= **/
//
// x[nx]:      OUTPUT: sampled value
// b[nx]:      INPUT:  canonical mean of the normal distribution
//             OUTPUT: mean of the normal distribution
// L[LT(nx)]:  lower triangle of the Cholesky decomposition of matrix Q,
//             that is, Q = L*t(L)
// nx:         dimension of the normal distribution
//
void
rmvnormC2006(double *x,  double *b,  const double *L,  const int *nx)
{
  static int i;
  static double *xP;
  static const double *bP;

  /** Solve L*w = b  **/
  AK_BLAS_LAPACK::chol_solve_forward(b, L, nx);

  /** Solve t(L)*mu = w **/
  AK_BLAS_LAPACK::chol_solve_backward(b, L, nx);

  /** Sample z ~ N(0, I) **/
  xP = x;
  for (i = 0; i < *nx; i++){
    *xP = norm_rand();
    xP++;
  }

  /** Solve t(L)*v = z **/
  AK_BLAS_LAPACK::chol_solve_backward(x, L, nx);

  /** Compute x = mu + v **/
  xP = x;
  bP = b;
  for (i = 0; i < *nx; i++){
    *xP += *bP;
    xP++;
    bP++;
  }

  return;
}


/***** R wrapper to rmvnorm2006, rmvnormQ2006, rmvnormC2006                             *****/
/**                                                                                        **/
/** ====================================================================================== **/
//
// x[nx*nrandom]:   OUTPUT: sampled values
// mub[nx]:         mean or canonical mean of the normal distribution
// QS[LT(nx)]:      lower triangle of either Sigma or Q = Sigma^{-1}
// err[1]:          error flag
// nx:              dimension of the normal distribution
// nrandom:         number of random numbers to be generated
// version:         0 = wrapper to rmvnorm2006
//                  1 = wrapper to rmvnormQ2006
//                  2 = wrapper to rmvnormC2006
//
extern "C"{
  void
  rmvnormR2006(double *x,  double *mub,  double *QS,  int *err,  const int *nx,  const int *nrandom,  const int *version)
  {
    try{
      GetRNGstate();

      double *xP = x;

      /** Cholesky decomposition of QS **/
      AK_BLAS_LAPACK::chol_dpptrf(QS, nx, err);
      if (*err){
        throw returnR("Error in Mvtdist3.cpp: rmvnormR2006. Supplied covariance/precision matrix is not positive definite", 1);
      }

      switch (*version){
      case 0:
        for (int i = 0; i < *nrandom; i++){
	  Mvtdist3::rmvnorm2006(xP, mub, QS, nx);
          xP += *nx;
        }
        break;

      case 1:
        for (int i = 0; i < *nrandom; i++){
          Mvtdist3::rmvnormQ2006(xP, mub, QS, nx);
          xP += *nx;
        }
        break;

      case 2:
        Mvtdist3::rmvnormC2006(xP, mub, QS, nx);
        xP += *nx;
        for (int i = 1; i < *nrandom; i++){
          Mvtdist3::rmvnormQ2006(xP, mub, QS, nx);
          xP += *nx;
        }
        break;        

      default:
        throw returnR("Error in Mvtdist3.cpp: rmvnormR2006. Unknown value of the argument version", 1);
      }

      PutRNGstate();
    }
    catch (returnR rr){
      PutRNGstate();
      return;
    }    
  } 
}

}  /*** end of namespace Mvtdist3 ***/

