// ************************************************************************************* //
// ********* Gspline_update_sigma.cpp ************************************************** //
// ************************************************************************************* //
// Member functions of class Gspline to update basis standard deviation
//   and some related functions
//   * for the first specification where knots depended on the basis standard deviation
//
// 03/12/2004: 'Gspline::update_sigma'
//             'Gspline::full_sigma_pars'
//             'full_sigma_logdens0'
//             'full_sigma_logdens3'
// 10/12/2004: overall intercept and scale of the G-spline (second specification) allowed
//
#include "Gspline.h"


// ***** Gspline::update_sigma: Update basis standard deviations
// =========================================================================================
void
Gspline::update_sigma(const double* regresResM,  const int* rM,  const int* nP,  const int* iter)
{
  static int j, j_, iter_nr, err_nr, overrelax, itmp;
  static double pars[4*_max_dim];
  static int ipars[1];
  static double slice[3];
  static double gx[2];
  static double tmp, gx0, horiz, dgx, newinvsigma2;

  const double* zeta_1;
  const double* sqrt_eta;
  const double* xi_half;
  
  /** Compute parameters of the full conditional distributions for all sigma's **/
  full_sigma_pars(pars, regresResM, rM, nP);

  /** Loop over dimensions **/
  zeta_1    = pars;
  sqrt_eta  = zeta_1 + 1;
  xi_half   = sqrt_eta + 1;
  for (j = 0; j < _dim; j++){
    if (_prior_for_sigma[j] == Fixed_) continue;
    overrelax = 1*((*iter/_k_overrelax_sigma[j]) != 0);

    /** Evaluate the full conditional distribution in the current point and sample the level defining the slice **/
    ipars[0] = (_prior_for_sigma[j] == SDUnif ? 1 : 0);
    full_sigma_logdens0(_invsigma2 + j, &gx0, zeta_1, ipars);
    horiz = gx0 - rexp(1);    
    
    /** Compute mode of the full conditional distribution, store it in slice[2] **/    
    if ((*zeta_1) <= 0) throw returnR("Zeta parameter for sigma update <= 1, is your sample size > 2?", 1);
    tmp = -(*xi_half) + sqrt((*xi_half)*(*xi_half) + 4*(*sqrt_eta)*(*sqrt_eta)*(*zeta_1));
    if (tmp < _epsilon) throw returnR("Trap in Gspline::update_sigma: Cannot allocate mode of the full conditional", 1);
    slice[2] = (4*(*zeta_1)*(*zeta_1))/(tmp*tmp);

    /** Give initial guesses for the interval defining the slice,                   **/
    /** starting left limit x0 of the interval must be such that g(x0) < horiz,     **/
    /** otherwise the Newton_raphson solver fails due to going below zero           **/
    /** Set initial right limit either to _invsigma2[j] or to mode + 2*std. dev.    **/
    /** but always to the right from the mode                                       **/
    if (_invsigma2[j] < slice[2]){
      dgx = (*zeta_1)/(slice[2]*slice[2]) + (*xi_half)/(2*slice[2]*sqrt(slice[2]));  /* l"(x_mode) */
      slice[1] = slice[2] + 2/sqrt(dgx);
    }
    else{
      slice[1] = _invsigma2[j];
    }
    slice[0] = 0.5*slice[2];
    full_sigma_logdens0(slice, gx, zeta_1, ipars);
    while (gx[0] >= horiz && slice[0] > _zero_invvariance){
      slice[0] *= 0.5;
      full_sigma_logdens0(slice, gx, zeta_1, ipars);
    }

    /* Compute the interval defining the slice */
    itmp = (slice[0] <= _zero_invvariance ? 1 : 0);
    for (j_=1; j_>=0; j_--){
      full_sigma_logdens3(slice+j_, gx+j_, &dgx, &tmp, zeta_1, ipars, 3);
      solver_newton_raphson(slice+j_, gx+j_, &dgx, &horiz, zeta_1, ipars, full_sigma_logdens3, 
                            &iter_nr, &_maxiter_solver_nr, &_toler_solver_nr, &_epsilon, &err_nr);
      if (err_nr >= 3){
        REprintf("\nerr_nr = %d\n", err_nr);
        REprintf("sigma[%d] = %f,  invsigma2[%d] = %f\n", j, _sigma[j], j, _invsigma2[j]);
        REprintf("mode = %f, horizontal = %f\n", slice[2], horiz);
        REprintf("zeta-1 = %f,  sqrt(eta) = %f, xi/2 = %f\n", *zeta_1, *sqrt_eta, *xi_half);
        throw returnR("Trap in Gspline::update_sigma: Unable to find an interval defining the slice", 1);
      }
    }
    if (ipars[0]){
      if (slice[0] <= zeta_1[3]) slice[0] = zeta_1[3];     /* left limit of the interval below the truncation limit       */
      if (slice[1] <= zeta_1[3]){                          /* also right limit of the interval below the truncation limit */
        _invsigma2[j] = zeta_1[3];
        _sigma[j] = 1/sqrt(_invsigma2[j]);
        _delta[j] = _c4delta[j] * _sigma[j];
        continue;
      }
    }

    /* Sample the new point */
    if (overrelax){
      Slice_sampler::ss_exact_overrelax(&newinvsigma2, slice, _invsigma2+j, &horiz, full_sigma_logdens0, zeta_1, ipars);
    }
    else{
      Slice_sampler::ss_exact_sample(&newinvsigma2, slice, gx, _invsigma2+j, &horiz, full_sigma_logdens0, zeta_1, ipars);
    }

    /* Update necessary quantities */
    _invsigma2[j] = newinvsigma2;
    _sigma[j]     = 1/sqrt(_invsigma2[j]);
    _delta[j]     = _c4delta[j] * _sigma[j];
    
    zeta_1    = xi_half + 2;
    sqrt_eta  = zeta_1 + 1;
    xi_half   = sqrt_eta + 1;
  }
  
  return;
}


// ***** Gspline::full_sigma_pars: Compute parameters of the full conditional distribution
//                                 of all sigma^{-2}
// =========================================================================================
//
// pars[4*dim] .............. computed parameters of the full conditional distribution for each sigma^{-2}
//          * raw parameters of the full conditional distribution are:
//                          zeta = zeta_0 + N/2		     
//		  	    eta  = eta_0 + 0.5*sum(y-gamma)^2
//			    xi   = c*sum[r*(y-gamma)]        
//          * we return:
//                        pars[0*j*4] ..... zeta - 1
//                        pars[1*j*4] ..... sqrt(eta)
//                        pars[2*j*4] ..... xi/2
//                        pars[3+j*4] ..... possibly truncation point 1/S^2
//
// regresResM[nP*_dim] ...... * in the case of i.i.d data, these are observations y_1,...,y_n,
//                            * if dim(y_i) = d > 1 then it is assumed that they are stored as y_1[1],...,y_1[d], ..., y_n[1], ..., y_n[d],
//                              i.e. y_i[j] = regresResM[i*d + j]
//                            * in the case of regression, these are y_i - beta'x_i - b_i'z_i
//                              - i.e. these residuals are not divided by the scale
// rM[nP] ................... component labels taking values 0, ..., total_length-1
//                            * sorted in the same way as regresResM
// nP ....................... total number of 'independent' observational vectors
//                   
void
Gspline::full_sigma_pars(double* pars,  const double* regresResM,  const int* rM,  const int* nP) const
{
  static int npars = 4;
  const double* y_obs;
  const int* rp;
  static int j, obs;
  static double y_gamma;

  /*** Perform this peace of code only when called the first time ***/
  static bool allFixed = true;
  static int jj = 0;
  while (allFixed && jj < _dim){
    if (_prior_for_sigma[j] != Fixed_) allFixed = false;
    jj++;
  }

  if (allFixed) return;

  /*** First store zeta, eta, xi in pars[j*npars+(0:2)]  ***/
  for (j = 0; j < _dim; j++){
    pars[j*npars+1] = pars[j*npars+2] = 0.0;
  }
  y_obs = regresResM;
  rp = rM;
  switch (_dim){
  case 1:
    for (obs = 0; obs < *nP; obs++){
      y_gamma = (*y_obs - _intcpt[0])/_scale[0] - _gamma[0];
      pars[1] += y_gamma*y_gamma;
      pars[2] += ((*rp) - _K[0])*y_gamma;
      y_obs++;
      rp++;
    }
    break;

  case 2:
    for (obs = 0; obs < *nP; obs++){
      y_gamma = (*y_obs - _intcpt[0])/_scale[0] - _gamma[0];
      pars[1] += y_gamma*y_gamma;
      pars[2] += (((*rp) % _length[0]) - _K[0])*y_gamma;
      y_obs++;

      y_gamma = (*y_obs - _intcpt[1])/_scale[1] - _gamma[1];
      pars[npars+1] += y_gamma*y_gamma;
      pars[npars+2] += (((*rp) / _length[0]) - _K[1])*y_gamma;
      y_obs++;
      rp++;
    }
    break;

  default:
    throw returnR("C++ Error: Gspline::full_sigma_pars not implemented for _dim > 2", 1);   
  }

  /*** Compute final values ***/
  for (j = 0; j < _dim; j++){
    pars[j*npars+1] *= 0.5;                        /* eta - eta_0, note that eta_0 = 0 for SDUnif prior */
    pars[j*npars+2] *= (_c4delta[j]/2);            /* xi/2                                              */

    switch (_prior_for_sigma[j]){
    case Fixed_:
      break;
    case Gamma:
      pars[j*npars]   = _prior_sigma[j*2] + (*nP)/2 - 1;       /* zeta - 1  */
      pars[j*npars+1] += _prior_sigma[j*2+1];                  /* eta       */
      pars[j*npars+1] = sqrt(pars[j*npars+1]);                 /* sqrt(eta) */
      break;
    case SDUnif:
      pars[j*npars]   = -0.5 + (*nP)/2 - 1;                    /* zeta - 1               */
      pars[j*npars+1] = sqrt(pars[j*npars+1]);                 /* sqrt(eta)              */
      pars[j*npars+3] = _prior_sigma[j*2+1];                   /* truncation limit 1/S^2 */
      break;
    default:
      throw returnR("C++ Error: Unknown prior appeared in Gspline::full_sigma_pars", 1);
    }
  }
  
  return;
}


// ***** full_sigma_logdens0: Compute log-density of a full conditional distribution    *****//
//                            of sigma^{-2}                                             *****//
//       
//   *** TRUNCATION: taken into account by      full_sigma_logdens0
//
// pars[4] ..... pars[0] = zeta - 1
//               pars[1] = sqrt(eta)
//               pars[2] = xi/2
//               pars[3] = truncation point 1/S^2
// ipars[1] .... ipars[0] = 0/1 is truncated?
//
void
full_sigma_logdens0(const double* x,  double* yu,  const double* pars,  const int* ipars)
{
  static double tmp, sqrt_x;

  if (ipars[0] && (*x) <= pars[3]){
    *yu = -FLT_MAX;
    return;
  }

  if (*x <= _epsilon){
    *yu = -FLT_MAX;
    return;
  }

  sqrt_x = sqrt(*x);
  tmp = pars[1]*sqrt_x - pars[2]/pars[1];
  *yu = pars[0]*log(*x) - tmp*tmp;
  return;
}


// ***** full_sigma_logdens3: Compute log-density (and possibly also some derivatives) of a full conditional distribution    *****//
//                            of sigma^{-2}                                                                                  *****//
//       
//   *** Curently: full_sigma_logdens3 is nowhere used
//
//   *** TRUNCATION: NOT taken into account by  full_sigma_logdens3
//
// pars[3] ..... pars[0] = zeta - 1
//               pars[1] = sqrt(eta)
//               pars[2] = xi/2
// ipars[1] .... not needed here (it is included here only to have the same prototype with related functions)
// what ........ 0 = compute l(x), l'(x), -l''(x)
//               1 = compute only l(x)
//               2 = compute only l'(x), -l''(x)
//               3 = compute only l(x), l'(x)
//
void
full_sigma_logdens3(const double* x,  double* yu,  double* ypu,  double* yppu,  const double* pars,  const int* ipars,  const int& what)
{  
  static double tmp, sqrt_x, invx, inv_sqrt_x, zeta_1_invx;

  sqrt_x = sqrt(*x);
  switch (what){
  case 0:
    tmp = pars[1]*sqrt_x - pars[2]/pars[1];
    *yu = pars[0]*log(*x) - tmp*tmp;

    invx = 1/(*x);
    zeta_1_invx = pars[0]*invx;
    inv_sqrt_x = 1/sqrt_x;
    *ypu = zeta_1_invx - pars[1]*pars[1] + pars[2]*inv_sqrt_x;

    *yppu = zeta_1_invx*invx + (pars[2]/2)*invx*inv_sqrt_x;
    return;

  case 1:
    tmp = pars[1]*sqrt_x - pars[2]/pars[1];
    *yu = pars[0]*log(*x) - tmp*tmp;
    return;

  case 2:
    invx = 1/(*x);
    zeta_1_invx = pars[0]*invx;
    inv_sqrt_x = 1/sqrt_x;
    *ypu = zeta_1_invx - pars[1]*pars[1] + pars[2]*inv_sqrt_x;

    *yppu = zeta_1_invx*invx + (pars[2]/2)*invx*inv_sqrt_x;
    return;

  case 3:
    tmp = pars[1]*sqrt_x - pars[2]/pars[1];
    *yu = pars[0]*log(*x) - tmp*tmp;

    invx = 1/(*x);
    zeta_1_invx = pars[0]*invx;
    inv_sqrt_x = 1/sqrt_x;
    *ypu = zeta_1_invx - pars[1]*pars[1] + pars[2]*inv_sqrt_x;
    return;

  default:
    throw returnR("C++ Error: incorrect 'what' in 'full_sigma_logdens3'", 1);
  }
}

