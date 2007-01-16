// Header for Gspline class
//   * a class to store G-spline and information to update it
//
//
// 07/10/2004: start working on it
// 20/02/2005: 'mean_univariate' added
//
#ifndef _GSPLINE_H_
#define _GSPLINE_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <iomanip>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "random.h"
#include "newton_raphson.h"
//#include "slice_sampler.h"
#include "Slice_sampler2.h"
//#include "ars.h"
#include "ARS2.h"
#include "GMRF_Gspline.h"
#include "GMRF_Gspline_Util.h"

enum neighborSystem {uniCAR, eight_neighbors, twelve_neighbors};            // types of neighboring systems for 'a' coefficients
enum priorForSigmaTypes {Fixed_, Gamma, SDUnif};                            // types of priors for sigma
enum priorForGammaTypes {_Fixed_, Normal};                                  // types of priors for gamma
enum a_sampling_scheme {Slice, ARS_quantile, ARS_mode, Block};              // sampling schemes for a

const int _max_dim = 2;                      // maximal dimension of the G-spline
const bool _center_a = false;                // if true, 'a' coefficients are maintained such that their average is always zero
                                             //  in that case, a[_izero] is not necessary equal to zero
const double _a_ceil = 10.0;                 // upper limit for 'a' coefficients (not applied if _center_a == true)
const double _amax_floor = 0.0;              // lower limit for maximal 'a' coefficient (not applied if _center_a == true)

/** constants for adaptive rejection sampling **/
const int _nabscis = 3;                      // number of starting abscissae for each ARS
const int _ns = 10;                          // maximal number of abscissae for each ARS
const int _n_ARS_step = 5;                   // number of sampled points within each ARS (the last one is taken as the sampled value)

/** constants for Newton-Raphson within ARS **/
const int _maxiter_nr = 10;                  // maximum number of NR iterations
const int _max_stephalf = 10;                // maximum number of step-halving steps
const double _toler_nr = 1e-3;               // tolerance to detect convergence of NR

/** constants for Newton-Raphson solver within the slice sampler **/
const int _maxiter_solver_nr = 10;           // maximum number of NR iterations
const double _toler_solver_nr = 1e-3;        // tolerance to detect convergence of NR

/** some other constants **/
const double _emax      = 64.0;                   // exp(-_emax) = 0.0
const double _epsilon   = exp(-_emax);
const double _exp_emax  = exp(_emax); 
const double _e_min_inf = 1e-300;                 // log(_e_min_inf) = -inf
const double _log_zero  = log(_e_min_inf);
const double _log_inf   = _emax;
const double _prob[3]   = {0.15, 0.50, 0.85};  // probs for quantiles of the upper hull to compute starting abscissae for the next iteration

const double _null_mass = 1e-6;              // area with this probability will be assigned a zero probability in some applications
                                             //   (point from such area would appear only every 1e6th iteration and I guess we will run 
                                             //   only rarely MCMC with more than 1 000 000 iterations)
const double _zero_invvariance = 1e-5;       // lower limit for inverse variance parameters (inverse scale or inverse sigma,
                                             //   used by Gspline_update_Scale and Gspline_update_sigma)

class Gspline
{
 private:
  int  _dim;                  // dimension of the G-spline  (only 1 or 2 will be allowed)                         1 
  int  _neighbor_system;      // neighboring system for 'a' coefficients                                          1
  bool _equal_lambda;         // are lambda's for both dimensions (if _dim == 2) equal?                           1

  int  _total_length;         // total length of the G-spline                                                     1
  int* _length;               // lengths of the G-spline in each dimension (it should be odd)                     _dim
  int* _K;                    // (_length[j] - 1)/2 for each dimension                                            _dim

  int  _total_izero;          // index of 'a' coefficient which is zero in a long vector                          1
  int* _izero;                // indeces of 'a' coefficients which are constantly zero (for each dimension)       _dim
                              // or which are used as a reference    

  int  _order;                // order of autoregression in the prior for 'a'                                     1
                              //  if equal to 0, a coefficients are assumed to be fixed
  double _log_null_w;         // 'a' coefficient which satisfies (a - a_max < _log_null_w) will be assumed to lead
                              //  to zero weight                                                                  1
                              // it is defined such that (_total_length*exp(_log_null_w))/(1+_total_length*exp(_log_null_w)) = _null_mass

  /** 'a' coefficients and related quantities **/
  double* _lambda;           // lambda precision parameters for each dimension                                   _dim
  double* _a;                // a coefficients (in a long array in column major order)                           _total_length
  double  _a_max;            // max(a)                                                                           1
  double* _expa;             // proportions of weights = exp(a) (in a long array in column major order)          _total_length
  double  _sumexpa;          // sum(exp(a))                                                                      1
  double** _sumexpa_margin;  // sum of exp(a) over all indeces up to one index                                   _dim, _length[j]
                             //  * effectively present only if _dim >= 2         
  double* _penalty;          // -0.5*sum (Delta a)^2 in the case of univariate CAR in each dimension              _dim
                             //  * if _equal_lambda then _penalty[0] gives the sum of row- and column-penalties
                             //  * if _equal_lambda then _penalty[1], ... are sometimes used as working space
  int     _k_effect;         // number of weights that correspond to non-zero mass                               1
  int*    _ind_w_effect;     // indeces of weights that correspond to non-zero mass                              _total_length

  /** variables for ARS/slice sampler update of a's **/
  double** _abscis;     // array of arrays of starting abscissae for ARS for each a                                    _total_length, _nabscis 
                        //   it is used to store intervals defining the slice when the slice sampler is used
  int* _iwv;            // integer working place for ARS                                                               7 + _ns
  double* _rwv;         // double working place for ARS                                                                9 + 6*(1+_ns)
  double* _hx;          // working array for ARS/slice sampler/update_sigma                                            _nabscis
  double* _hpx;         // working array for ARS/slice sampler/update_sigma                                            _nabscis
  int _type_update_a;   // type of the a update (slice sampler, two types of ARS or in 1 block using MH)               1
  int _k_overrelax_a;   // how often (every kth) will be performed a normal update of a instead of overrelaxed one     1
                        // (applicable only for a slice sampler, ignored otherwise)

  /** variables for slice sampler update of sigma's/scale's **/
  int _k_overrelax_sigma[_max_dim];     // how often (every kth) will be performed a normal update instead of overrelaxed one
  int _k_overrelax_scale[_max_dim];     // how often (every kth) will be performed a normal update instead of overrelaxed one

  /** parameters specifying priors **/
  int _prior_for_lambda[_max_dim];      // type of prior for lambda(s)                                                            _max_dim
  double _prior_lambda[2*_max_dim];     // prior shape and rate for update of lambda(s) (if uniform prior on (0, S)               2*_max_dim
                                        //   for lambda[j]^{-1/2} is used then we store 1/S^2 in _prior_lambda[2*j+1]

  int _prior_for_gamma[_max_dim];       // type of prior for gamma(s)                                                             _max_dim
  double _prior_gamma[2*_max_dim];      // prior mean and inverted variance for gamma(s) (overall mean in each dimension)         2*_max_dim

  int _prior_for_sigma[_max_dim];       // type of prior for sigma(s)                                                             _max_dim
  double _prior_sigma[2*_max_dim];      // prior shape and rate for update of sigma2 (if uniform prior on (0, S)
                                        //   for sigma[j] is used then we store 1/S^2 in _prior_sigma[2*j+1]

  int _prior_for_intcpt[_max_dim];      // type of prior for intercept(s) (2nd specification)                                     _max_dim
  double _prior_intcpt[2*_max_dim];     // prior mean and inverted variance for intercepts                                        2*_max_dim

  int _prior_for_scale[_max_dim];       // type of prior for scale(s) (2nd specification)                                         _max_dim
  double _prior_scale[2*_max_dim];      // prior shape and rate for update of scale2 (if uniform prior on (0, S)
                                        //   for scale[j] is used then we store 1/S^2 in _prior_scale[2*j+1]


  /** mixture parameters **/
  double* _gamma;       // middle knot for each dimension                                                  _dim
  double* _invsigma2;   // basis sigma^{-2} for each dimension                                             _dim
  double* _sigma;       // basis sigma for each dimension                                                  _dim
  double* _c4delta;     // 'c' coefficient to compute the distance between two knots for each dimension    _dim
  double* _delta;       // distance between two knots for each dimension                                   _dim
  double* _intcpt;      // intercept for each dimension (2nd specification)                                _dim
  double* _invscale2;   // scale^{-2} for each dimension (2nd specification)                               _dim
  double* _scale;       // scale for each dimension (2nd specification)                                    _dim

  /** quantities needed for joint update of a coefficients **/
  /** !!! Currently only implemented for _dim = 1 !!!      **/
  /** !!! Set to NULL for _dim > 1 !!!                     **/
  double *_w;                // weights (in a long array in column major order)                           _total_length
  double _minw;              // minimal weight                                                            1
  double *_Q;                // Q = t(D)*D                                                                LT(_total_length)
  double *_Da;               // D*a in first _total_length - _order places                                _total_length
  double *_Qa;               // Q*a = t(D)*D*a                                                            _total_length
  int *_diffOper;            // vector defining the difference operator                                   _order + 1 (if _dim = 1)
  int _constraint;           // type of the constraint (see GMRF_Gspline_Util.h for possible values)

  double _par_rscale[6];     // parameters for the lambda proposal according to Rue and Held 
                             // actually not used in this version

  int _LTna;                  // = LT(_total_length)
  int _LTna_1;                // = LT(_total_length - 1)
  int _nworkML;
  int _nworka;
  int _nworkGMRF;
  double *_workML;           // working arrays, see GMRF_Gspline.cpp for how long they should be
  double *_worka;
  double *_workGMRF;

 public:

// ======================================================================================
// ***** Gspline.cpp: basic functions for the class Gspline
// ======================================================================================

/***** Constructors and destructors *****/
  Gspline();
  Gspline(const int* parmI,   const double* parmD);
  Gspline(const Gspline& gg);
  Gspline& operator=(const Gspline& gg);
  ~Gspline();

/*****  Access functions to various components of the G-spline *****/
  inline int
  dim() const { return _dim;}

  inline int
  K(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::K(j).", 1);
    return _K[j];
  }

  inline bool
  equal_lambda() const { return _equal_lambda;}


  inline int
  total_length() const { return _total_length;}

  inline int
  length(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::length(j).", 1);
    return _length[j];
  }

  inline int
  total_izero() const { return _total_izero;}

  inline int
  izero(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::izero(j).", 1);
    return _izero[j];
  }

  inline int
  order() const { return _order;}

  inline double
  log_null_w() const { return _log_null_w;}

  inline double
  lambda(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::lambda(j).", 1);
    return _lambda[j];
  }

  inline double*
  lambdaP() const { return _lambda;}

  inline double
  gamma(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::gamma(j).", 1);
    return _gamma[j];
  }

  inline double*
  gammaP() const { return _gamma;}

  inline double
  intcpt(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::intcpt(j).", 1);
    return _intcpt[j];
  }

  inline double*
  intcptP() const { return _intcpt;}

  inline double
  invsigma2(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::invsigma2(j).", 1);
    return _invsigma2[j];
  }

  inline double
  sigma(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::sigma(j).", 1);
    return _sigma[j];
  }

  inline double*
  sigmaP() const { return _sigma;}

  inline double
  invscale2(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::invscale2(j).", 1);
    return _invscale2[j];
  }

  inline double
  scale(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::scale(j).", 1);
    return _scale[j];
  }

  inline double*
  scaleP() const { return _scale;}

  inline double
  c4delta(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::c4delta(j).", 1);
    return _c4delta[j];
  }

  inline double
  delta(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::delta(j).", 1);
    return _delta[j];
  }

  inline double*
  deltaP() const { return _delta;}

  inline double
  a_max() const { return _a_max;}

  inline double*
  a_maxP() { return &_a_max;}

  inline double
  sumexpa() const { return _sumexpa;}

  inline double*
  sumexpaP() { return &_sumexpa;}

  inline int
  k_effect() const { return _k_effect;}

  inline double*
  aP() const { return _a;}

  inline double*
  expaP() const { return _expa;}

  inline double
  penalty(const int& j) const
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::penalty(j).", 1);
    return _penalty[j];
  }

  inline double*
  penaltyP() const { return _penalty;}

  inline double
  suma() const {
    double *aP = _a;
    double _suma_ = 0.0;
    for (int i = 0; i < _total_length; i++){
      _suma_ += *aP;
      aP++;
    }
    return _suma_;
  }


/***** Indexing 'operators' *****/
  inline double
  a(const int& ia) const
  {
    if (ia < 0 || ia >= _total_length) throw returnR("C++ Error: Incorrect ia in Gspline:a(ia)", 1);
    return _a[ia];
  }

  inline double
  expa(const int& ia) const
  {
    if (ia < 0 || ia >= _total_length) throw returnR("C++ Error: Incorrect ia in Gspline:expa(ia)", 1);
    return _expa[ia];
  }  
  
  inline double
  w(const int& ia) const
  {
    if (ia < 0 || ia >= _total_length) throw returnR("C++ Error: Incorrect ia in Gspline:w(ia)", 1);
    return _expa[ia]/_sumexpa;
  }

  inline double
  mu_component(const int& j, const int& ia) const     /* return jth component of the mu vector corresponding to a[ia] */
  {
    if (ia < 0 || ia >= _total_length) throw returnR("C++ Error: Incorrect ia in Gspline:mu_component(j, ia)", 1);
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::mu_component(j, ia).", 1);
    switch (_dim){
    case 1:
      return _gamma[0] + (ia - _K[0])*_delta[0];
    case 2:
      if (j == 0) return _gamma[0] + ((ia % _length[0]) -_K[0])*_delta[0];
      else        return _gamma[1] + ((ia / _length[0]) -_K[1])*_delta[1];
    default:
      throw returnR("C++ Error: Incorrect use of Gspline::mu_component function", 1);
    }
  }

  inline int ind_w_effect
  (const int& i) const
  {
    if (i < 0 || i >= _k_effect) throw returnR("C++ Error: Incorrect i in Gspline:ind_w_effect(i)", 1);
    return _ind_w_effect[i];
  }  

  inline double
  a(const int& k1, const int& k2) const   /* parameters kj, j=1,2 are assumed to take values 0, 1, ..., _length[j]-1 */
  {
    if (_dim != 2) throw returnR("C++ Error: Incorrect use of Gspline:a(k1, k2)", 1);
    if (k1 < 0 || k1 >= _length[0]) throw returnR("C++ Error: Incorrect k1 in Gspline:a(k1, k2)", 1);
    if (k2 < 0 || k2 >= _length[1]) throw returnR("C++ Error: Incorrect k2 in Gspline:a(k1, k2)", 1);
    return _a[k2 * _length[0] + k1];
  }

  inline double
  expa(const int& k1, const int& k2) const   /* parameters kj, j=1,2 are assumed to take values 0, 1, ..., _length[j]-1 */
  {
    if (_dim != 2) throw returnR("C++ Error: Incorrect use of Gspline:expa(k1, k2)", 1);
    if (k1 < 0 || k1 >= _length[0]) throw returnR("C++ Error: Incorrect k1 in Gspline:expa(k1, k2)", 1);
    if (k2 < 0 || k2 >= _length[1]) throw returnR("C++ Error: Incorrect k2 in Gspline:expa(k1, k2)", 1);
    return _expa[k2 * _length[0] + k1];
  }

  inline double
  w(const int& k1, const int& k2) const   /* parameters kj, j=1,2 are assumed to take values 0, 1, ..., _length[j]-1 */
  {
    if (_dim != 2) throw returnR("C++ Error: Incorrect use of Gspline:w(k1, k2)", 1);
    if (k1 < 0 || k2 >= _length[0]) throw returnR("C++ Error: Incorrect k1 in Gspline:w(k1, k2)", 1);
    if (k1 < 0 || k2 >= _length[1]) throw returnR("C++ Error: Incorrect k2 in Gspline:w(k1, k2)", 1);
    return _expa[k2 * _length[0] + k1]/_sumexpa;
  }

  inline double
  sumexpa_margin(const int& j, const int& k) const   /* parameter k is assumed to take values 0, 1, ..., _length[j]-1 */
  {
    if (_dim < 2) throw returnR("C++ Error: Incorrect use of Gspline::sumexpa_margin(j, k)", 1);
    if (k < 0 || k >= _length[j]) throw returnR("C++ Error: Incorrect k in Gspline:sumexpa_margin(j, k)", 1);
    return _sumexpa_margin[j][k];
  }

  inline double
  abscis(const int& ia, const int& k) const
  {
    if (ia < 0 || ia >= _total_length) throw returnR("C++ Error: Incorrect ia in Gspline:abscis(ia, k)", 1);
    if (k < 0 || k >= _nabscis) throw returnR("C++ Error: Incorrect k in Gspline:abscis(ia, k)", 1);
    return _abscis[ia][k];
  }

  inline double
  mu(const int& j, const int& k) const      /* parameter k is assumed to take values 0, 1, ..., _length[j]-1 */
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::mu(j, k).", 1);
    if (k >= _length[j] || k < 0) throw returnR("C++ Error: Incorrect k in Gspline::mu(j, k).", 1);
    return _gamma[j] + (k - _K[j])*_delta[j];
  }

  /** Compute all knots in dimension j and return them in an array muj **/
  inline void
  muP(const int& j, double* muj) const 
  {
    if (j >= _dim || j < 0) throw returnR("C++ Error: Incorrect j in Gspline::muP(j, muj).", 1);
    double* knot = muj;
    (*knot) = _gamma[j] - _K[j]*_delta[j]; 
    for (int i = 1; i < _length[j]; i++){
      knot++;
      (*knot) = *(knot-1) + _delta[j];
    }
    return;
  }

/***** Other functions *****/
  void
  print() const;

  void
  Gspline2initArray(int* parmI,  double* parmD) const;

  void
  a2expa();

  void
  a2expa_total_length();

  void
  change_a(const double* newa, const int& ia);

  void
  update_k_effect();

  void
  update_a_max();

  void
  update_a_max_block();

  void
  update_a_max_center_and_k_effect();

  void
  update_a_max_center_and_k_effect2006();

  void
  moments(double* EY, double* varY) const;

  void
  mean_univariate(double* EY) const;

  void
  adjust_intcpt(double* adj);

// ===============================================================================
// ***** Gspline_update_a.cpp: Functions for update of 'a' coefficients    ***** //
// ===============================================================================
  void
  find_start_abscis(const int* ia);

  inline void
  full_a_pars(const int* ija, double* mean, double* invvar) const
  {
    switch (_neighbor_system){
    case uniCAR:           this->full_a_pars_uniCAR(ija, mean, invvar);           break;
    case eight_neighbors:  this->full_a_pars_eight_neighbors(ija, mean, invvar);  break;
    case twelve_neighbors: this->full_a_pars_twelve_neighbors(ija, mean, invvar); break;
    default:               throw returnR("C++ Error: Strange _neighbor_system in Gspline::full_a_pars", 1);
    }
  }
  
  void
  full_a_pars_uniCAR(const int* ija, double* mean, double* invvar) const;

  void
  full_a_pars_eight_neighbors(const int* ija, double* mean, double* invvar) const;

  void
  full_a_pars_twelve_neighbors(const int* ija, double* mean, double* invvar) const;

  void
  find_eval_abscis(const int& ia, const double* a_pars, const int* a_ipars);

  void 
  check_abscis(const int& ia, const double* a_pars, const int* a_ipars);

  void
  sample_a_by_slice(double* newa, const int& ia, const double* a_pars, const int* a_ipars, const int* overrelax);

  void
  sample_a_by_ARS(double* newa, const int& ia, const double* a_pars, const int* a_ipars);

  void
  update_a(const int* ija,  const int* a_ipars,  const int* overrelax);

  void
  update_alla_lambda(const int* mixtureNM,  int* a_ipars,  const int* iter);

// =============================================================================================
// ***** Gspline_update_lambda.cpp: Functions for update of precision parameters lambda ***** //
// =============================================================================================
  inline void
  penalty() const
  {
    switch (_neighbor_system){
    case uniCAR:           this->penalty_uniCAR();           break;
    case eight_neighbors:  this->penalty_eight_neighbors();  break;
    case twelve_neighbors: this->penalty_twelve_neighbors(); break;
    default:               throw returnR("C++ Error: Strange _neighbor_system in Gspline::penalty", 1);
    }
  }

  void
  penalty_uniCAR() const;

  void
  penalty_eight_neighbors() const;

  void
  penalty_twelve_neighbors() const;

  void
  update_lambda();

// ====================================================================================================
// ***** Gspline_update_sigma.cpp: Functions for update of basis std. deviations sigma         ***** //
// ====================================================================================================
  void
  full_sigma_pars(double* pars,  const double* regresResM, const int* rM, const int* nP) const;

  void
  update_sigma(const double* regresResM,  const int* rM,  const int* nP, const int* iter);

// ===========================================================================================================
// ***** Gspline_update_Scale.cpp: Functions for update of the G-spline scale (2nd specification)     ***** //
// ===========================================================================================================
  void
  update_Scale(const double* regresResM,  const int* rM,  const int* nP,  const int* iter);

  void
  full_Scale_pars(double* pars,  const double* regresResM,  const int* rM,  const int* nP) const;


// ============================================================================================
// ***** Gspline_update_gamma.cpp: Functions for update of overall means gamma         ***** //
// ============================================================================================
  void
  update_gamma(const double* regresResM, const int* rM, const int* nP);


// ======================================================================================================
// ***** Gspline_update_Intcpt.cpp: Functions for update an intercept (2nd specification)        ***** //
// =====================================================================================================
  void
  update_Intcpt(const double* regresResM,  const int* rM,  const int* nP);

};    // end of the class Gspline

/***** Additional functions related to class Gspline *****/
void
full_a_logdens0(const double* ai,  double* yu,  const double* a_pars,  const int* a_ipars);

void
full_a_logdens(const double* ai,  double* yu,  double* ypu,  const double* a_pars,  const int* a_ipars);

void
full_a_logdens2(const double* ai,  double* yu,  double* ypu,  double* yppu,  const double* a_pars,  const int* a_ipars);

void
full_a_logdens3(const double* ai,  double* yu,  double* ypu,  double* yppu,  const double* a_pars,  const int* a_ipars,  const int& what);

void
diff(const int& order, const int& na, double* Da);

void
full_sigma_logdens0(const double* x,  double* yu,  const double* pars,  const int* ipars);

void
full_sigma_logdens3(const double* x,  double* yu,  double* ypu,  double* yppu,  const double* pars,  const int* ipars,  const int& what);

#endif
