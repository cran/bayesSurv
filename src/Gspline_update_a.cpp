// ********************************************************************************* //
// ********* Gspline_update_a.cpp ************************************************** //
// ********************************************************************************* //
// Member functions of class Gspline and some related functions
//   for update of 'a' coefficients
// 
// 07/10/2004: 'Gspline::update_a'
//             'Gspline::update_alla'
//             'Gspline::find_start_abscis'
//             'Gspline::full_a_pars_uniCAR'
//             'Gspline::full_a_pars_eight_neighbors'
//             'Gspline::full_a_pars_twelve_neighbors' (not yet implemented)
//             'Gspline::find_eval_abscis'
//             'Gspline::check_abscis'
//             'Gspline::sample_a_by_slice'
//             'Gspline::sample_a_by_ARS'
//             'full_a_logdens0'
//             'full_a_logdens'
//             'full_a_logdens2'
//             'full_a_logdens3'
// 22/02/2005: 'Gspline::update_a' extended such that centering of a's is possible after update of each a
//             'Gspline::update_a' rewritten such that izero is ignored but it is checked that a's do not exceed
//               a given value (_a_ceil)
// 29/11/2006:  Function Gspline::update_alla changed to Gspline::update_alla_lambda which updates both a's and lambda(s)
//
#include "Gspline.h"


/***** update_a *****/
//
// Update one particular 'a' coefficient
//
// ija .............. index of 'a' that we want to update
//                    either one dimensional (if _dim == 1)
//                    or two dimensional (if _dim == 2)
// a_ipars[2] ....... a_ipars[0] = number of all observations
//                    a_ipars[1] = number of observations currently belonging to the (i,j)th component
// overrelax ........ 1/0 indicating whether overrelaxation is to be used (used only when the slice sampler is used,
//                    ignored otherwise)
//
void
Gspline::update_a(const int* ija,  const int* a_ipars,  const int* overrelax)
{
  static int ia, i;
  static double a_pars[4];
  static double newa;

  switch (_dim){
  case 1:
    ia = ija[0];
    break;
  case 2:
    ia = ija[1]*_length[0] + ija[0];
    break;
  default:
    throw returnR("C++ Error: Strange _dim in Gspline::update_a", 1);
  }

  full_a_pars(ija, a_pars + 0, a_pars + 1);          /* compute mean and inv. variance of [a[ia] | a[-ia], lambda] */
  a_pars[2] = _expa[ia];
  a_pars[3] = _sumexpa;

  /*** Find mode of the full conditional (if necessary)                                              ***/
  /***  store mode in _abscis[ia][1]                                                                 ***/
  /*** Compute either starting abscissae for ARS or initial guesses for interval defining the slice  ***/
  switch (_type_update_a){
  case Slice:
  case ARS_mode:
    find_eval_abscis(ia, a_pars, a_ipars);
    break;

  case ARS_quantile:    /** Starting abscissae are taken as quantiles of an upper hull from the previous iteration **/
                        /** Evaluate the function to sample from in starting abscissae                             **/
    for (i = 0; i < _nabscis; i++){
      full_a_logdens(_abscis[ia] + i, _hx + i, _hpx + i, a_pars, a_ipars);            
    }
    break;
  default: throw returnR("C++ Error: Unimplemented _type_update_a appeared in Gsplie::update_a", 1);
  }

  /** Check whether starting abscissae/initial guesses for the interval defining the slice lie on correct size of the mode   **/
  check_abscis(ia, a_pars, a_ipars);

  switch (_type_update_a){  
  case Slice:
    sample_a_by_slice(&newa, ia, a_pars, a_ipars, overrelax);
    break;

  case ARS_quantile:
  case ARS_mode:
    sample_a_by_ARS(&newa, ia, a_pars, a_ipars);
    break;    

  default:
    throw returnR("C++ Error: Unknown _type_update_a inGspline::update_a", 1);
  }

  change_a(&newa, ia);

    /***** Code used to keep _a[_izero] equal to zero *****/
    /***** * removed 22/02/2005                       *****/
    //    if (ia == _total_izero){     /* subtract newa from all a's to keep a[_total_izero] equal to 0 */
    //      _a[ia] = newa;
    //      for (i = 0; i < _total_length; i++) _a[i] -= newa;
    //      a2expa();                  /* recalculate exp(a), sum(exp(a)) etc. */
    //    }
    //    else{
    //      change_a(&newa, ia);       /* recalculate exp(a), sum(exp(a)), assign new 'a' to _a array */
    //    }

  return;
}  /** end of the function Gspline::update_a **/


/***** update_alla_lambda *****/
//
// Update all G-spline coefficients and also smoothing hyperparameter(s) lambda
//  Sequence of updating (if not in 1 block): try to minimize autocorrelation
//    * e.g. _dim = 1, _order = 3, update a[0], a[4], a[8], ...
//                                        a[1], a[5], a[9], ...
//                                        a[2], a[6], a[10], ...
//                                        a[3], a[7], a[11], ...
//           _dim = 2, _order = 3, update a[0,0], a[0,4], a[0,8], ...
//                                        a[4,0], a[4,4], a[4,8], ...
//                                        a[8,0], a[8,4], a[8,8],...
//                                        .....
//                                        a[1,0], a[1,4], a[1,8], ...
//                                        a[5,0], a[5,4], a[5,8], ...
//                                        a[9,0], a[9,4], a[9,8],...
//                                        .....
//                                        a[2,0], a[2,4], a[2,8], ...
//                                        a[6,0], a[6,4], a[6,8], ...
//                                        a[10,0], a[10,4], a[10,8],...
//                                        .....
//                                        a[3,0], a[3,4], a[3,8], ...
//                                        a[7,0], a[7,4], a[7,8], ...
//                                        a[11,0], a[11,4], a[11,8],...
//                                        .....
//                                        a[0,1], a[0,5], a[0,9], ...
//                                        a[4,1], a[4,5], a[4,9], ...
//                                        a[8,1], a[8,5], a[8,9],...
//                                        .....
//
//
// a_ipars[2] ......... a_ipars[0] = number of all observations
//                      a_ipars[1] = working place to store current mixtureNM
// mixtureNM .......... numbers of observations belonging to each component (long vector)  [_total_length]
// iter ............... number of iteration (to determine whether overrelaxation will be used)
//
void
Gspline::update_alla_lambda(const int* mixtureNM,  int* a_ipars,  const int* iter)
{
  static const int _ZERO_ = 0;
  static int accept;

  if (_order == 0) return;     // a's are fixed, and there is no lambda

  static int k0, k1;
  static int ija[2];
  static int ia;
  static int overrelax;
  ija[0] = ija[1] = 0;

  overrelax = 1*((*iter/_k_overrelax_a) != 0);

  switch (_dim){

  case 1:
    if (_type_update_a < Block){     /* Slice, ARS_quantile, ARS_mode */
      update_lambda();
      for (k0 = 0; k0 <= _order; k0++){
        for (ija[0] = k0; ija[0] < _length[0]; ija[0] += (_order + 1)){
          a_ipars[1] = mixtureNM[ija[0]];
          update_a(ija + 0, a_ipars, &overrelax);
        }
      }
      update_a_max_center_and_k_effect2006();      /** Original 'update_a_max_center_and_k_effect()' replaced on 01/12/2006 **/
      penalty();
    }
    else{                                /* _type_update_a == Block */
      //Rprintf("\n_a: ");
      //AK_BLAS_LAPACK::printArray(_a, _total_length);
      //Rprintf("_expa (sum(exp(a)=%g)): ", _sumexpa);
      //AK_BLAS_LAPACK::printArray(_expa, _total_length);
      //Rprintf("_w (_minw = %g): ", _minw);
      //AK_BLAS_LAPACK::printArray(_w, _total_length);
      //Rprintf("_Da: ");
      //AK_BLAS_LAPACK::printArray(_Da, _total_length);
      //Rprintf("_Qa (penalty=%g): ", *_penalty);
      //AK_BLAS_LAPACK::printArray(_Qa, _total_length);

      GMRF_Gspline::update(&accept, _a, _lambda, _expa, &_sumexpa, _w, &_minw, _Da, _Qa, _penalty, _workML, _worka, _workGMRF,
                           mixtureNM, _prior_for_lambda, _prior_lambda, _par_rscale, _Q, &_order, _diffOper, 
                           &GMRF_Gspline::_null_mass, &_constraint, _izero, &_total_length, a_ipars, &_ZERO_);
      if (accept) update_a_max_block();
    }
    return;  /** end of _dim == 1 **/

  case 2:
    update_lambda();
    for (k1 = 0; k1 <= _order; k1++){
      for (k0 = 0; k0 <= _order; k0++){
        for (ija[0] = k0; ija[0] < _length[0]; ija[0] += (_order + 1)){
          for (ija[1] = k1; ija[1] < _length[1]; ija[1] += (_order + 1)){
            ia = ija[1]*_length[0] + ija[0];
            a_ipars[1] = mixtureNM[ia];
            update_a(ija + 0, a_ipars, &overrelax);
          }
        }
      }
    }
    update_a_max_center_and_k_effect();  
    penalty();
    return;  /** end of _dim == 2 **/

  default:
    throw returnR("C++ Error: Strange _dim in Gspline::update_Gspline", 1);
  }  /** end of switch(_dim) **/

}  /** end of the function update_alla **/


/***** find_start_abscis *****/
// Find starting abscissae for a particular a
//  * use mean and mean +- 3*sd of the distribution [a[ia] | a[-ia], lambda]
//    as starting abscissae
//  * do not check whether they lie on both sides of the mode,
//    this will be done always before sampling is done
//
// ia ...... index of a for which we want to find the abscissae
//
void
Gspline::find_start_abscis(const int* ia)
{
  double a_pars[4];
  int ija[2];
  a_pars[2] = _expa[*ia];
  a_pars[3] = _sumexpa;
  switch (_dim){
  case 1:
    ija[0] = (*ia);
    break;
  case 2:
    ija[0] = (*ia) % _length[0];
    ija[1] = (*ia) / _length[0];
    break;
  }
  full_a_pars(ija + 0, a_pars + 0, a_pars + 1);
  double three_sd = 3/sqrt(a_pars[1]);
  if (_nabscis != 3) throw returnR("Dear Arnost, please update Gspline::find_start_abscis() function after changing _nabscis ;-)", 1);
  _abscis[*ia][0] = a_pars[0] - three_sd;
  _abscis[*ia][1] = a_pars[0];
  _abscis[*ia][2] = a_pars[0] + three_sd;
  
  return;
}


/***** full_a_pars_uniCAR *****/
//
// Compute mean and precision of the full conditional distribution of a[ija] in the case 
//  of univariate CAR in each dimension
//
// ija ...... index of 'a' for which we want to compute full conditional distribution
//            either one dimensional (if _dim == 1)
//            or two dimensional (if _dim == 2)
// mean ..... computed mean
// invvar ... computed inverse variance
//
void
Gspline::full_a_pars_uniCAR(const int* ija, double* mean, double* invvar) const
{
  int ia;
  int nr, nc;
  double mean1_2, mean2_1;
  double invvar1_2, invvar2_1;

  switch (_dim){

/*** dimension = 1 ***/
  case 1:
    ia = ija[0];
    if (ia < 0 || ia >= _total_length) throw returnR("C++ Error: Incorrect ija in Gspline:full_a_pars_uniCAR", 1);

    switch (_order){
    case 1:
      if (ia >= 1 && ia <= _length[0] - 2){
        *mean = (_a[ia-1] + _a[ia+1])/2;
        *invvar = 2*_lambda[0];
      }      
      else{     // ia = 0 or _length - 1
        if (ia == 0) *mean = _a[1];
        else         *mean = _a[_length[0]-2];
        *invvar = _lambda[0];
      }
      break;

    case 2:
      if (ia >= 2 && ia <= _length[0] - 3){
        *mean = (-_a[ia-2] + 4*_a[ia-1] + 4*_a[ia+1] - _a[ia+2])/6;
        *invvar = 6*_lambda[0];
      }
      else
        if (ia == 1 || ia == _length[0] - 2){
          if (ia == 1) *mean = (2*_a[0] + 4*_a[2] - _a[3])/5;
          else         *mean = (-_a[_length[0]-4] + 4*_a[_length[0]-3] + 2*_a[_length[0]-1])/5;
          *invvar = 5*_lambda[0];
        }
        else{     // ia = 0 or _length - 1
          if (ia == 0) *mean = 2*_a[1] - _a[2];
          else         *mean = -_a[_length[0]-3] + 2*_a[_length[0]-2];
          *invvar = _lambda[0];
        }
      break;

    case 3:
      if (ia >= 3 && ia <= _length[0] - 4){
        *mean = (_a[ia-3] - 6*_a[ia-2] + 15*_a[ia-1] + 15*_a[ia+1] - 6*_a[ia+2] + _a[ia+3])/20;
        *invvar = 20*_lambda[0];
      }
      else
        if (ia == 2 || ia == _length[0] - 3){
          if (ia == 2) *mean = (-3*_a[0] + 12*_a[1] + 15*_a[3] - 6*_a[4] + _a[5])/19;
          else         *mean = (_a[_length[0]-6] - 6*_a[_length[0]-5] + 15*_a[_length[0]-4] + 12*_a[_length[0]-2] - 3*_a[_length[0]-1])/19;
          *invvar = 19*_lambda[0];
        }
        else
          if (ia == 1 || ia == _length[0] - 2){
            if (ia == 1) *mean = (3*_a[0] + 12*_a[2] - 6*_a[3] + _a[4])/10; 
            else         *mean = (_a[_length[0]-5] - 6*_a[_length[0]-4] + 12*_a[_length[0]-3] + 3*_a[_length[0]-1])/10;
            *invvar = 10*_lambda[0];
          }
          else{     // ia = 0 or _length - 1
            if (ia == 0) *mean = 3*_a[1] - 3*_a[2] + _a[3];
            else         *mean = _a[_length[0]-4] - 3*_a[_length[0]-3] + 3*_a[_length[0]-2];
            *invvar = _lambda[0];
          }      
      break;

    default:
      throw returnR("C++ Error: Unimplemented _order appeared in Gspline::full_a_pars_uniCAR.", 1);
    }
    break;  /** end of _dim == 1 **/

/*** dimension = 2 ***/
  case 2:
    nr = _length[0];
    nc = _length[1];
    //    ia = ija[1]*nr + ija[0];
    if (ija[0] < 0 || ija[0] >= nr || ija[1] < 0 || ija[1] >= nc) throw returnR("C++ Error: Incorrect ija in Gspline:full_a_pars_uniCAR", 1);

    switch (_order){
    case 1:
      /** first conditional (with fixed column) **/
      if (ija[0] >= 1 && ija[0] <= nr - 2){
        mean1_2 = (_a[ija[1]*nr+ija[0]-1] + _a[ija[1]*nr+ija[0]+1])/2;
        invvar1_2 = 2*_lambda[0];
      }      
      else{     // ija[0] = 0 or nr - 1
        if (ija[0] == 0) mean1_2 = _a[ija[1]*nr+1];
        else             mean1_2 = _a[ija[1]*nr+nr-2];
        invvar1_2 = _lambda[0];
      }
      /** second conditional (with fixed row) **/
      if (ija[1] >= 1 && ija[1] <= nc - 2){
        mean2_1 = (_a[(ija[1]-1)*nr+ija[0]] + _a[(ija[1]+1)*nr+ija[0]])/2;
        invvar2_1 = 2*_lambda[1];
      }      
      else{     // ija[1] = 0 or nc - 1
        if (ija[1] == 0) mean2_1 = _a[1*nr+ija[0]];
        else             mean2_1 = _a[(nc-2)*nr+ija[0]];
        invvar2_1 = _lambda[1];
      }
      break;  /** end of _order == 2 **/

    case 2:
      /** first conditional (with fixed column) **/
      if (ija[0] >= 2 && ija[0] <= nr - 3){
        mean1_2 = (-_a[ija[1]*nr+ija[0]-2] + 4*_a[ija[1]*nr+ija[0]-1] + 4*_a[ija[1]*nr+ija[0]+1] - _a[ija[1]*nr+ija[0]+2])/6;
        invvar1_2 = 6*_lambda[0];
      }
      else
        if (ija[0] == 1 || ija[0] == nr - 2){
          if (ija[0] == 1) mean1_2 = (2*_a[ija[1]*nr+0] + 4*_a[ija[1]*nr+2] - _a[ija[1]*nr+3])/5;
          else             mean1_2 = (-_a[ija[1]*nr+nr-4] + 4*_a[ija[1]*nr+nr-3] + 2*_a[ija[1]*nr+nr-1])/5;
          invvar1_2 = 5*_lambda[0];
        }
        else{     // ija[0] = 0 or nr - 1
          if (ija[0] == 0) mean1_2 = 2*_a[ija[1]*nr+1] - _a[ija[1]*nr+2];
          else             mean1_2 = -_a[ija[1]*nr+nr-3] + 2*_a[ija[1]*nr+nr-2];
          invvar1_2 = _lambda[0];
        }
      /** second conditional (with fixed row) **/
      if (ija[1] >= 2 && ija[1] <= nc - 3){
        mean2_1 = (-_a[(ija[1]-2)*nr+ija[0]] + 4*_a[(ija[1]-1)*nr+ija[0]] + 4*_a[(ija[1]+1)*nr+ija[0]] - _a[(ija[1]+2)*nr+ija[0]])/6;
        invvar2_1 = 6*_lambda[1];
      }
      else
        if (ija[1] == 1 || ija[1] == nc - 2){
          if (ija[1] == 1) mean2_1 = (2*_a[0*nr+ija[0]] + 4*_a[2*nr+ija[0]] - _a[3*nr+ija[0]])/5;
          else             mean2_1 = (-_a[(nc-4)*nr+ija[0]] + 4*_a[(nc-3)*nr+ija[0]] + 2*_a[(nc-1)*nr+ija[0]])/5;
          invvar2_1 = 5*_lambda[1];
        }
        else{     // ija[1] = 0 or nc - 1
          if (ija[1] == 0) mean2_1 = 2*_a[1*nr+ija[0]] - _a[2*nr+ija[0]];
          else             mean2_1 = -_a[(nc-3)*nr+ija[0]] + 2*_a[(nc-2)*nr+ija[0]];
          invvar2_1 = _lambda[1];
        }
      break;  /** end of _order == 2 **/

    case 3:
      /** first conditional (with fixed column) **/
      if (ija[0] >= 3 && ija[0] <= nr - 4){
        mean1_2 = (_a[ija[1]*nr+ija[0]-3] - 6*_a[ija[1]*nr+ija[0]-2] + 15*_a[ija[1]*nr+ija[0]-1] + 
                   15*_a[ija[1]*nr+ija[0]+1] - 6*_a[ija[1]*nr+ija[0]+2] + _a[ija[1]*nr+ija[0]+3])/20;
        invvar1_2 = 20*_lambda[0];
      }
      else
        if (ija[0] == 2 || ija[0] == nr - 3){
          if (ija[0] == 2) mean1_2 = (-3*_a[ija[1]*nr+0] + 12*_a[ija[1]*nr+1] + 15*_a[ija[1]*nr+3] - 6*_a[ija[1]*nr+4] + _a[ija[1]*nr+5])/19;
          else             mean1_2 = (_a[ija[1]*nr+nr-6] - 6*_a[ija[1]*nr+nr-5] + 15*_a[ija[1]*nr+nr-4] + 12*_a[ija[1]*nr+nr-2] - 3*_a[ija[1]*nr+nr-1])/19;
          invvar1_2 = 19*_lambda[0];
        }
        else
          if (ija[0] == 1 || ija[0] == nr - 2){
            if (ija[0] == 1) mean1_2 = (3*_a[ija[1]*nr+0] + 12*_a[ija[1]*nr+2] - 6*_a[ija[1]*nr+3] + _a[ija[1]*nr+4])/10; 
            else             mean1_2 = (_a[ija[1]*nr+nr-5] - 6*_a[ija[1]*nr+nr-4] + 12*_a[ija[1]*nr+nr-3] + 3*_a[ija[1]*nr+nr-1])/10;
            invvar1_2 = 10*_lambda[0];
          }
          else{     // ija[0] = 0 or nr - 1
            if (ija[0] == 0) mean1_2 = 3*_a[ija[1]*nr+1] - 3*_a[ija[1]*nr+2] + _a[ija[1]*nr+3];
            else             mean1_2 = _a[ija[1]*nr+nr-4] - 3*_a[ija[1]*nr+nr-3] + 3*_a[ija[1]*nr+nr-2];
            invvar1_2 = _lambda[0];
          }      
      /** second conditional (with fixed row) **/
      if (ija[1] >= 3 && ija[1] <= nc - 4){
        mean2_1 = (_a[(ija[1]-3)*nr+ija[0]] - 6*_a[(ija[1]-2)*nr+ija[0]] + 15*_a[(ija[1]-1)*nr+ija[0]] + 
                   15*_a[(ija[1]+1)*nr+ija[0]] - 6*_a[(ija[1]+2)*nr+ija[0]] + _a[(ija[1]+3)*nr+ija[0]])/20;
        invvar2_1 = 20*_lambda[1];
      }
      else
        if (ija[1] == 2 || ija[1] == nc - 3){
          if (ija[1] == 2) mean2_1 = (-3*_a[0*nr+ija[0]] + 12*_a[1*nr+ija[0]] + 15*_a[3*nr+ija[0]] - 6*_a[4*nr+ija[0]] + _a[5*nr+ija[0]])/19;
          else             mean2_1 = (_a[(nc-6)*nr+ija[0]] - 6*_a[(nc-5)*nr+ija[0]] + 15*_a[(nc-4)*nr+ija[0]] + 
                                      12*_a[(nc-2)*nr+ija[0]] - 3*_a[(nc-1)*nr+ija[0]])/19;
          invvar2_1 = 19*_lambda[1];
        }
        else
          if (ija[1] == 1 || ija[1] == nc - 2){
            if (ija[1] == 1) mean2_1 = (3*_a[0*nr+ija[0]] + 12*_a[2*nr+ija[0]] - 6*_a[3*nr+ija[0]] + _a[4*nr+ija[0]])/10; 
            else             mean2_1 = (_a[(nc-5)*nr+ija[0]] - 6*_a[(nc-4)*nr+ija[0]] + 12*_a[(nc-3)*nr+ija[0]] + 3*_a[(nc-1)*nr+ija[0]])/10;
            invvar2_1 = 10*_lambda[1];
          }
          else{     // ija[1] = 0 or nc - 1
            if (ija[1] == 0) mean2_1 = 3*_a[1*nr+ija[0]] - 3*_a[2*nr+ija[0]] + _a[3*nr+ija[0]];
            else             mean2_1 = _a[(nc-4)*nr+ija[0]] - 3*_a[(nc-3)*nr+ija[0]] + 3*_a[(nc-2)*nr+ija[0]];
            invvar2_1 = _lambda[1];
          }      
      break;  /** end of _order == 3 **/

    default:
      throw returnR("C++ Error: Unimplemented _order appeared in Gspline::full_a_pars_uniCAR.", 1);
    }  /** end of switch (_order) **/

    *invvar = invvar1_2 + invvar2_1;        
    *mean = (invvar1_2*mean1_2 + invvar2_1*mean2_1) / (*invvar);

    break;  /** end of _dim == 2 **/
 
  default:
    throw returnR("C++ Error: Strange _dim in Gspline::full_a_pars_uniCAR", 1);
  }  /** end of switch (_dim) **/

  return;
}


/***** full_a_pars_eight_neighbors *****/
//
// Compute mean and precision of the full conditional distribution of a[ija] in the case 
//  of eight neighbors system and local quadratic smoothing
//
// * only _lambda[0] is used as a common lambda parameter
//
// ija ...... index of 'a' for which we want to compute full conditional distribution
//            either one dimensional (if _dim == 1)
//            or two dimensional (if _dim == 2)
// mean ..... computed mean
// invvar ... computed inverse variance
//
void
Gspline::full_a_pars_eight_neighbors(const int* ija, double* mean, double* invvar) const
{
  if (_dim != 2) throw returnR("C++ Error: Strange _dim appeares in Gspline::full_a_pars_eight_neighbors", 1);
  int nr = _length[0];
  int nc = _length[1];
  if (ija[0] < 0 || ija[0] >= nr || ija[1] < 0 || ija[1] >= nc) throw returnR("C++ Error: Incorrect ija in Gspline:full_a_pars_eight_neighbors", 1);

  if (ija[0] > 0 && ija[0] < nr-1 && ija[1] > 0 && ija[1] < nc-1){     /* not on the edges */
    *mean = (2*(_a[(ija[1]-1)*nr+ija[0]] + _a[(ija[1]+1)*nr+ija[0]] + _a[ija[1]*nr+ija[0]-1] + _a[ija[1]*nr+ija[0]+1]) - 
	       (_a[(ija[1]-1)*nr+ija[0]-1] + _a[(ija[1]-1)*nr+ija[0]+1] + _a[(ija[1]+1)*nr+ija[0]-1] + _a[(ija[1]+1)*nr+ija[0]+1]))/4;
    *invvar = 4*_lambda[0];
  }
  else  /* on one of the four edges */
    if (ija[0] == 0){  /* upper edge */
      if (ija[1] == 0){  /* (0, 0) upper left corner */
        *mean = _a[1*nr+0] + _a[0*nr+1] - _a[1*nr+1];
        *invvar = _lambda[0];
      }
      else
        if (ija[1] == nc - 1){  /* (0, nc-1) upper right corner */
          *mean = _a[(nc-2)*nr+0] + _a[(nc-1)*nr+1] - _a[(nc-2)*nr+1];
          *invvar = _lambda[0];
        }
        else{  /* (0; 1, ..., nc-2) upper edge but not the corner */
          *mean = (2*_a[ija[1]*nr+1] + (_a[(ija[1]-1)*nr+0] + _a[(ija[1]+1)*nr+0]) - (_a[(ija[1]-1)*nr+1] + _a[(ija[1]+1)*nr+1]))/2;
          *invvar = 2*_lambda[0];
        }
    }  /* end of the upper edge */
    else  
      if (ija[0] == nr - 1){  /* lower edge */
        if (ija[1] == 0){  /* (nr-1, 0) lower left corner */
          *mean = _a[1*nr+nr-1] + _a[0*nr+nr-2] - _a[1*nr+nr-2];
          *invvar = _lambda[0];
        }
        else
          if (ija[1] == nc - 1){  /* (nr-1, nc-1) lower right corner */
            *mean = _a[(nc-2)*nr+nr-1] + _a[(nc-1)*nr+nr-2] - _a[(nc-2)*nr+nr-2];
            *invvar = _lambda[0];
          }
          else{  /* (nr-1; 1, ..., nc-2) lower edge but not the corner */
            *mean = (2*_a[ija[1]*nr+nr-2] + (_a[(ija[1]-1)*nr+nr-1] + _a[(ija[1]+1)*nr+nr-1]) - (_a[(ija[1]-1)*nr+nr-2] + _a[(ija[1]+1)*nr+nr-2]))/2;
            *invvar = 2*_lambda[0];
          }
      }  /* end of the lower edge */
      else  
        if (ija[1] == 0){  /* left edge but not the corner */
          *mean = (2*_a[1*nr+ija[0]] + (_a[0*nr+ija[0]-1] + _a[0*nr+ija[0]+1]) - (_a[1*nr+ija[0]-1] + _a[1*nr+ija[0]+1]))/2;
          *invvar = 2*_lambda[0];
        }
        else{  /* right edge but not the corner */
          *mean = (2*_a[(nc-2)*nr+ija[0]] + (_a[(nc-1)*nr+ija[0]-1] + _a[(nc-1)*nr+ija[0]+1]) - (_a[(nc-2)*nr+ija[0]-1] + _a[(nc-2)*nr+ija[0]+1]))/2;
          *invvar = 2*_lambda[0];
        }       
  return;
}

void
Gspline::full_a_pars_twelve_neighbors(const int* ija, double* mean, double* invvar) const
{
  if (_dim != 2) throw returnR("C++ Error: Strange _dim appeares in Gspline::full_a_pars_eight_neighbors", 1);
  int nr = _length[0];
  int nc = _length[1];
  if (ija[0] < 0 || ija[0] >= nr || ija[1] < 0 || ija[1] >= nc) throw returnR("C++ Error: Incorrect ija in Gspline:full_a_pars_eight_neighbors", 1);

  if (ija[0] > 1 && ija[0] < nr-2 && ija[1] > 1 && ija[1] < nc-2){     /* not on the boundaries */
    *mean = (2*(_a[(ija[1]-1)*nr+ija[0]] + _a[(ija[1]+1)*nr+ija[0]] + _a[ija[1]*nr+ija[0]-1] + _a[ija[1]*nr+ija[0]+1]) +
	     (_a[(ija[1]-1)*nr+ija[0]-1] + _a[(ija[1]-1)*nr+ija[0]+1] + _a[(ija[1]+1)*nr+ija[0]-1] + _a[(ija[1]+1)*nr+ija[0]+1]) -
             (_a[(ija[1]-2)*nr+ija[0]] + _a[(ija[1]+2)*nr+ija[0]] + _a[ija[1]*nr+ija[0]-2] + _a[ija[1]*nr+ija[0]+2]))/8;    
    *invvar = 8*_lambda[0];
  }

  throw returnR("C++ Error: Not yet implemented (Gspline::full_a_pars_twelve_neighbors)", 1);
  return;
}


/****** find_eval_abscis: Find starting abscissae for the ARS or initial guess for the interval defining the slice ******/
//
// it modifies: _abscis[ia]:
//    on OUTPUT: _abscis[ia][0] = mode of full conditional - 2*sd
//               _abscis[ia][1] = mode of full conditional
//               _abscis[ia][2] = mode of full conditional + 2*sd
//
// ia ........... index of 'a' that will be updated using computed abscissae
// a_pars[4]..... double parameters for the full conditional distribution of a
// a_ipars[2] ... integer parameters for the full conditional distribution of a
//
void
Gspline::find_eval_abscis(const int& ia, const double* a_pars, const int* a_ipars)
{
  static double hppx;
  static int err_nr, iter_nr;

  _abscis[ia][1] = _a[ia];
  full_a_logdens3(_abscis[ia] + 1, _hx + 1, _hpx + 1, &hppx, a_pars, a_ipars, 0);
  try{
    newton_raphson(_abscis[ia] + 1, _hx + 1, _hpx + 1, &hppx, a_pars, a_ipars, full_a_logdens3, 
      	     &iter_nr, &_maxiter_nr, &_max_stephalf, &_toler_nr, &_epsilon, &err_nr);
  }
  catch(returnR){
    this->print();
    REprintf("Trap in find_eval_abscis: newton_raphson failed\n");
    throw;  
  }
  if (err_nr >= 3){
    REprintf("err_nr = %d\n", err_nr);
    REprintf("a = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             _a[ia], a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    this->print();
    throw returnR("Trap in Gspline::update_a: Unable to find a mode of the full conditional distribution", 1);
  }
  if (hppx <= _epsilon) hppx = _epsilon;
  hppx = 2/sqrt(hppx);                      /* approx. 2*std. deviation of the full conditional */
  _abscis[ia][0] = _abscis[ia][1] - hppx;
  _abscis[ia][2] = _abscis[ia][1] + hppx;
  full_a_logdens(_abscis[ia] + 0, _hx + 0, _hpx + 0, a_pars, a_ipars);
  full_a_logdens(_abscis[ia] + 2, _hx + 2, _hpx + 2, a_pars, a_ipars);

  return;
}


/***** check_abscis: Check whether starting abscissae/initial guesses for the interval defining the slice lie on correct size of the mode   **/
//
void 
Gspline::check_abscis(const int& ia, const double* a_pars, const int* a_ipars)
{
  static double step_left, step_right;
  static bool left_bad, right_bad;
  step_left = _abscis[ia][1] - _abscis[ia][0];
  step_right = _abscis[ia][_nabscis-1] - _abscis[ia][_nabscis-2];
  left_bad = right_bad = true;
  while (left_bad){
    if (_hpx[0] < _epsilon){
      _abscis[ia][0] -= step_left;
      full_a_logdens(_abscis[ia] + 0, _hx + 0, _hpx + 0, a_pars, a_ipars);
    }
    else
      left_bad = false;
  }
  while (right_bad){
    if (_hpx[_nabscis-1] > -_epsilon){
      _abscis[ia][_nabscis-1] += step_right;
      full_a_logdens(_abscis[ia] + _nabscis-1, _hx + _nabscis-1, _hpx + _nabscis-1, a_pars, a_ipars);
    }
    else
      right_bad = false;
  }
  return;
}


/***** sample_a_by_slice: Sample new a using a slice sampler *****/
//
// ASSUMPTION: _abscis[ia][1] = mode of the full conditional distribution
//             _abscis[ia][0] = mode - 2*sd
//             _abscis[ia][2] = mode + 2*sd
//
void
Gspline::sample_a_by_slice(double* newa, const int& ia, const double* a_pars, const int* a_ipars, const int* overrelax)
{
  static double horiz;
  static int i, iter_nr, err_nr;

  /*** Reshuffle _abscis[ia] such that the interval defining the slice will be at the beginnnig ***/
  /***  and evaluate h in current point - store it in _hx[2]                                    ***/  
  _abscis[ia][1] = _abscis[ia][2];
  _hx[1] = _hx[2];
  _hpx[1] = _hpx[2];
  full_a_logdens0(_a+ia, _hx+2, a_pars, a_ipars);    

  /*** Sample the horizontal level defining the slice (on log-scale), store it in horiz ***/
  horiz = _hx[2] - rexp(1);   

  /*** Find the interval defining the slice ***/
  for (i = 0; i < 1; i++){
    try{
      solver_newton_raphson(_abscis[ia] + i, _hx + i, _hpx + i, &horiz, a_pars, a_ipars, full_a_logdens3, 
                            &iter_nr, &_maxiter_solver_nr, &_toler_solver_nr, &_epsilon, &err_nr);
    }
    catch(returnR){
      this->print();
      REprintf("Trap in sample_a_by_slice: solver_newton_raphson failed\n");
      throw;
    }
    if (err_nr >= 3){
      REprintf("err_nr = %d\n", err_nr);
      REprintf("a = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
               _a[ia], a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
      this->print();
      throw returnR("Trap in Gspline::update_a: Unable to find an interval defining the slice", 1);
    }
  }

  /*** Sample the new point ***/
  if (*overrelax){
    Slice_sampler::ss_exact_overrelax(newa, _abscis[ia], _a+ia, &horiz, full_a_logdens0, a_pars, a_ipars);
  }
  else{
    Slice_sampler::ss_exact_sample(newa, _abscis[ia], _hx, _a+ia, &horiz, full_a_logdens0, a_pars, a_ipars);
  }

  return;
}


/***** sample_a_by_ARS: Sample new a using the ARS *****/
//
void
Gspline::sample_a_by_ARS(double* newa, const int& ia, const double* a_pars, const int* a_ipars)
{
  static int i, ifault, r_zero;
  static double hlb, hub;

  /***  Initialize arrays for ARS, ifault can be either 0 or 5, other values (1, 2, 3, 4) are not possible  ***/
  ifault = 1;
  ARS::initial_(&_ns, &_nabscis, &_emax, _abscis[ia], _hx, _hpx, &ZERO_INT, &hlb, &ZERO_INT, &hub, &ifault, _iwv, _rwv);
  if (ifault > 0){    /*  numerical non-log-concavity detected -> use slice sampler instead */
    sample_a_by_slice(newa, ia, a_pars, a_ipars, &ZERO_INT);
    return;
//    REprintf("ifault=%d;  a=%e, a_pars=%e, %e, %e, %e, a_ipars=%d, %d\n", 
//             ifault, _a[ia], a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
//    throw returnR("Trap in ARS: Numerical non-log-concavity (initial_) detected when updating 'a' coefficients, use slice sampler instead", 1);
  }

  /***  Sample  the new point ***/
  for (i = 0; i < _n_ARS_step; i++){
    r_zero = 0;
    ifault = 6;    
    while (ifault == 6){    /* while random number generator generates zero */
      ARS::sample_(_iwv, _rwv, full_a_logdens, a_pars, a_ipars, newa, &ifault);
      switch (ifault){
      case 5:    /*  numerical non-log-concavity detected -> use slice sampler instead */
        sample_a_by_slice(newa, ia, a_pars, a_ipars, &ZERO_INT);
        return;
//          REprintf("ifault=%d;  a=%e, a_pars=%e, %e, %e, %e, a_ipars=%d, %d\n", 
//		   ifault, _a[ia], a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
//          throw returnR("Trap in ARS: Numerical non-log-concavity (sample_) detected when updating 'a' coefficients, use slice sampler instead", 1);
      case 6:
        r_zero++;
        Rprintf("Warning: Random number generator generated zero during ARS.\n");
        if (r_zero >= 10) throw returnR("Trap in ARS: Too many zeros generated by the random number generator", 1);
        break;
      case 7:
        throw returnR("Trap in ARS: Numerical instability detected by sample_", 1);
      }
    }
  }

  /***  Compute starting abscissae for the next iteration of MCMC (if necessary)  ***/
  if (_type_update_a == ARS_quantile) ARS::quantile_(_iwv, _rwv, &_nabscis, _prob + 0, _abscis[ia], &ZERO_INT);
  return;
}



/***** full_a_logdens0  *****/
/***** full_a_logdens  *****/
/***** full_a_logdens2 *****/
/***** full_a_logdens3 *****/
// 
// full_a_logdens0: Compute only log-density
// full_a_logdens:  Compute log-density and its first derivative of the log of the full conditional of a
// full_a_logdens2: Compute also  minus second derivative of log-density
// full_a_logdens3: Compute selectively log-density or/and first and minus second derivative of log-density
//
// ai ......... new a, point where it should be evaluated
// yu ......... log-density
// ypu ........ derivative of the log-density
// yppu ....... minus second derivative of the log-density
// a_pars ..... a_pars[0] = mean of [a[ia] | a[-ia], lambda]
//              a_pars[1] = inverse variance of [a[ia] | a[-ia], lambda]
//              a_pars[2] = exp(old ai) 
//              a_pars[3] = sum(exp(old a)) (evaluated at current a vector)
// a_ipars .... a_ipars[0] = N     = number of all observations
//              a_ipars[1] = N[ia] = number of observations in a component ia
//
void
full_a_logdens0(const double* ai,  double* yu,  const double* a_pars,  const int* a_ipars)
{
  static double new_expai, new_sumexpa, a_min_A;
  if (*ai >= _log_inf){
    //REprintf("\na = %e leads to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //         *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    //throw returnR("Trap in full_a_logdens0.", 1);
    //Rprintf("\na = %e may lead to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //        *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    new_expai   = _exp_emax;
    new_sumexpa = _exp_emax;
  }
  else{
    new_expai   = exp(*ai);
    new_sumexpa = a_pars[3] - a_pars[2] + new_expai;   
  }

  a_min_A = (*ai) - a_pars[0];

  *yu = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;

  if (!R_finite(*yu)){
    REprintf("\na = %e, yu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in full_a_logdens0, NaN is not allowed.", 1);
  }
  return;
}

void
full_a_logdens(const double* ai,  double* yu,  double* ypu,  const double* a_pars,  const int* a_ipars)
{
  static double new_expai, new_sumexpa, a_min_A, new_wi;
  if (*ai >= _log_inf){
    //REprintf("\na = %e leads to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //         *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    //throw returnR("Trap in full_a_logdens.", 1);
    //Rprintf("\na = %e may lead to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //        *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    new_expai   = _exp_emax;
    new_sumexpa = _exp_emax;
  }
  else{
    new_expai   = exp(*ai);
    new_sumexpa = a_pars[3] - a_pars[2] + new_expai;   
  }

  a_min_A = (*ai) - a_pars[0];
  new_wi  = new_expai / new_sumexpa;

  *yu  = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
  *ypu = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;

  if (!R_finite(*yu)){
    REprintf("\na = %e, yu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in full_a_logdens, NaN is not allowed.", 1);
  }
  if (!R_finite(*ypu)){
    REprintf("\na = %e, yu = %e, ypu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, *ypu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in full_a_logdens, NaN is not allowed.", 1);
  }
  return;
}

void
full_a_logdens2(const double* ai,  double* yu,  double* ypu,  double* yppu,  const double* a_pars,  const int* a_ipars)
{
  static double new_expai, new_sumexpa, a_min_A, new_wi;

  if (*ai >= _log_inf){
    //REprintf("\na = %e leads to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //         *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    //throw returnR("Trap in full_a_logdens2.", 1);
    //Rprintf("\na = %e may lead to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //        *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    new_expai   = _exp_emax;
    new_sumexpa = _exp_emax;
  }
  else{
    new_expai   = exp(*ai);
    new_sumexpa = a_pars[3] - a_pars[2] + new_expai;   
  }

  a_min_A = (*ai) - a_pars[0];
  new_wi  = new_expai / new_sumexpa;

  *yu   = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
  *ypu  = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;
  *yppu = a_ipars[0]*new_wi*(1 - new_wi) + a_pars[1];

  if (!R_finite(*yu)){
    REprintf("\na = %e, yu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in full_a_logdens2, NaN is not allowed.", 1);
  }
  if (!R_finite(*ypu)){
    REprintf("\na = %e, yu = %e, ypu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, *ypu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in full_a_logdens2, NaN is not allowed.", 1);
  }
  if (!R_finite(*yppu)){
    REprintf("\na = %e, yu = %e, ypu = %e, yppu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, *ypu, *yppu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in full_a_logdens2, NaN is not allowed.", 1);
  }
  return;
}

//
// what .... 0 = compute l(x), l'(x), -l''(x)
//           1 = compute only l(x)
//           2 = compute only l'(x), -l''(x)
//           3 = compute only l(x), l'(x)
//
void
full_a_logdens3(const double* ai,  double* yu,  double* ypu,  double* yppu,  const double* a_pars,  const int* a_ipars,  const int& what)
{
  static double new_expai, new_sumexpa, a_min_A, new_wi;

  if (*ai >= _log_inf){
    //Rprintf("\na = %e may lead to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //        *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    new_expai   = _exp_emax;
    new_sumexpa = _exp_emax;
  }
  else{
    new_expai   = exp(*ai);
    new_sumexpa = a_pars[3] - a_pars[2] + new_expai;
  }
  a_min_A = (*ai) - a_pars[0];
  new_wi = new_expai / new_sumexpa;

  switch (what){
  case 0:
    *yu   = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
    *ypu  = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;
    *yppu = a_ipars[0]*new_wi*(1 - new_wi) + a_pars[1];
    return;
  case 1:
    *yu = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
    return;
  case 2:
    *ypu  = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;
    *yppu = a_ipars[0]*new_wi*(1 - new_wi) + a_pars[1];
    return;
  case 3:
    *yu  = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
    *ypu = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;
    return;
  default:
    throw returnR("C++ Error: incorrect 'what' in 'full_a_logdens3'", 1);
  }
}
