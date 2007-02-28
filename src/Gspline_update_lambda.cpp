// ************************************************************************************** //
// ********* Gspline_update_lambda.cpp ************************************************** //
// ************************************************************************************** //
// Member function of class Gspline to update precision parameters 'lambda'
//
// 18/10/2004: 'Gspline::update_lambda'
//             'Gspline::penalty_uniCAR'
//             'Gspline::penalty_eight_neighbors'
//             'Gspline::penalty_twelve_neighbors' (not yet implemented)
//             'diff'
//             'diff_in_R'
//
#include "Gspline.h"


/***** update_lambda *****/
void
Gspline::update_lambda()
{
  if (_order == 0) return;     // a's are fixed -> there is no random lambda

  int k;

  double shape, rate, scale;
  if (_dim == 1){                    /** Make also correction for the order of the penalty **/
    switch (_prior_for_lambda[0]){
    case GMRF_Gspline_Util::_Fixed_:
      return;

    case GMRF_Gspline_Util::_Gamma_:
      shape = _prior_lambda[0] + 0.5*(_total_length - _order + 1);
      rate = _prior_lambda[1] - _penalty[0];
      if (rate <= 0) throw returnR("Gspline::update_lambda: Trap in update of lambda (non-positive rate parameter)", 1);
      scale = 1/rate;
      _lambda[0] = rgamma(shape, scale);
      return;

    case GMRF_Gspline_Util::_SDUnif_:
      shape = 0.5*(_total_length - _order);
      rate = -_penalty[0];
      rltruncGamma(_lambda, &shape, &rate, _prior_lambda + 1, &ONE_INT, &ZERO_INT);
      return;

    default:     
      throw returnR("C++ Error: Unimplemented prior for lambda appeared in Gspline::update_lambda", 1);
    }
  }
  else{                             /** Correction for the order of the penalty is not done!!! **/
    for (k = 0; k <= (_dim-1)*(!_equal_lambda); k++){
      switch (_prior_for_lambda[k]){
      case GMRF_Gspline_Util::_Fixed_:
        break;

      case GMRF_Gspline_Util::_Gamma_:
        shape = _prior_lambda[2*k] + 0.5*(_total_length - 1);
        rate = _prior_lambda[2*k+1] - _penalty[k];
        if (rate <= 0) throw returnR("Gspline::update_lambda: Trap in update of lambda (non-positive rate parameter)", 1);
        scale = 1/rate;
        _lambda[k] = rgamma(shape, scale);
        break;

      case GMRF_Gspline_Util::_SDUnif_:
        shape = 0.5*(_total_length - 2);
        rate = -_penalty[k];
        rltruncGamma(_lambda + k, &shape, &rate, _prior_lambda + 2*k + 1, &ONE_INT, &ZERO_INT);           
        break;

      default:     
        throw returnR("C++ Error: Unimplemented prior for lambda appeared in Gspline::update_lambda", 1);
      }
    }
  }
}


/***** penalty_uniCAR *****/
//
// Compute -0.5*sum (Delta a)^2 in the case of univariate CAR in each dimension
//
// penalty ....... either one-component or two-component array
//                 * if _equal_lambda then penalty[0] gives the sum of row- and column-penalties
//
void
Gspline::penalty_uniCAR() const
{
  int ia, col, row;
  double* Da;

  switch (_dim){

    /*** dimension = 1 ***/
  case 1:
    Da = new double[_length[0]];
    if (Da == NULL) throw returnR("C++ Error: Could not allocate working memory", 1);
    for (ia = 0; ia < _length[0]; ia++) Da[ia] = _a[ia];
    diff(_order, _length[0], Da);
    _penalty[0] = 0.0;
    for (ia = 0; ia < _length[0] - _order; ia++) _penalty[0] += (Da[ia] * Da[ia]);
    _penalty[0] *= (-0.5);
    break;

  /*** dimension = 2 ***/
  case 2:
    Da = new double[(_length[0] > _length[1] ? _length[0] : _length[1])];
    if (Da == NULL) throw returnR("C++ Error: Could not allocate working memory", 1);

    /* penalty over rows with fixed columns */
    _penalty[0] = 0.0;
    for (col = 0; col < _length[1]; col++){
      for (row = 0; row < _length[0]; row++) Da[row] = _a[col*_length[0] + row];
      diff(_order, _length[0], Da);
      for (row = 0; row < _length[0] - _order; row++) _penalty[0] += (Da[row] * Da[row]);
    }
    _penalty[0] *= (-0.5);

    /* penalty over cols with fixed rowss */
    _penalty[1] = 0.0;
    for (row = 0; row < _length[0]; row++){
      for (col = 0; col < _length[1]; col++) Da[col] = _a[col*_length[0] + row];
      diff(_order, _length[1], Da);
      for (col = 0; col < _length[1] - _order; col++) _penalty[1] += (Da[col] * Da[col]);
    }
    _penalty[1] *= (-0.5);

    if (_equal_lambda)  _penalty[0] += _penalty[1];
    break;

  default:
    throw returnR("C++ Error: Strange _dim in Gspline::penalty_uniCAR", 1);
  }  /** end of switch (_dim)  **/

  delete [] Da;
  return;
}


/***** penalty_eight_neighbors *****/
//
// Penalty for the eight neighbors system
//  computed via sum of weighted squared pairwise differences
//
// penalty ..... one-component pointer
//
void
Gspline::penalty_eight_neighbors() const
{
  int i, j;
  if (_dim != 2) throw returnR("C++ Error: Strange _dim appeares in Gspline::penalty_eight_neighbors", 1);
  int nr = _length[0];     /* this should be always at least 5 */
  int nc = _length[1];     /* this should be always at least 5 */
  _penalty[0] = 0.0;

  // Neighbors of sites (i, 0)
  // ===========================
  /* neighbors of site (0, 0) */
  _penalty[0] += (_a[0]-_a[1])*(_a[0]-_a[1]) + (_a[0]-_a[nr])*(_a[0]-_a[nr]) - (_a[0]-_a[nr+1])*(_a[0]-_a[nr+1]);

  /* neighbors (not yet included) of sites (i, 0), i=1,...,nr-2 */
  for (i = 1; i <= nr-2; i++){
    _penalty[0] += 2*(_a[i]-_a[nr+i])*(_a[i]-_a[nr+i]) + (_a[i]-_a[i+1])*(_a[i]-_a[i+1]) 
                - (_a[i]-_a[nr+i-1])*(_a[i]-_a[nr+i-1]) - (_a[i]-_a[nr+i+1])*(_a[i]-_a[nr+i+1]);
  }

  /* neighbors (not yet included) of site (nr-1, 0) */
  _penalty[0] += (_a[nr-1]-_a[nr+nr-1])*(_a[nr-1]-_a[nr+nr-1]) - (_a[nr-1]-_a[nr+nr-2])*(_a[nr-1]-_a[nr+nr-2]);


  // Neighbors (not yet included) of sites (i, j), j=1,...,nc-2
  // ==========================================================
  for (j = 1; j <= nc-2; j++){
    /* neighbors (not yet included) of site (0, j) */
    _penalty[0] += 2*(_a[j*nr]-_a[j*nr+1])*(_a[j*nr]-_a[j*nr+1]) + (_a[j*nr]-_a[(j+1)*nr])*(_a[j*nr]-_a[(j+1)*nr])
                - (_a[j*nr]-_a[(j+1)*nr+1])*(_a[j*nr]-_a[(j+1)*nr+1]);

    /* neighbors (not yet included) of sites (i, j), i=1,...,nr-2 */
    for (i = 1; i <= nr-2; i++){
      _penalty[0] += 2*((_a[j*nr+i]-_a[(j+1)*nr+i])*(_a[j*nr+i]-_a[(j+1)*nr+i]) + (_a[j*nr+i]-_a[j*nr+i+1])*(_a[j*nr+i]-_a[j*nr+i+1]))
        	  - (_a[j*nr+i]-_a[(j+1)*nr+i-1])*(_a[j*nr+i]-_a[(j+1)*nr+i-1]) - (_a[j*nr+i]-_a[(j+1)*nr+i+1])*(_a[j*nr+i]-_a[(j+1)*nr+i+1]);
    }

    /* neighbors (not yet included) of site (nr-1, j) */
    _penalty[0] += (_a[j*nr+nr-1]-_a[(j+1)*nr+nr-1])*(_a[j*nr+nr-1]-_a[(j+1)*nr+nr-1]) - (_a[j*nr+nr-1]-_a[(j+1)*nr+nr-2])*(_a[j*nr+nr-1]-_a[(j+1)*nr+nr-2]);
  }


  // Neighbors (not yet included) of sites (i, nc-1)
  // ===============================================
  /* neighbors (not yet included) of sites (i, nc-1), i=0,...,nr-2 */
  for (i = 0; i <= nr-2; i++){
    _penalty[0] += (_a[(nc-1)*nr+i]-_a[(nc-1)*nr+i+1])*(_a[(nc-1)*nr+i]-_a[(nc-1)*nr+i+1]);
  }

  _penalty[0] *= (-0.5);
  return;
}


/***** penalty_twelve_neighbors *****/
//
// Penalty for the twelve neighbors system
//
void
Gspline::penalty_twelve_neighbors() const
{
  throw returnR("C++ Error: Not yet implemented (penalty::twelve_neighbors)", 1);
}


/***** diff *****/
//
// Compute ordered differences (diff_in_R is a version which can be called directly from R)
//
// order ..... order of the differences (at least 1)
// na ........ length of the whole vector 'a'
// Da ........  INPUT: vector for which differences should be computed
//              OUTPUT: computed differences (at places 0, ..., la-1-order)
//
void
diff(const int& order, const int& na, double* Da)
{
  int i;
  if (order < 0 || order > na-1) throw returnR("C++ Error: order must be >= 0 & <= length(a) in diff", 1);

  if (order == 0) return;
  else{
    for (i = 1; i < na; i++) Da[i-1] = Da[i] - Da[i-1];
    diff(order-1, na-1, Da);
  }
  return;
}


#ifdef __cplusplus
extern "C" {
#endif
void diff_in_R(int* order, int* na, double* Da)
{
  try{
    int orderR = *order;
    int naR = *na;
    diff(orderR, naR, Da);
    return;
  }
  catch (returnR){
    return;
  }
}
#ifdef __cplusplus
	}
#endif
