// ************************************************************************************** //
// ********* Gspline_update_gamma.cpp *************************************************** //
// ************************************************************************************** //
// Member function of class Gspline to update overall intercept 'gamma'
//   * for the first specification with variable knots that depend on the basis standard deviation
//
// 19/10/2004: 'Gspline::update_gamma'
// 10/12/2004: overall intercept and scale of the G-spline (second specification) allowed
//
#include "Gspline.h"

// ***** Gspline::update_gamma *****
//
// regresResM[nP*_dim] ...... * in the case of i.i.d data, these are observations y_1,...,y_n,
//                            * if dim(y_i) = d > 1 then it is assumed that they are stored as y_1[1],...,y_1[d], ..., y_n[1], ..., y_n[d],
//                              i.e. y_i[j] = regresResM[i*d + j]
//                            * in the case of regression, these are y_i - beta'x_i - b_i'z_i
// rM[nP] ................... component labels taking values 0, ..., total_length-1
//                            * sorted in the same way as regresResM
// nP ....................... total number of 'independent' observational vectors
//
void
Gspline::update_gamma(const double* regresResM, const int* rM, const int* nP)
{
  static int j, obs;
  static double fullmean[_max_dim];            // means of the full conditional
  static double fullscale[_max_dim];           // for scale of the full conditional
  const double* y_obs;
  const int* rp;

  y_obs = regresResM;
  rp = rM;
  switch (_dim){
  case 1:
    switch (_prior_for_gamma[0]){
    case _Fixed_:
      break;

    case Normal:
      fullscale[0] = 1 / ((*nP) * (_invsigma2[0]) + _prior_gamma[1]);                  // = var(gamma|...)
      fullmean[0] = 0.0;
      for (obs = 0; obs < *nP; obs++){
        fullmean[0] += (*y_obs - _intcpt[0])/_scale[0] - ((*rp) - _K[0])*_delta[0];
        y_obs++;
        rp++;
      }
      fullmean[0] *= _invsigma2[0];
      fullmean[0] += _prior_gamma[0] * _prior_gamma[1];
      fullmean[0] *= fullscale[0];
      fullscale[0] = sqrt(fullscale[0]);
      _gamma[0] = rnorm(fullmean[0], fullscale[0]);
      break;

    default:     
      throw returnR("C++ Error: Unimplemented prior for gamma appeared in Gspline::update_gamma", 1);
    }
    break;

  case 2:
    for (j = 0; j < _dim; j++) fullscale[j] = 1 / ((*nP) * (_invsigma2[j]) + _prior_gamma[2*j+1]);                  // = var(gamma|...)
    fullmean[0] = fullmean[1] = 0.0;
    for (obs = 0; obs < *nP; obs++){
      fullmean[0] += (*y_obs - _intcpt[0])/_scale[0] - (((*rp) % _length[0]) - _K[0])*_delta[0];
      y_obs++;

      fullmean[1] += (*y_obs - _intcpt[1])/_scale[1] - (((*rp) / _length[0]) - _K[1])*_delta[1];
      y_obs++;
      rp++;
    }
    for (j = 0; j < _dim; j++){
      switch (_prior_for_gamma[j]){
      case _Fixed_:
        break;

      case Normal:
        fullmean[j] *= _invsigma2[j];
        fullmean[j] += _prior_gamma[2*j] * _prior_gamma[2*j+1];
        fullmean[j] *= fullscale[j];
        fullscale[j] = sqrt(fullscale[j]);
        _gamma[j] = rnorm(fullmean[j], fullscale[j]);
        break;

      default:     
        throw returnR("C++ Error: Unimplemented prior for gamma appeared in Gspline::update_gamma", 1);
      }
    } 
    break;

  default:
    throw returnR("C++ Error: Gspline::update_gamma not implemented for _dim > 2", 1);
  }

  return;
}
