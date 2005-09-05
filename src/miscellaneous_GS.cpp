// Some helping functions for a Bayesian estimation of the density
//  using the G-splines
//
// 23/10/2004: 'findClosestKnot'
//
#include "miscellaneous_GS.h"

extern "C"{

// ====================================================================================
// findClosestKnot: For each observation, find the index of the closest knot which is 
//                  smaller than (or equal to) the observation
//                  (except when the observation is smaller than the first knot)
//           * to be used in bayesHistogram.priorInit function in R
// ====================================================================================
//
// index[nobs] ....... computed indeces of closest knots (on scale 0, ..., nknot-1)
// knot[nknot] ....... SORTED knots
// obs[nobs] ......... observations
//
void
findClosestKnot(int* index,  const double* knot,  const double* obs,  const int* nknot,  const int* nobs)
{
  int i, left, ind, right;
  for (i = 0; i < *nobs; i++){
    if (obs[i] >= knot[*nknot - 1]){
      index[i] = *nknot - 1;
      continue;
    }
    if (obs[i] <= knot[0]){
      index[i] = 0;
      continue;
    }
    left = 0;
    right = *nknot - 1;
    while (right-left > 1){
      ind = (right + left)/2;
      if (obs[i] < knot[ind]) right = ind;
      else                    left = ind;
    }
    index[i] = left;
  }
  return;
}  /** end of function find closest knots **/

}  /** end of extern "C" **/
