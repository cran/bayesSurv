/*** structRandomEff32.cpp ***/
//
#ifndef _STRUCT_RANDOM_EFF_32_H_
#define _STRUCT_RANDOM_EFF_32_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BLAS_LAPACK.h"
#include "Gspline.h"
#include "Mvtdist3.h"

namespace RandomEff32 {

  /*** Pointers _nwithinCl, _d, _b  point to external data  ***/
typedef struct re {
  int _nRandom;          // dimension of the random effect                                   (init makes it 2)
  int _lD;               // length of the lower triangle of the matrix _nRandom x _nRandom 

  int  _nCluster;            // number of clusters  
  const int *_nwithinCl;     // [_nCluster] number of observations within each cluster

  double *_d;          // [_nCluster] random intercepts for onset                                 (pointer to external data)
  double *_b;          // [_nCluster] random intercepts for time-to-event                         (pointer to external data)

  double *_D;          // lower triangle of the covariance matrix of random effects               (pointer to external data)
  double _Di[3];       // lower triangle of the inverse covariance matrix of the random effects
  double _detD;        // determinant of D

  double _priorDF;     // prior degrees of freedom for the Wishart prior
  double _priorSi[3];  // lower triangle of the prior inverse scale matrix of the Wishart prior

  double _propVar[3];    // working space to store (inverse) variance of the full conditional distribution
  double _propMean[2];   // working space to store mean of the full conditional distributions
  double _propValue[2];  // working space to store the proposed value of the random effects

  double _propDF;          // degrees of freedom of the Wishart full conditional distribution
  double _propSi[3];       // working space to store the (inverse) scale matrix of the Wishart full conditional distribution
  double _workWishart[8];  // working space for Mvtnorm3::rwishart
} RE;

void
init(RandomEff32::RE *data,  double *dVal,   double *bVal,  double *parD,  const int *pardI,  const int *parbI);

void
update(RandomEff32::RE *data,  
       double *regResOnset,     double *regResTime,
       const int *nP,
       const Gspline *gg_zeta,  double** const mu_zeta,    const int *rM_zeta,
       const Gspline *gg_eps,   double** const mu_eps,     const int *rM_eps);

void
updateAfterChangeD(RandomEff32::RE *data);

void
predict_db(RandomEff32::RE *data);

}  /** end of the namespace RandomEff32 **/

#endif
