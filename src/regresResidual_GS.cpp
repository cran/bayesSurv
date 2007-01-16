// Various functions to compute regression residuals 
//  * used by functions that work with G-splines
//
// 12/01/2005: 'regresRes_GS'
// 12/12/2006: 'regresRes_GS2006'
// 24/01/2005: 'linPred_GS'

#include "regresResidual_GS.h"

// ***** regresRes_GS *****
//
// regResid[nobs] ................. computed regression residuals = Y - x'beta - z'b
// YsA[nobs] ...................... observations
// bg ............................. object holding information on regression part of the model
// bA[ncluster*bg->nRandom()] ..... current values of random effects
// XA[nobs*bg->nbeta()] ........... design matrix stored in an array in ROW major order
//                                  i.e. XA[0],...,XA[nbeta-1] is a covariate vector for the first observation etc.
// nincluster[ncluster] ........... numbers of observations per cluster (observations are assumed to be sorted according to clusters)
// nobs ........................... number of observations
// ncluster ....................... number of clusters
//
void
regresRes_GS(double* regResid,  const double* YsA,      const BetaGamma* bg,  const double* bA,  
             const double* XA,  const int* nincluster,  const int* nobs,      const int* ncluster)
{

  int i, j, k;
  double* regRes = regResid;
  const double* yobs   = YsA;
  const double* bobs   = bA;
  const double* xrow   = XA;

  /*** No regression  ***/
  if (!bg->nFixed() && !bg->nRandom()){
    for (i = 0; i < *nobs; i++){
      *regRes = *yobs;
      regRes++;
      yobs++;
    }
    return;
  }

  /*** Only fixed effects, no random effects ***/
  if (!bg->nRandom()){
    for (i = 0; i < *nobs; i++){
      *regRes = *yobs;
      for (j = 0; j < bg->nbeta(); j++){
        *regRes -= xrow[j] * bg->beta(j);
      }
      regRes++;
      yobs++;
      xrow += bg->nbeta();
    }
    return;
  }

  /*** There are some random effects ***/
  for (i = 0; i < *ncluster; i++){
    for (j = 0; j < nincluster[i]; j++){
      *regRes = *yobs;
      if (bg->randomIntcpt()) *regRes -= bobs[0];
      for (k = 0; k < bg->nbeta(); k++){
        if (bg->indbA(k) == -1) *regRes -= xrow[k] * bg->beta(k);
        else                    *regRes -= xrow[k] * bobs[bg->indbA(k)];
      }
      regRes++;
      yobs++;
      xrow += bg->nbeta();
    }
    bobs += bg->nRandom();
  }
  return;
}


// ***** regresRes_GS2006 *****
//
// Improved version of regresRes
//
// regResid[nobs] ................. computed regression residuals = Y - x'beta - z'b
// YsA[nobs] ...................... observations
// bg ............................. object holding information on regression part of the model
// bA[ncluster*bg->nRandom()] ..... current values of random effects
// XA[nobs*bg->nbeta()] ........... design matrix stored in an array in ROW major order
//                                  i.e. XA[0],...,XA[nbeta-1] is a covariate vector for the first observation etc.
// nincluster[ncluster] ........... numbers of observations per cluster (observations are assumed to be sorted according to clusters)
// nobs ........................... number of observations
// ncluster ....................... number of clusters
//
void
regresRes_GS2006(double* regResid,  const double* YsA,      const BetaGamma* bg,  const double* bA,  
                 const double* XA,  const int* nincluster,  const int* nobs,      const int* ncluster)
{

  int i, j, k;
  double* regRes = regResid;
  const double *yobs   = YsA;
  const double *bobs   = bA;
  const double *xP     = XA;
  const double *beta;
  const int *nClP      = nincluster;

  /*** No regression  ***/
  if (!bg->nFixed() && !bg->nRandom()){
    for (i = 0; i < *nobs; i++){
      *regRes = *yobs;
      regRes++;
      yobs++;
    }
    return;
  }

  /*** Only fixed effects, no random effects ***/
  if (!bg->nRandom()){
    for (i = 0; i < *nobs; i++){
      beta    = bg->betaP();
      *regRes = *yobs;      
      for (j = 0; j < bg->nbeta(); j++){
        *regRes -= (*xP) * (*beta);
        xP++;
        beta++;
      }
      regRes++;
      yobs++;
    }
    return;
  }

  /*** There are some random effects ***/
  for (i = 0; i < *ncluster; i++){
    for (j = 0; j < *nClP; j++){
      beta    = bg->betaP();
      *regRes = *yobs;
      if (bg->randomIntcpt()) *regRes -= bobs[0];
      for (k = 0; k < bg->nbeta(); k++){
        if (bg->indbA(k) == -1) *regRes -= (*xP) * (*beta);
        else                    *regRes -= (*xP) * bobs[bg->indbA(k)];
        xP++;
        beta++;
      }
      regRes++;
      yobs++;
    }
    bobs += bg->nRandom();
    nClP++;
  }
  return;
}


// ***** linPred_GS *****
//
// * used by 'predictive_GS'
//
// linPred[nobs] .................. computed linear predictors = x'beta + z'b
// bg ............................. object holding information on regression part of the model
// bA[ncluster*bg->nRandom()] ..... current values of random effects
// XA[nobs*bg->nbeta()] ........... design matrix stored in an array in ROW major order
//                                  i.e. XA[0],...,XA[nbeta-1] is a covariate vector for the first observation etc.
// nincluster[ncluster] ........... numbers of observations per cluster (observations are assumed to be sorted according to clusters)
// nobs ........................... number of observations
// ncluster ....................... number of clusters
//
void
linPred_GS(double* linPred,   const BetaGamma* bg,    const double* bA,  
           const double* XA,  const int* nincluster,  const int* nobs,      const int* ncluster)
{

  int i, j, k;
  double* lin_pred = linPred;
  const double* bobs   = bA;
  const double* xrow   = XA;

  /*** No regression  ***/
  if (!bg->nFixed() && !bg->nRandom()){
    for (i = 0; i < *nobs; i++){
      *lin_pred = 0.0;
      lin_pred++;
    }
    return;
  }

  /*** Only fixed effects, no random effects ***/
  if (!bg->nRandom()){
    for (i = 0; i < *nobs; i++){
      *lin_pred = 0.0;
      for (j = 0; j < bg->nbeta(); j++){
        *lin_pred += xrow[j] * bg->beta(j);
      }
      lin_pred++;
      xrow += bg->nbeta();
    }
    return;
  }

  /*** There are some random effects ***/
  for (i = 0; i < *ncluster; i++){
    for (j = 0; j < nincluster[i]; j++){
      *lin_pred = 0.0;
      if (bg->randomIntcpt()) *lin_pred += bobs[0];
      for (k = 0; k < bg->nbeta(); k++){
        if (bg->indbA(k) == -1) *lin_pred += xrow[k] * bg->beta(k);
        else                    *lin_pred += xrow[k] * bobs[bg->indbA(k)];
      }
      lin_pred++;
      xrow += bg->nbeta();
    }
    bobs += bg->nRandom();
  }
  return;
}




