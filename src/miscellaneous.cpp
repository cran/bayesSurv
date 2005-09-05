// Smaller helping function for the Bayesian AFT model

// 24/11/2003: start working on it

#include "miscellaneous.h"

using namespace std;

// ================================================================================
// Function to compute numbers of observations belonging to each mixture component
// ** overloaded function with two prototypes
// ================================================================================
void
giveMixtureN(int* mixtureNM,                   // result (at least *kP x 1)
             const int* kP,                    // number of mixture components
             const int* rM,                    // vector of pertainence indicators  (*nP x 1)
             const int* nP                     // number of observations             
             )
{
  for (int j = 0; j < *kP; j++) mixtureNM[j] = 0;
  for (int obs = 0; obs < *nP; obs++){
      mixtureNM[rM[obs]]++;
  }

  return;
}   // end of function giveMixtureN


// ================================================================================
void
giveMixtureN(int* mixtureNM,           // result (at least *kP x 1)
             const int* kP,            // number of mixture components
             const List<int>* invrM    // pointers to lists indicating observations belonging to different mixture components
             )
{
  for (int j = 0; j < *kP; j++) mixtureNM[j] = invrM[j].length();

  return;
}   // end of function giveMixtureN


// ================================================================================
// invsigma2 -> sigma
// ================================================================================
void
giveSigmaAndInvsigma2(double* sigmaM,  double* invsigma2M,  const double* sigma2M,  const int* kP)
{
  int j;

  for (j = 0; j < *kP; j++){
    if (sigma2M[j] <= 0.0){
      sigmaM[j] = 0.0;
      invsigma2M[j] = FLT_MAX;
    }
    else{
      sigmaM[j] = sqrt(sigma2M[j]);
      invsigma2M[j] = 1 / sigma2M[j];
    }
  }
  
  return;
}


// ================================================================================
// Function to compute regression residuals
//  r = y - x'beta - z'b
//
// ** overloaded function with six prototypes
//
// ================================================================================
void
regresResidual(double* regresResA,
               const double* YsA,    const double* betaA,     const double* bA,
               const double* XA,     const int* clusteriA,    const int* randomIntP,   
               const int* indbA,     const int* nP,           const int* nXP,         const int* nrandomP)
{
  int i, j;

  for (i = 0; i < *nP; i++){
    regresResA[i] = YsA[i];
    if (*randomIntP) regresResA[i] -= bA[(*nrandomP)*clusteriA[i]];     // subtract a random intercept
    for (j = 0; j < *nXP; j++){
      if (indbA[j] == -1)
        regresResA[i] -= XA[(*nP)*j + i] * betaA[j];
      else
        regresResA[i] -= XA[(*nP)*j + i] * bA[(*nrandomP)*clusteriA[i] + indbA[j]];
    }
  }

  return;
}


// ================================================================================
// ***** regresPredictor *****
//
// Function to compute (x'beta + z'b)
//
// ================================================================================
void 
regresPredictor(double* regresPredA,
                const double* betaA,  const double* bA,
                const double* XA,     const int* clusteriA,    const int* randomIntP,   
                const int* indbA,     const int* nP,           const int* nXP,         const int* nrandomP)
{
  int i, j;

  for (i = 0; i < *nP; i++){
    regresPredA[i] = 0.0;
    if (*randomIntP) regresPredA[i] += bA[(*nrandomP)*clusteriA[i]];     // subtract a random intercept
    for (j = 0; j < *nXP; j++){
      if (indbA[j] == -1)
        regresPredA[i] += XA[(*nP)*j + i] * betaA[j];
      else
        regresPredA[i] += XA[(*nP)*j + i] * bA[(*nrandomP)*clusteriA[i] + indbA[j]];
    }
  }

  return;
}



// ================================================================================
// Function to update regression residuals when y changes
//  r = y - x'beta - z'b
// ================================================================================
//
// YsA ...... old value of Y
// newYsA ... new value of Y
//
void
regresResidual(double* regresResA,
               const double* YsA,   const double* newYsA,
               const int* nP)
{
  int i;

  for (i = 0; i < *nP; i++){
    regresResA[i] -= (YsA[i] - newYsA[i]);
  }

  return;
}


// ================================================================================
// Function to update regression residuals when beta changes
// NOTE: changes in beta's that correspond to random effects do not have any impact
//       to regression residuals
//  r = y - x'beta - z'b
// ================================================================================
//
// betaA ...... array of ALL old betas
// newbetaA ... array of ALL new betas ONLY                     
// indnewA .... indeces of new betas in the whole beta vector   length = nnewP
//
void
regresResidual(double* regresResA,
               const double* betaA,   const double* newbetaA,   const int* indnewA,  
               const int* nnewP,      const double* XA,         const int* indbA,         
               const int* nP)
{
  int i, j;

  for (i = 0; i < *nP; i++){
    for (j = 0; j < *nnewP; j++){
      if (indbA[indnewA[j]] == -1) regresResA[i] += XA[(*nP)*indnewA[j] + i] * (betaA[indnewA[j]] - newbetaA[indnewA[j]]);
    }
  }

  return;
}


// ================================================================================
// Function to update regression residuals when b (random effects) changes
//  * change of same random effects for all clusters
//  r = y - x'beta - z'b
// ================================================================================
//
// bA ........ array of ALL old values of b                 length = (nclusterP * nrandomP)
// newbA ..... array of NEW values of b ONLY                length = (nclusterP * nnewP)
// indnewA ... indeces of new b's in the whole sequence
//
void
regresResidual(double* regresResA,
               const double* bA,      const double* newbA,  const int* indnewA,
               const int* nnewP,      const double* XA,     const int* clusteriA,    const int* randomIntP,   
               const int* indbinXA,   const int* nP,        const int* nXP,          const int* nrandomP)
{
  int i, j;

  int start = 0;
  for (i = 0; i < *nP; i++){
    if (*randomIntP && (indnewA[0] == 0)){ 
      start = 1;
      regresResA[i] += (bA[(*nrandomP)*clusteriA[i]] - newbA[(*nnewP)*clusteriA[i]]);
    }
    for (j = start; j < *nnewP; j++){
      regresResA[i] += XA[(*nP)*indbinXA[indnewA[j]] + i] * (bA[(*nrandomP)*clusteriA[i] + indnewA[j]] - newbA[(*nnewP)*clusteriA[i] + j]);
    }
  }

  return;
}


// ===================================================================================================
// Function to update regression residuals when all b (random effects) changes for a specific cluster
//  r = y - x'beta - z'b
// ===================================================================================================
//
// bA ....... array of ALL old values of b                          (nrandomP x nclusterP)
// bclA ..... array of new values of b for a specific cluster only  (nrandomP)
// cl ....... index of a cluster for which b changed                (1)
// indobs ... indeces of observations belonging to a cluster cl     
//
void
regresResidual(double* regresResA,
               const double* bA,         const double* bclA,   const int* cl,
               const List<int>* indobs,  const double* XA,     const int* randomIntP,
               const int* indbinXA,      const int* nP,        const int* nXP,         const int* nrandomP)
{
  int nobs = indobs->length();
  int i, j, obs;

  for (i = 0; i < nobs; i++){
    obs = (*indobs)[i];
    if (*randomIntP) 
      regresResA[obs] += (bA[(*nrandomP)*(*cl)] - bclA[0]);
    for (j = *randomIntP; j < *nrandomP; j++){
      regresResA[obs] += XA[(*nP)*indbinXA[j] + obs] * (bA[(*nrandomP)*(*cl) + j] - bclA[j]);
    }
  }

  return;
}

// ===================================================================================================
// Function to update regression residuals when some b (random effects) changes for a specific cluster
//  r = y - x'beta - z'b
// ===================================================================================================
//
// bA ....... array of ALL old values of b                             (nrandomP x nclusterP)
// bclA ..... array of BOTH new values of b and unchanged values of b  
//            for a specifiv cluster                                   (nrandomP)
// indnewA .. indeces of changed b's                                   (nnewP)
// nnewP .... number of changed b's                                    (1)
// cl ....... index of a cluster for which b changed                   (1)
// indobs ... indeces of observations belonging to a cluster cl     
//
void
regresResidual(double* regresResA,
               const double* bA,         const double* bclA,   
               const int* indnewA,       const int* nnewP,     const int* cl,
               const List<int>* indobs,  const double* XA,     const int* randomIntP,
               const int* indbinXA,      const int* nP,        const int* nXP,         const int* nrandomP)
{
  int nobs = indobs->length();
  int i, j, obs;

  int start = 0;
  for (i = 0; i < nobs; i++){
    obs = (*indobs)[i];
    if (*randomIntP && (indnewA[0] == 0)){ 
      start = 1;
      regresResA[obs] += (bA[(*nrandomP)*(*cl)] - bclA[0]);
    }
    for (j = start; j < *nnewP; j++){
      regresResA[obs] += XA[(*nP)*indbinXA[indnewA[j]] + obs] * (bA[(*nrandomP)*(*cl) + indnewA[j]] - bclA[indnewA[j]]);
    }
  }

  return;
}

