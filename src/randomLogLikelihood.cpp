// Function to compute a log-likelihood of random effects
//  in AFT model with random effects being (multivariate) normal
//
// = log(prod_{i=1}^N p(b_i|D, gamma))
//     where gamma is a subvector of beta
//                 is a vector of means of b's (except the random intercept)
//
// 03/02/2004: start working on it
// 13/03/2004: allow for general mean of the random intercept
//
#include "randomLogLikelihood.h"

using namespace std;

// ***** randomLogLikelihood *****
// Compute the value of the randomLogLikelihood
//   either completely from scratch,
//   or update it when values of random effects in several clusters changed
//
// ** function with two prototypes
//
// OUTPUT PARAMETERS:
// 
// randomloglik ........ overall value of the log-likelihood                      (1)
//                       (updated appropriately)
// randomllcl .......... value of the log-likelihood separate for each cluster    (nclusterP)
//        if nclusterP == nclusterUpd, everything is recalculated from scratch,
//          otherwise only appropriate components of randomllcl are updated and using their old values
//          randomloglik is updated
//
// INPUT PARAMETERS:
//
// clusterUpd .......... indeces of clusters for which the log-likelihood contributions are to be updated
//                       -> it is ignored if nclusterUpd = nclusterP
// nclusterUpd ......... number of clusters for which the log-likelihood is to be updated
// nclusterP ........... total number of clusters
// bM .................. values of random effects for all clusters
// betaM ...............
// Dcm .................
// indbinXA ............
// Eb0 ................. mean of the random intercept
//
void
randomLogLikelihood(double* randomloglik,    double* randomllcl,
                    const int* clusterUpd,   const int* nclusterUpd,  const int* nclusterP,
                    const double* bM,        const double* betaM,     const covMatrix* Dcm,
                    const double* Eb0,       const int* indbinXA)
{
  int cl, cll, j;
  int dimb = Dcm->nrow();

  // If D matrix is singular or non-positive definite, assign to all components -FLT_MAX
  //   irrespective of which components the user wants to update
  if (Dcm->rank() < dimb || Dcm->det <= 0){
    *randomloglik = -FLT_MAX;
    for (cl = 0; cl < *nclusterP; cl++) randomllcl[cl] = -FLT_MAX;
    return;
  }

  double* bmingamma = new double[dimb];
  double helpd;
  
  double half_log_det_p_logsqrt2pi = -0.5 * log(Dcm->det) - dimb * LOG_SQRT_2PI;

  if (*nclusterUpd == *nclusterP){
    *randomloglik = 0.0;                              // reset it completely
    for (cl = 0; cl < *nclusterP; cl++){
      for (j = 0; j < dimb; j++){     
        bmingamma[j] = bM[dimb*cl + j] - (indbinXA[j] < 0 ? (*Eb0) : betaM[indbinXA[j]]);
      }
      randomllcl[cl] = half_log_det_p_logsqrt2pi; 
      axMxa(&helpd, bmingamma, Dcm->icovm, &ONE_INT, &dimb, &dimb, Dcm->diagI);
      randomllcl[cl] += (-0.5)*helpd;
      *randomloglik += randomllcl[cl];
    }
  }

  else{
    for (cll = 0; cll < *nclusterUpd; cll++){
      cl = clusterUpd[cll];
      *randomloglik -= randomllcl[cl];                 // subtract the old value
      for (j = 0; j < dimb; j++){     
        bmingamma[j] = bM[dimb*cl + j] - (indbinXA[j] < 0 ? 0.0 : betaM[indbinXA[j]]);
      }
      randomllcl[cl] = half_log_det_p_logsqrt2pi; 
      axMxa(&helpd, bmingamma, Dcm->icovm, &ONE_INT, &dimb, &dimb, Dcm->diagI);
      randomllcl[cl] += (-0.5)*helpd;
      *randomloglik += randomllcl[cl];                 // add the new value
    }
  }  

  delete [] bmingamma;
  return;

}    // end of the function randomLogLikelihood



// ***** randomLogLikelihood *****
// Update the value of the randomLogLikelihood
//   when values of the random effect changed for only one cluster
//
// OUTPUT PARAMETERS:
// 
// randomloglik ........ overall value of the log-likelihood                      (1)
//                       (updated appropriately)
// randomllcl .......... value of the log-likelihood separate for each cluster    (nclusterP)
//              * only one component is updated
//
// INPUT PARAMETERS:
//
// clusterUpd .......... index of the cluster for which the random effects changed (1)
// nclusterP ........... number of all clusters                    
// bclM ................ values of random effects for the cluster that changed
// betaM ...............
// Dcm .................
// indbinXA ............
//
void
randomLogLikelihood(double* randomloglik,   double* randomllcl,
                    const int* clusterUpd,  const int* nclusterP, 
                    const double* bclM,     const double* betaM,    const covMatrix* Dcm,
                    const double* Eb0,      const int* indbinXA)
{
  int cl, j;
  int dimb = Dcm->nrow();

  // If D matrix is singular or non-positive definite, assign to all components -FLT_MAX
  //   irrespective of which components the user wants to update
  if (Dcm->rank() < dimb || Dcm->det <= 0){
    *randomloglik = -FLT_MAX;
    for (cl = 0; cl < *nclusterP; cl++) randomllcl[cl] = -FLT_MAX;
    return;
  }

  double* bmingamma = new double[dimb];
  double helpd;
  
  double half_log_det_p_logsqrt2pi = -0.5 * log(Dcm->det) - dimb * LOG_SQRT_2PI;

  cl = *clusterUpd;
  *randomloglik -= randomllcl[cl];                 // subtract the old value
  for (j = 0; j < dimb; j++){     
    bmingamma[j] = bclM[j] - (indbinXA[j] < 0 ? (*Eb0) : betaM[indbinXA[j]]);
  }
  randomllcl[cl] = half_log_det_p_logsqrt2pi; 
  axMxa(&helpd, bmingamma, Dcm->icovm, &ONE_INT, &dimb, &dimb, Dcm->diagI);
  randomllcl[cl] += (-0.5)*helpd;
  *randomloglik += randomllcl[cl];                 // add the new value

  delete [] bmingamma;
  return;
}

