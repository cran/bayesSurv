// Function to update the allocations
//   Gibbs step is used

// 27/11/2003: start woking on it
// 13/03/2004: general mean of the random intercept allowed
// 19/03/2004: updateAlloc that does not update invrM and mixtureNM was added
//             (to be used by 'predictive' functions)
//

#include "bayessurvreg.h"

using namespace std;

// ======================================================================
//
// OUTPUT PARAMETERS:
//
// rM ...........
// invrM ........
// mixtureNM ....
//
// INPUT PARAMETERS:
//
// wM ........... (at least kP x 1)
// muM .......... (at least kP x 1)
// invsigma2M ... (at least kP x 1)
// kP ........... curent number of mixture components
// regresResM ... (nP x 1)
// Eb0 .......... (1)
// randomIntP ... (1)
// nP ........... number of observations
//
void
updateAlloc(int* rM,                List<int>* invrM,          int* mixtureNM, 
            const double* wM,       const double* muM,         const double* invsigma2M,
            const int* kP,          const double* regresResM,  const double* Eb0,
            const int* randomIntP,  const int* nP)
{
  int obs, j;
  double intcptadd;
  if (*kP == 1){
    for (obs = 0; obs < *nP; obs++) rM[obs] = 0;
    return;
  }

  double u;

  if (*randomIntP) intcptadd = *Eb0;
  else             intcptadd = 0.0;  

  // Reset the lists of observations belonging to the components (sufficient first *kP lists)
  //  and reset numbers of observations belonging to the mixture components
  for(j = 0; j < *kP; j++){
    invrM[j] = List<int>();
    mixtureNM[j] = 0;
  }

  // Update the allocations
  double*  winvsigmaM = new double[*kP];
  double* probsM = new double[*kP];          // probabilities of non-zero components only
  double* cumprobsM = new double[*kP];
  int* origInd = new int[*kP];               // correspondence between indeces of non-zero components and all indeces
  for (j = 0; j < *kP; j++) winvsigmaM[j] = wM[j] * sqrt(invsigma2M[j]);

  for (obs = 0; obs < *nP; obs++){   
    probsM[0] = winvsigmaM[0] * exp(-0.5 * invsigma2M[0] * (regresResM[obs] - muM[0] + intcptadd)*(regresResM[obs] - muM[0] + intcptadd));
    cumprobsM[0] = probsM[0];
    for (j = 1; j < *kP; j++){
      probsM[j] = winvsigmaM[j] * exp(-0.5 * invsigma2M[j] * (regresResM[obs] - muM[j] + intcptadd)*(regresResM[obs] - muM[j] + intcptadd));
      cumprobsM[j] = cumprobsM[j - 1] + probsM[j];
    }   
    discreteSampler(rM + obs, cumprobsM, kP, &ONE_INT, &ONE_INT, &ZERO_INT);
    invrM[rM[obs]].addNode(obs);
    mixtureNM[rM[obs]]++;
  }

  delete [] winvsigmaM;
  delete [] probsM;
  delete [] cumprobsM;
  delete [] origInd;
  return;
}   // end of function updateAlloc


// ======================================================================
// second prototype:
//
void
updateAlloc(int* rM,                
            const double* wM,       const double* muM,         const double* invsigma2M,
            const int* kP,          const double* regresResM,  const double* Eb0,
            const int* randomIntP,  const int* nP)
{
  int obs, j;
  double intcptadd;
  if (*kP == 1){
    for (obs = 0; obs < *nP; obs++) rM[obs] = 0;
    return;
  }

  double u;

  if (*randomIntP) intcptadd = *Eb0;
  else             intcptadd = 0.0;  

  // Update the allocations
  double*  winvsigmaM = new double[*kP];
  double* probsM = new double[*kP];          // probabilities of non-zero components only
  double* cumprobsM = new double[*kP];
  int* origInd = new int[*kP];               // correspondence between indeces of non-zero components and all indeces
  for (j = 0; j < *kP; j++) winvsigmaM[j] = wM[j] * sqrt(invsigma2M[j]);

  for (obs = 0; obs < *nP; obs++){   
    probsM[0] = winvsigmaM[0] * exp(-0.5 * invsigma2M[0] * (regresResM[obs] - muM[0] + intcptadd)*(regresResM[obs] - muM[0] + intcptadd));
    cumprobsM[0] = probsM[0];
    for (j = 1; j < *kP; j++){
      probsM[j] = winvsigmaM[j] * exp(-0.5 * invsigma2M[j] * (regresResM[obs] - muM[j] + intcptadd)*(regresResM[obs] - muM[j] + intcptadd));
      cumprobsM[j] = cumprobsM[j - 1] + probsM[j];
    }   
    discreteSampler(rM + obs, cumprobsM, kP, &ONE_INT, &ONE_INT, &ZERO_INT);
  }

  delete [] winvsigmaM;
  delete [] probsM;
  delete [] cumprobsM;
  delete [] origInd;
  return;
}   // end of function updateAlloc

