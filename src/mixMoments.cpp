// Function to compute mixture mean and standard deviation
//
// 13/03/2004
//
#include "mixMoments.h"

using namespace std;

//
// PARAMETERS:
// 
// mixMomentM ........ mixMomentM[0] = mixture mean
//                     mixMomentM[1] = mixture standard deviation
// onlySD ............ recalculate only standard deviation
//
void
mixMoments(double* mixMomentM,  const int* kP,             const double* wM,
           const double* muM,   const double* invsigma2M,  const bool onlySD)
{
  int j;
  double sigma2;

  if (!onlySD) mixMean(mixMomentM + 0, kP, wM, muM);

  mixMomentM[1] = 0.0;
  for (j = 0; j < *kP; j++){
    if (invsigma2M[j] > 0.0) sigma2 = 1/invsigma2M[j];
    else{                    mixMomentM[1] = FLT_MAX;  return;}
    mixMomentM[1] += wM[j] * (muM[j]*muM[j] + sigma2);
  }
  mixMomentM[1] -= mixMomentM[0] * mixMomentM[0];
  if (mixMomentM[1] >= 0.0) mixMomentM[1] = sqrt(mixMomentM[1]);
  else                      mixMomentM[1] = 0.0;

  return;
}


void
mixMean(double* mixMeanP,  const int* kP,  const double* wM,  const double* muM)
{
  int j;

  *mixMeanP = wM[0]*muM[0];
  for (j = 1; j < *kP; j++) *mixMeanP += wM[j]*muM[j];  

  return;
}
