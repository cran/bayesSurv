// Functions to update quantiles and cumulative sums
//  and to compute means at the end

// Only cumulative sums are updated, quantiles cannot be computed adaptively 
//   from the beginning (first I thought that it would be possible...)

//  (used primarily in 'predictive' function)

// 23/03/2004

#include "bayessurvreg.h"

using namespace std;

// =================================================
// cumsumQuantile1:
// =================================================
//
// PARAMETERS:
//
// quant[nobs][nquant + 1] ...... array with current quantiles,
//                                quant[*][nquant] = cumulative sums
// newval[nobs][iter] ........... last sampled values
// nquant ....................... number of quantiles to be computed
// nobs ......................... number of observartions/parameters for which
//                                quantiles are to be computed
// iter ......................... index of last sampled value
// 
void
cumsumQuantile1(double** quant,  double** newval,  const int nquant,  const int nobs,  
                const int iter)
{
  int obs, i;
  for (obs = 0; obs < nobs; obs++){
    quant[obs][nquant] += newval[obs][iter];
  }
  return;
}


// =================================================
// cumsumQuantile2:
// =================================================
//
// PARAMETERS:
//
// quant[nobs][ngridM[obs]][nquant + 1] ..... array with current quantiles
//                                            quant[*][*][nquant] = cumulative sums
//
// newval[nobs][ngridM[obs]][iter]........... last sampled values 
// nquant ........................ number of quantiles to be computed
// nobs .......................... number of observartions/parameters for which
//                                 quantiles are to be computed
// ngridM ........................ number of grid values for each observation/parameter
// iter .......................... index of last sampled value
// 
void
cumsumQuantile2(double*** quant,  double*** newval,     const int nquant,  const int nobs,  const int* ngridM,  
                const int iter)
{
  int obs, i, j;
  for (obs = 0; obs < nobs; obs++){
    for (j = 0; j < ngridM[obs]; j++){
      quant[obs][j][nquant] += newval[obs][j][iter];
    }
  }
  return;
}


// ====================================================
// meanQuantile1:
// ====================================================
//
// PARAMETERS:
//
// quant[nobs][nquant] ......... on INPUT: array for quantiles and cumulative sums
//                               on OUTPUT: quantiles in quant[obs][0,...,nquant-1]
//                                          means in quant[obs][nquant]
// sample[nobs][sampleSize] ... sampled values
// quantileA[nquant] ........................ quantiles to be computed
// indquant1[nquant]......................... first indeces of observations in the whole chain 
//                                            corresponding to desired quantiles
// indquant2[nquant] ........................ second indeces of observations in the whole chain 
//                                            corresponding to desired quantiles
// nobs ........... length of quant
// nquant ......... length of each quant[obs]
// sampleSize ..... sample size
//
void
meanQuantile1(double** quant,           double** sample,   
              const double* quantileA,  const int* indquant1,  const int* indquant2,
              const int nobs,           const int nquant,      const int sampleSize)
{
  int obs, i;
  if (sampleSize <= 0) throw returnR("C++ Error: sample size = 0 when computing empirical mean.", 1);

  for (obs = 0; obs < nobs; obs++){
    quant[obs][nquant] /= sampleSize;
    sort(sample[obs], sample[obs] + sampleSize);
    for (i = 0; i < nquant; i++) 
      quant[obs][i] = quantileA[i]*sample[obs][indquant1[i]] + (1 - quantileA[i])*sample[obs][indquant2[i]];    
  }
 
  return;
}


// ====================================================
// meanQuantile2:
// ====================================================
//
// PARAMETERS:
//
// quant[nobs][ngridM[obs]][nquant] ......... on INPUT: array for quantiles and cumulative sums
//                                            on OUTPUT: quantiles in quant[obs][grid.val][0,...,nquant-1]
//                                                       means in quant[obs][grid.val][nquant]
// sample[nobs][ngridM[obs]][sampleSize] ... sampled values
// quantileA[nquant] ........................ quantiles to be computed
// indquant1[nquant]......................... first indeces of observations in the whole chain 
//                                            corresponding to desired quantiles
// indquant2[nquant] ........................ second indeces of observations in the whole chain 
//                                            corresponding to desired quantiles
// nobs ........... length of quant
// ngridM[nobs] ... length of quant[obs]
// nquant ......... length of each quant[obs][grid.val]
// sampleSize ..... sample size
//
void
meanQuantile2(double*** quant,          double*** sample,  
              const double* quantileA,  const int* indquant1,  const int* indquant2,
              const int nobs,           const int* ngridM,     const int nquant,      const int sampleSize)
{
  int obs, i, j;
  if (sampleSize <= 0) throw returnR("C++ Error: sample size = 0 when computing empirical mean.", 1);

  for (obs = 0; obs < nobs; obs++){
    for (j = 0; j < ngridM[obs]; j++){    
      quant[obs][j][nquant] /= sampleSize;
      sort(sample[obs][j], sample[obs][j] + sampleSize);
      for (i = 0; i < nquant; i++) 
        quant[obs][j][i] = quantileA[i]*sample[obs][j][indquant1[i]] + (1 - quantileA[i])*sample[obs][j][indquant2[i]];    
    }
  }
 
  return;
}


