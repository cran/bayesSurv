// Functions to update quantiles and cumulative sums
//  and to compute means at the end

// Only cumulative sums are updated, quantiles cannot be computed adaptively 
//   from the beginning (first I thought that it would be possible...)

//  (used primarily in 'predictive' function)
//  * most of these functions used to use complete sorting, which is not too optimal,
//    especially for large data...
//  * later on I have rewritten it using partial sorting routines from R

// 23/03/2004: 'cumsumQuantile1'
//             'cumsumQuantile2'
//             'meanQuantile1'
//             'meanQuantile2'
// 07/02/2005: rewritten using partial sorting
//             'resetAverage'
//             'cumsum2average'
//             'value2quantile'
//
#include "quantile.h"

using namespace std;

// =================================================
// cumsumQuantile1: only update the cumulative sum
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
  int obs;
  for (obs = 0; obs < nobs; obs++){
    quant[obs][nquant] += newval[obs][iter];
  }
  return;
}


// =================================================
// cumsumQuantile2: only update the cumulative sum
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
  int obs, j;
  for (obs = 0; obs < nobs; obs++){
    for (j = 0; j < ngridM[obs]; j++){
      quant[obs][j][nquant] += newval[obs][j][iter];
    }
  }
  return;
}


// ====================================================
// meanQuantile1: compute mean and neeed quantiles
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

    for (i = 0; i < nquant; i++){
      rPsort(sample[obs], sampleSize, indquant1[i]);
      if (indquant2[i] != indquant1[i]){
        rPsort(sample[obs], sampleSize, indquant2[i]);
        quant[obs][i] = quantileA[i]*sample[obs][indquant1[i]] + (1 - quantileA[i])*sample[obs][indquant2[i]];
      }
      else{
        quant[obs][i] = sample[obs][indquant1[i]];
      }
    }

    /** Old code using complete sorting:  **/
    //    sort(sample[obs], sample[obs] + sampleSize);
    //    for (i = 0; i < nquant; i++) 
    //      quant[obs][i] = quantileA[i]*sample[obs][indquant1[i]] + (1 - quantileA[i])*sample[obs][indquant2[i]];

  }
 
  return;
}


// ====================================================
// meanQuantile2: compute mean and needed quantiles
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
    Rprintf("\n observ. %d", obs);
    for (j = 0; j < ngridM[obs]; j++){    
      quant[obs][j][nquant] /= sampleSize;

      for (i = 0; i < nquant; i++){
        rPsort(sample[obs][j], sampleSize, indquant1[i]);
        if (indquant2[i] != indquant1[i]){
          rPsort(sample[obs][j], sampleSize, indquant2[i]);
          quant[obs][j][i] = quantileA[i]*sample[obs][j][indquant1[i]] + (1 - quantileA[i])*sample[obs][j][indquant2[i]];
        }
        else{
          quant[obs][j][i] = sample[obs][j][indquant1[i]];
        }
      }

      /** Old code using complete sorting:  **/
      //    sort(sample[obs][j], sample[obs][j] + sampleSize);
      //      for (i = 0; i < nquant; i++) 
      //        quant[obs][j][i] = quantileA[i]*sample[obs][j][indquant1[i]] + (1 - quantileA[i])*sample[obs][j][indquant2[i]];    
    }
    Rprintf("  Done.");   
  }
  Rprintf("\n");
 
  return;
}


// =================================================
// resetAverage
// =================================================
// * used in predictive_GS
//
// Reset values used to store averages
//
// average ........ [ngrid[0] + ... + ngrid[nobs-1]]
//     on OUTPUT: filled with 0
// nobs ........... number of observations
// ngrid .......... number of gridpoints per observation
// predict ........ if 0, nothing is performed by this function
//
void
resetAverage(double* average,  const int* nobs,  const int* ngrid,  const int* predict)
{
  if (!(*predict)) return;

  int i, ix;
  double* tmpdP = average;
  for (i = 0; i < *nobs; i++){
    for (ix = 0; ix < ngrid[i]; ix++){
      *tmpdP = 0.0;  
      tmpdP++;
    }
  }
  return;
}

// =================================================
// cumsum2average
// =================================================
// * used in predictive_GS
//
// Compute averages from cumulative sums
//
// average ........ [ngrid[0] + ... + ngrid[nobs-1]]
//     on INPUT:  cumulative sums
//     on OUTPUT: averages
// sampleSize ..... sample size
// nobs ........... number of observations
// ngrid .......... number of gridpoints per observation
// predict ........ if 0, nothing is performed by this function
//
void
cumsum2average(double* average,  const int* sampleSize,  const int* nobs,  const int* ngrid,  const int* predict)
{
  if (!(*predict)) return;

  int i, ix;
  double* tmpdP = average;
  for (i = 0; i < *nobs; i++){
    for (ix = 0; ix < ngrid[i]; ix++){
      *tmpdP /= (*sampleSize);  
      tmpdP++;
    }
  }
  return;
}


// =================================================
// values2quantile
// =================================================
// * used in predictive_GS
//
// Compute quantiles for predictive quantities
//
// value ......... [shift_pointer*(ngrid[0] + ... + ngrid[nobs-1])]
//                 * values from which the quantiles are to be computed
//                 * stored in a long array as in predictive_GS
// quantile ...... [nquant*(ngrid[0] + ... + ngrid[nobs-1])]
//                 * computed quantiles
//                 * stored in a long array as in predictive_GS
// probs ......... [nquant] probabilities of quantiles
// indquant1 ..... [nquant] first indeces of observations in the whole chain 
//                          corresponding to desired quantiles
// indquant2 ..... [nquant] second indeces of observations in the whole chain 
//                          corresponding to desired quantiles
// nquant
// sampleSize ....... usually M_now
// nobs
// ngrid
// predict
// shift_pointer .... usually M_now_max
//
void
value2quantile(double* value,          double* quantile,  
               const double* probs,    const int* indquant1,  const int* indquant2,  const int* nquant,
               const int* sampleSize,  const int* nobs,       const int* ngrid,      const int* predict,    const int* shift_pointer)
{
  if (!(*predict)) return;

  int i, iq, ix;
  double* valP = value;
  double* quantP = quantile;
  for (i = 0; i < *nobs; i++){
    Rprintf("\n observ. %d", i);
    for (ix = 0; ix < ngrid[i]; ix++){
      for (iq = 0; iq < *nquant; iq++){
        rPsort(valP, *sampleSize, indquant1[iq]);
        if (indquant2[iq] != indquant1[iq]){
          rPsort(valP, *sampleSize, indquant2[iq]);
          quantP[iq*ngrid[i]+ix] = probs[iq]*valP[indquant1[iq]] + (1 - probs[iq])*valP[indquant2[iq]];
        }
        else{
          quantP[iq*ngrid[i]+ix] = valP[indquant1[iq]];
        }
      }
      valP += (*shift_pointer);
    }
    quantP += ngrid[i]*(*nquant);
    Rprintf("  Done.");
  }
  Rprintf("\n");

  return;
}

