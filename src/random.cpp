// Additional random generators and related functions

// 27/11/2003: start working on it
// 02/08/2004: left-truncated gamma distribution, random generator added

#include "bayessurvreg.h"

extern "C"{

using namespace std;

// ********** rltruncGamma ************
// Sample from left-truncated gamma distribution
//
// PARAMETERS:
//   shape .......... shape of the gamma distribution
//   rate ........... rate (= 1/scale) of the gamma distribution
//   minx ........... a value > 0 indicating where the gamma distribution is to be truncated
//   n .............. number of random numbers that should be returned in the vector x
//   callFromR ...... logical true = if this function is called directly from R (random seed must be got and put)
//   
void
rltruncGamma(double* x,  
             const double* shape,  const double* rate,   const double* minx,
             const int* n,         const int* callFromR)
{
  int i;
  double u;

  if (*callFromR) GetRNGstate();

  double scale = 1/(*rate);
  double Flower = pgamma(*minx, *shape, scale, 1, 0);
  if (Flower >= 1 - NORM_ZERO)   // truncation time irrealistic large => sampled values = minx
    for (i = 0; i < *n; i++) x[i] = *minx;
  else
    if (Flower <= NORM_ZERO)     // truncation time = 0, sample from gamma distribution
      for (i = 0; i < *n; i++) x[i] = rgamma(*shape, scale);
    else
      for (i = 0; i < *n; i++){
        u = runif(0, 1) * (1 - Flower) + Flower;
        x[i] = qgamma(u, *shape, scale, 1, 0);
      }

  if (*callFromR) PutRNGstate();
  return;
}

// ********** discreteUniformSampler **********
// Sample from a discrete distribution on 0, 1, ..., *kP - 1
//   where P(Y = j) = 1/*kP for all j = 0, 1, ..., *kP - 1
//
// ASSUMPTIONS:
//   * kP > 0
//
// I do not check this assumption in the C++ code!!!
//
// PARAMETERS:
//
// kP .................... length of the uniform distribution
// n ..................... number of random variates to be sampled 
// callFromR ............. logical true = if this function is called directly from R (random seed must be got and put)
//
// RETURN:
//
// sampledj .............. array of sampled values
//
void
discreteUniformSampler(int* sampledj,
                       const int* kP,
                       const int* nP,
                       const int* callFromR) 
{
  int i;
  double u;

  if (*kP <= 1){
    for (i = 0; i < *nP; i++) sampledj[i] = 0;
    return;
  }

  if (*callFromR) GetRNGstate();

  for (i = 0; i < *nP; i++){
     u = runif(0, 1);
     sampledj[i] = findUniformIndex(u, 0, *kP - 1, *kP);
  }

  if (*callFromR) PutRNGstate();
  return;
}


// ********** discreteSampler **********
// Sample from a discrete distribution on 0, 1, ..., *kP - 1
//   when either cumulative proportions are given or proportions are given,
//   i.e. if (cumul){
//          P(Y = 0) \propto propA[0]
//          P(Y = j) \propto propA[j] - propA[j-1], j = 1, ..., *kP - 1
//        }
//        else{
//          P(Y = j) \propto  propA[j], j = 0, ..., *kp - 1
//        }
//
// ASSUMPTIONS:
//   * if (cumul)
//       (i)  \exist j propA[j] > 0 
//       (ii) propA[0] <= ... <= propA[*kP - 1]
//   * if (!cumul)
//       (ii) \exist j propA[j] > 0  
//       (i)  propA[j] => 0 for all j  
//
// I do not check these assumptions in the C++ code!!!
//
// PARAMETERS:
//
// propA ................. array of either proportions or cumulative proportions
// kP .................... length of the array propA
// n ..................... number of random variates to be sampled 
// cumul ................. logical, true if the array propA gives cumulative proportions
// callFromR ............. logical true = if this function is called directly from R (random seed must be got and put)
//
// RETURN:
//
// sampledj .............. array of sampled values
//
void
discreteSampler(int* sampledj,
                double* propA,
                const int* kP,
                const int* nP,
                const int* cumul,
                const int* callFromR)
{
  int i, j;
  double u;

  if (*kP <= 1){
    for (i = 0; i < *nP; i++) sampledj[i] = 0;
    return;
  }

  if (*callFromR) GetRNGstate();

  // Compute cumulative proportions if necessary
  if (!(*cumul)){
    for (j = 1; j < *kP; j++)
      propA[j] += propA[j-1];
  }

  // Check whether all cumulative probabilities are distinct
  int* origInd = new int[*kP];                // correspondence between original and non-zero indeces
  double* nonZeropropA = new double[*kP];     // cumulative probabilities for non-zero components only
  int kNonZero = 1;
  j = 0;
  while (propA[j] <= ZERO) j++;
  origInd[0] = j;
  nonZeropropA[0] = propA[j];
  for (int m = j + 1; m < *kP; m++){
    if (propA[m] - propA[m - 1] <= ZERO) continue;
    origInd[kNonZero] = m;
    nonZeropropA[kNonZero] = propA[m]; 
    kNonZero++;
  }

  int sampledInd;
  if (kNonZero == 1) for (i = 0; i < *nP; i++) sampledj[i] = origInd[0];
  else               for (i = 0; i < *nP; i++){
                        u = runif(0, nonZeropropA[kNonZero - 1]);
                        sampledInd = findIndex(u, 0, kNonZero - 1, nonZeropropA);
                        sampledj[i] = origInd[sampledInd];
                     }

  if (*callFromR) PutRNGstate();
  delete [] origInd;
  delete [] nonZeropropA;
  return;  

}   // end of function discreteSampler


// ********** findUniformIndex **********
// Find the index of the smallest element of ValuesA = (1/k, 2/k, ..., k/k)
//  which is higher than u 
// The search is restricted to ValuesA[startInd], ..., ValuesA[endInd]
//
// This function is to be used by 'discreteUniformSampler' routine.
// That's why it relies on the following assumptions:
//   * startInd < endInd
//   * 0 <= u <= 1
//
// I do not check these assumptions since it would be wasting of CP time.
//
int 
findUniformIndex(const double u,
                 const int startInd,
                 const int endInd,
                 const int k)
{
  if (startInd == endInd - 1){
    if (u <= double(startInd+1)/double(k)) return startInd;
    else                                   return endInd;
  } 
  else{
    int midInd = int(ceil(0.5 * (startInd + endInd)));
    int Index;
    if (u <= double(midInd+1)/double(k)) Index = findUniformIndex(u, startInd, midInd, k);
    else                                 Index = findUniformIndex(u, midInd, endInd, k); 
    return Index;
  }
}   // end of function findIndex



// ********** findIndex **********
// Find the index of the smallest element of ValuesA
//  which is higher than u 
// The search is restricted to ValuesA[startInd], ..., ValuesA[endInd]
//
// This function is to be used by 'discreteSampler' routine.
// That's why it relies on the following assumptions:
//   * startInd < endInd
//   * ValuesA[0] > 0
//   * ValuesA[startInd] < ... < ValuesA[endInd]
//   * 0 <= u <= valuesA[endInd]
//
// I do not check these assumptions since it would be wasting of CP time.
//
int 
findIndex(const double u,
          const int startInd,
          const int endInd,
          const double* ValuesA)
{
  if (startInd == endInd - 1){
    if (u <= ValuesA[startInd]) return startInd;
    else                        return endInd;
  } 
  else{
    int midInd = int(ceil(0.5 * (startInd + endInd)));
    int Index;
    if (u <= ValuesA[midInd]) Index = findIndex(u, startInd, midInd, ValuesA);
    else                      Index = findIndex(u, midInd, endInd, ValuesA); 
    return Index;
  }
}   // end of function findIndex

}   // end of extern "C"
