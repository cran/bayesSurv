// Stuff for working with multivariate distributions
// 
// 20/01/2004: multivariate normal
// 09/02/2004: wishart

#include "mvtdist.h"

extern "C"{

using namespace std;

// ===================================================================
// Compute L*A*L' where
//   A is a symmetric matrix
//   L is a lower triangular matrix
// ===================================================================
//
// LAL ......... result (its lower triangle)
//               -> array of length 0.5*nA*(nA+1)
// L ........... lower triangle of L
//               -> array of length 0.5*nA*(nA+1)
// A ........... lower triangle of A
//               -> array of length 0.5*nA*(nA+1)
// nA .......... dimension of both L and A
// diagI ...... indeces of diagonal elements of the covariance matrix in the array A
//               -> array of length nA,
//               diagI[j] = (j * (2*(*nA) - j + 1)) / 2
//
void
LxMxtL(double* LAL,     const double* L,   const double* A,
       const int* nA,   const int* diagI)
{
  int j, i, k, minij;

  // First, compute L*A
  double** LA = new double*[*nA];  
  for (j = 0; j < *nA; j++){
    LA[j] = new double[*nA];               // jth column of L*A
    for (i = 0; i < *nA; i++){
      minij = (i <= j ? i  : j);
      LA[j][i] = L[i] * A[j];
      for (k = 1; k <= minij; k++){
        LA[j][i] += L[diagI[k] + i - k] * A[diagI[k] + j - k];
      }
      for (k = minij + 1; k <= i; k++){
        LA[j][i] += L[diagI[k] + i - k] * A[diagI[j] + k - j];
      }
    }
  }
  
  // Second, compute L*A*L'
  for (j = 0; j < *nA; j++){
    for (i = j; i < *nA; i++){
      LAL[diagI[j] + i - j] = LA[0][i] * L[j];
      for (k = 1; k <= j; k++){
        LAL[diagI[j] + i - j] += LA[k][i] * L[diagI[k] + j - k];        
      }
    }
  }

  for (j = 0; j < *nA; j++) delete [] LA[j];
  delete [] LA;   
  return;

}    // end of the function LxMxtL


// ===================================================================
// Compute L'*A*L where
//   A is a symmetric matrix
//   L is a lower triangular matrix
// ===================================================================
//
// LAL ......... result (its lower triangle)
//               -> array of length 0.5*nA*(nA+1)
// L ........... lower triangle of L
//               -> array of length 0.5*nA*(nA+1)
// A ........... lower triangle of A
//               -> array of length 0.5*nA*(nA+1)
// nA .......... dimension of both L and A
// diagI ...... indeces of diagonal elements of the covariance matrix in the array A
//               -> array of length nA,
//               diagI[j] = (j * (2*(*nA) - j + 1)) / 2
//
void
tLxMxL(double* LAL,     const double* L,   const double* A,
       const int* nA,   const int* diagI)
{
  int j, i, k;

  // First, compute L'*A
  double** LA = new double*[*nA];  
  for (j = 0; j < *nA; j++){
    LA[j] = new double[*nA];              // jth column of L'*A
    for (i = 0; i <= j; i++){
      LA[j][i] = L[diagI[i]] * A[diagI[i] + j - i];
      for (k = i + 1; k <= j; k++){
        LA[j][i] += L[diagI[i] + k - i] * A[diagI[k] + j - k];
      }
      for (k = j + 1; k < *nA; k++){
        LA[j][i] += L[diagI[i] + k - i] * A[diagI[j] + k - j];
      }
    }

    for (i = j + 1; i < *nA; i++){
      LA[j][i] = L[diagI[i]] * A[diagI[j] + i - j];
      for (k = i + 1; k < *nA; k++){
        LA[j][i] += L[diagI[i] + k - i] * A[diagI[j] + k - j];
      }
    }
  }
  
  // Second, compute L'*A*L
  for (j = 0; j < *nA; j++){
    for (i = j; i < *nA; i++){
      LAL[diagI[j] + i - j] = LA[j][i] * L[diagI[j]];
      for (k = j + 1; k < *nA; k++){
        LAL[diagI[j] + i - j] += LA[k][i] * L[diagI[j] + k - j];        
      }
    }
  }

  for (j = 0; j < *nA; j++) delete [] LA[j];
  delete [] LA;   
  return;

}    // end of the function tLxMxL


// ====================================================================
// Compute A*a where A is a SYMMETRIC matrix stored as
//   a lower triangle in column major order in a one-dimensional array
//   and a a column vector 
// ====================================================================
// 
// Ma ......... result (array of length nA)
// a .......... vector of length at least nA
// A .......... a lower triangle of a symmetric matrix stored in a column major order
//              (array of length 0.5*nA*(nA+1))
// inda ....... indeces of a that correspond to the columns and rows of matrix A
//              -> array of length nA, containing a subset of elements 0, 1, ..., na-1
//              * this array is ignored if na == nA
// na ......... length of a
// nA ......... dimension of A
// diagI ...... indeces of diagonal elements of the covariance matrix in the array A
//               -> array of length nA,
//               diagI[j] = (j * (2*(*nA) - j + 1)) / 2
//
void
Mxa(double* Ma,     const double* a,  const double* A,  const int* inda,  
    const int* na,  const int* nA,    const int* diagI)
{
  int i, j;

  if (*na == *nA){
    for (j = 0; j < *nA; j++){
      Ma[j] = A[diagI[j]] * a[j];
      for (i = j + 1; i < *nA; i++) Ma[j] += A[diagI[j] + i - j] * a[i];
      for (i = 0; i < j; i++)       Ma[j] += A[diagI[i] + j - i] * a[i];
    }
  }
  else{
    for (j = 0; j < *nA; j++){
      Ma[j] = A[diagI[j]] * a[inda[j]]; 
      for (i = j + 1; i < *nA; i++) Ma[j] += A[diagI[j] + i - j] * a[inda[i]]; 
      for (i = 0; i < j; i++)       Ma[j] += A[diagI[i] + j - i] * a[inda[i]]; 
    }
  }

  return;
}


// ====================================================================
// Compute A*a where A is a submatrix of a SYMMETRIC matrix stored as
//   a lower triangle in column major order in a one-dimensional array
//   and a a column vector 
// (used in 'GIBBSproposalMeanRandom')
// ====================================================================
// 
// Ma ......... result (array of length na)
// a .......... vector of length na
// A .......... a lower triangle of a symmetric matrix stored in a column major order
//              (array of length 0.5*nA*(nA+1))
// indA ....... columns of A that correspond to the vector a
//              -> array of length na, containing a subset of elements 0, 1, ..., nA-1
//              * this array is ignored if na == nA
// na ......... length of a
// nA ......... dimension of A, nA >= na
// diagI ...... indeces of diagonal elements of the covariance matrix in the array A
//               -> array of length nA,
//               diagI[j] = (j * (2*(*nA) - j + 1)) / 2
//
void
Mxa2(double* Ma,     const double* a,  const double* A,  const int* indA,  
     const int* na,  const int* nA,    const int* diagI)
{
  int i, j;

  if (*na == *nA){
    for (j = 0; j < *nA; j++){
      Ma[j] = A[diagI[j]] * a[j];
      for (i = j + 1; i < *nA; i++) Ma[j] += A[diagI[j] + i - j] * a[i];
      for (i = 0; i < j; i++)       Ma[j] += A[diagI[i] + j - i] * a[i];
    }
  }
  else{
    for (j = 0; j < *na; j++){
      Ma[j] = A[diagI[indA[j]]] * a[j]; 
      for (i = 0; i < *na; i++){
        if (indA[i] > indA[j] ) Ma[j] += A[diagI[indA[j]] + indA[i] - indA[j]] * a[i];
        if (indA[i] < indA[j] ) Ma[j] += A[diagI[indA[i]] + indA[j] - indA[i]] * a[i];        
      }
    }
  }

  return;
}


// ====================================================================
// Compute W*a where W is an oblong taken from a SYMMETRIC matrix A
//  and a is a vector
// * W does not contain any diagonal elements of the matrix A
// (used in 'GIBBSproposalMeanRandom')
// ====================================================================
// 
// Wa[nrow] ........... result
// a[na] .............. vector a
// A[0.5*nA*(nA+1)] ... matrix A stored as a lower triangle
// indr[nrow] ......... indeces of rows of A that define W
// indc[na] ........... indeces of columns of A that define W
// na ................. length of a
// nA ................. dimension of A
// nrow ............... length of the result, number of rows taken from A to get W
// diagI[nA] .......... indeces of diagonal elements of A in an array A
//
void
Wxa(double* Wa,   const double* a,   const double* A,  const int* indr,  const int* indc,
    const int* na,  const int* nA,  const int* nrow,   const int* diagI)
{
  int i, j;

  for (i = 0; i < *nrow; i++){
    Wa[i] = 0.0;
    for (j = 0; j < *na; j++){
      if (indr[i] >= indc[j]) Wa[i] += A[diagI[indc[j]] + indr[i] - indc[j]] * a[j];
      else                    Wa[i] += A[diagI[indr[i]] + indc[j] - indr[i]] * a[j];
    }
  }
  
  return;
}



// ====================================================================
// Compute a'*A*a where A is a SYMMETRIC matrix stored as
//   a lower triangle in column major order in a one-dimensional array
//   and a a column vector 
// ====================================================================
// 
// aMa ........ result (scalar)
// a .......... vector of length at least nA
// A .......... a lower triangle of a symmetric matrix stored in a column major order
//              (array of length 0.5*nA*(nA+1))
// inda ....... indeces of a that correspond to the columns and rows of matrix A
//              -> array of length nA, containing a subset of elements 0, 1, ..., na-1
//              * this array is ignored if na == nA
// na ......... length of a
// nA ......... dimension of a
// diagI ...... indeces of diagonal elements of the covariance matrix in the array A
//               -> array of length nA,
//               diagI[j] = (j * (2*(*nA) - j + 1)) / 2
//
void
axMxa(double* aMa,    const double* a,  const double* A,  const int* inda,  
      const int* na,  const int* nA,    const int* diagI)
{
  int i, j;

  *aMa = 0.0;

  if (*na == *nA)
    for (i = 0; i < *nA; i++){
      (*aMa) += a[i] * A[diagI[i]] * a[i];
      for (j = i + 1; j < *nA; j++)
        (*aMa) += 2 * a[i] * A[diagI[i] + j - i] * a[j];
    }
  else  
    for (i = 0; i < *nA; i++){
      (*aMa) += a[inda[i]] * A[diagI[i]] * a[inda[i]];
      for (j = i + 1; j < *nA; j++)
        (*aMa) += 2 * a[inda[i]] * A[diagI[i] + j - i] * a[inda[j]];
    }

  return;
}   // end of the function axMxa


// ================================================================
// Proportion of the density of a multivariate normal distribution
//
// * this function evaluates only exp(-0.5*(x-mu)'var^{-1}(x-mu))
//   or the logarithm of this expression
// ================================================================
//
// OUTPUT PARAMETERS:
// 
// dens ........... array with evaluated density,               of length  nP
//
// INPUT PARAMETERS:
//
// x ........... array of values at each the density should be evaluated      
//               -> array of length nx * nP
//               (first nx values is the first vector at which the mvt normal density is to be evaluated, etc.)
// mean ........ mean of the normal distribution,                             
//               -> array of length nmean, arbitrary 
// vari ........ inversion of the covariance matrix of the mvt normal distribution
//               (its lower triangle stored in column major order), in the form as returned by the functions 'cholesky' and 'chinv',
//               array of length 0.5*nxrepl*(nxrepl+1)
// indx ........ indeces of 'mean' that corresponds to the mvt norm distribution we want to evaluate
//               (for the case that the mean of the desired distribution is just a subvector of the array 'mean')
//               * this array is ignored if nx = nmean
//               * only entries from 0 to nmean - 1 should appear in this array
//               * repeated entries are allowed (e.g. for the case that all components of the mvt distribution have the same mean)
//               * if nx != nmean, 
//                 -> this array should be of length nx
// indxrepl .... indeces of x that are to be considered for the evaluation
//                 -> array of length nxrepl containing subset of values 0, 1, ... nx - 1
// nx .......... length of each x
// nmean ....... length of the array 'mean'
// nxrepl ...... dimension of the multivariate normal distribution
// nP .......... number of the vectors that should be sampled
// logP ........ 0/1, do I want a logarithm of the value or not
//
void
dmvtnorm(double* dens,      const double* x,      const double* mean,  const double* vari,  
         const int* indx,   const int* indxrepl, 
         const int* nx,     const int* nmean,     const int* nxrepl,    
         const int* nP,     const int* diagI,     const int* logP)
{
  int sv, i;

  double* xminmu = new double[*nxrepl];
  for (sv = 0; sv < *nP; sv++){
    if ((*nx) == (*nmean)) for (i = 0; i < *nxrepl; i++) xminmu[i] = x[(*nx)*sv + indxrepl[i]] - mean[indxrepl[i]];
    else                   for (i = 0; i < *nxrepl; i++) xminmu[i] = x[(*nx)*sv + indxrepl[i]] - mean[indx[indxrepl[i]]];
    axMxa(dens + sv, xminmu, vari, indxrepl, nxrepl, nxrepl, diagI);
    dens[sv] *= (-0.5);    
    if (!(*logP)) dens[sv] = exp(dens[sv]);
  }

  delete [] xminmu;
  return;  
}   // end of the function dmvtnorm


// =============================================================
// Sample from a multivariate normal distribution,
//  using a Cholesky decomposition of the covariance matrix
// =============================================================
//
// OUTPUT PARAMETERS:
// 
// x ........... array with sampled values,                      of length nx * nP
//
// INPUT PARAMETERS:
//
// mean ........ mean of the normal distribution,                 
//               -> array of length nmean, arbitrary
// L ........... Cholesky decomposition of the covariance matrix of the normal distribution
//               (its lower triangle stored in column major order), in the form as returned by the function 'cholesky',
//               -> array of length 0.5*nxrepl*(nxrepl+1)
// indx ........ indeces of 'mean' that corresponds to x
//               (for the case that the mean of desired x is just a subvector of the array 'mean')
//               * this array is ignored if nx = nmean
//                 OR                    if nxrepl = nmean
//               * only entries from 0 to nmean - 1 should appear in this array
//               * repeated entries are allowed (e.g. for the case that all components of x have the same mean)
//               * if nx != nmean AND nxrepl != nmean, 
//                  -> array should be of length nx
// indxrepl .... indeces of the components of x vector that are to be sampled, the rest is left unchanged,
//               -> array of length nxrepl 
//               * if additionally nx != nmean AND nxrepl = nmean 
//                                 then Ex[indxrepl[0]] = mean[0], Ex[indxrepl[1]] = mean[1] etc.
// nx .......... length of an array x divided by nP 
// nmean ....... length of the array 'mean'
// nxrepl ...... dimension of the multivariate normal distribution
// nP .......... number of the vectors that should be sampled
// diagI ....... indeces of diagonal elements of the covariance matrix in the array varChol
//               -> array of length nxrepl,
//               diagI[j] = (j * (2*(*nxrepl) - j + 1)) / 2
// callFromR ... 0/1 indicating whether the function is called from R
//
void
rmvtnorm(double* x,        const double* mean,   const double* L,  
         const int* indx,  const int* indxrepl,
         const int* nx,    const int* nmean,     const int* nxrepl,      
         const int* nP,    const int* diagI,     const int* callFromR)
{
  int i, j, sv;

  if (*callFromR) GetRNGstate();

  double u;

  for (sv = 0; sv < *nP; sv++){

    if ((*nx) == (*nmean)) 
      for (i = 0; i < *nxrepl; i++) x[(*nx)*sv + indxrepl[i]] = mean[indxrepl[i]];
    else
      if ((*nxrepl) == (*nmean))
        for (i = 0; i < *nxrepl; i++) x[(*nx)*sv + indxrepl[i]] = mean[i];
      else                   
        for (i = 0; i < *nxrepl; i++) x[(*nx)*sv + indxrepl[i]] = mean[indx[indxrepl[i]]];

    for (j = 0; j < *nxrepl; j++){
      u = rnorm(0, 1);
      for (i = j; i < *nxrepl; i++){
         x[(*nx)*sv + indxrepl[i]] += u * L[diagI[j] + i - j];         
      } 
    }        
  }

  if (*callFromR) PutRNGstate();

  return;
}   // end of the function rmvtnorm


// =============================================================
// Sample from a multivariate normal distribution,
//  using an inverse of the Cholesky decomposition 
//  of the inverse of the covariance matrix
// =============================================================
//
// OUTPUT PARAMETERS:
// 
// x ........... array with sampled values,                      of length nx * nP
//
// INPUT PARAMETERS:
//
// mean ........ mean of the normal distribution,                 
//               -> array of length nmean, arbitrary
// iLi ......... inverse of the Cholesky decomposition of the inverse 
//               of the covariance matrix of the normal distribution
//               (its lower triangle stored in column major order), in the form as returned by the function 'cholesky',
//               -> array of length 0.5*nxrepl*(nxrepl+1)
// indx ........ indeces of 'mean' that corresponds to x
//               (for the case that the mean of desired x is just a subvector of the array 'mean')
//               * this array is ignored if nx = nmean
//                 OR                    if nxrepl = nmean
//               * only entries from 0 to nmean - 1 should appear in this array
//               * repeated entries are allowed (e.g. for the case that all components of x have the same mean)
//               * if nx != nmean AND nxrepl != nmean, 
//                  -> array should be of length nx
// indxrepl .... indeces of the components of x vector that are to be sampled, the rest is left unchanged,
//               -> array of length nxrepl 
//               * if additionally nx != nmean AND nxrepl = nmean 
//                                 then Ex[indxrepl[0]] = mean[0], Ex[indxrepl[1]] = mean[1] etc.
// nx .......... length of an array x divided by nP 
// nmean ....... length of the array 'mean'
// nxrepl ...... dimension of the multivariate normal distribution
// nP .......... number of the vectors that should be sampled
// diagI ....... indeces of diagonal elements of the covariance matrix in the array varChol
//               -> array of length nxrepl,
//               diagI[j] = (j * (2*(*nxrepl) - j + 1)) / 2
// callFromR ... 0/1 indicating whether the function is called from R
//
void
rmvtnorm2(double* x,        const double* mean,   const double* iLi,  
          const int* indx,  const int* indxrepl,
          const int* nx,    const int* nmean,     const int* nxrepl,      
          const int* nP,    const int* diagI,     const int* callFromR)
{
  int i, j, sv;

  if (*callFromR) GetRNGstate();

  double u;

  for (sv = 0; sv < *nP; sv++){

    if ((*nx) == (*nmean)) 
      for (i = 0; i < *nxrepl; i++) x[(*nx)*sv + indxrepl[i]] = mean[indxrepl[i]];
    else
      if ((*nxrepl) == (*nmean))
        for (i = 0; i < *nxrepl; i++) x[(*nx)*sv + indxrepl[i]] = mean[i];
      else                   
        for (i = 0; i < *nxrepl; i++) x[(*nx)*sv + indxrepl[i]] = mean[indx[indxrepl[i]]];

    for (i = *nxrepl - 1; i >= 0; i--){
      u = rnorm(0, 1);
      for (j = i; j >= 0; j--){
        x[(*nx)*sv + indxrepl[j]] += u * iLi[diagI[j] + i - j];
      }
    }
  }

  if (*callFromR) PutRNGstate();

  return;
}   // end of the function rmvtnorm2


// =============================================================
// Sample from a multivariate uniform distribution
//   with independent components
// =============================================================
//
// OUTPUT PARAMETERS:
// 
// x ........... array with sampled values,                      of length nx * nP
//
// INPUT PARAMETERS:
//
// mean ........ means of each uniform distribution,                 array of length nmean
// halfRange ... half ranges of each uniform distribution            array of length nmean
// indx ........ indeces of 'mean' and 'halfRange' that corresponds to x
//               (for the case that the mean of desired x is just a subvector of the array 'mean')
//               * this array is ignored if nx = nmean
//               * only entries from 0 to nmean - 1 should appear in this array
//               * repeated entries are allowed (e.g. for the case that all components of x have the same mean)
//               * if nx != nmean, this array should be of length nx
// indxrepl .... indeces of the components of x vector that are to be sampled, the rest is left unchanged,
//               -> array of length nxrepl, with elements from 0, 1, ..., nx - 1
// nx .......... length of an array x divided by nP 
// nmean ....... length of the array 'mean' and 'halfRange'
// nxrepl ...... dimension of the multivariate uniform distribution
// nP .......... number of the vectors that should be sampled
// callFromR ... 0/1 indicating whether the function is called from R
//
void
rmvtiunif(double* x,        const double* mean,   const double* halfRange,  
          const int* indx,  const int* indxrepl,  
          const int* nx,    const int* nmean,     const int* nxrepl,
          const int* nP,    const int* callFromR)
{
  int i, sv;

  if (*callFromR) GetRNGstate();

  double u;

  if ((*nx) == (*nmean))
    for (sv = 0; sv < *nP; sv++)
      for (i = 0; i < *nxrepl; i++){
        u = rnorm(0, 1);
        x[(*nx)*sv + indxrepl[i]] = mean[indxrepl[i]] - halfRange[indxrepl[i]] + 2 * halfRange[indxrepl[i]] * u; 
      }
  else
    for (sv = 0; sv < *nP; sv++)
      for (i = 0; i < *nxrepl; i++){
        u = rnorm(0, 1);
        x[(*nx)*sv + indxrepl[i]] = mean[indx[indxrepl[i]]] - halfRange[indx[indxrepl[i]]] + 2 * halfRange[indx[indxrepl[i]]] * u; 
      }

  if (*callFromR) PutRNGstate();

  return;
}    // end of the function rmvtiunif


// =============================================================
// Sample from a Wishart distribution (nu, S), version 1
//   * decomposition of a scale matrix is given
//   * parametrization as in Gelman et al. (2004), appendix A
// =============================================================
//
// OUTPUT PARAMETERS:
//
// w ..... lower triangles of the sampled matrices with the Wishart distribution
//
// INPUT PARAMETERS:
//
// p ........... dimension of the Wishart distribution
// nu .......... degrees of freedom of the Wishart distribution
//               * it does not have to be an integer
//               * it must only satisfy: nu > p - 1
// L ........... decomposition of the scale matrix S such that S = L*L',
//               * only its lower triangle stored in an array
//               * it is assumed that the upper triangle contains only zeros
// diagI ....... indeces of diagonal elements in arrays
// nP .......... number of the matrices that should be sampled
// callFromR ... 0/1 indicating whether the function is called from R
//
void
rwishart(double* w,
         const int* p,      const double* nu,   const double* L,
         const int* diagI,  const int* nP,      const int* callFromR)
{
  try{
    if (*nu <= *p - 1) returnR("C++ Error: Incorrect degrees of freedom for a Wishart distribution.", -1);
    if (*callFromR) GetRNGstate();

    int i, j, k, sv;
    double shape, scale;

    if (*p == 1){          // Wishart(nu, S) = gamma(nu/2, 1/2S) = gamma(shape = nu/2, scale = 2S)
      shape = *nu/2;
      scale = 2*(*L)*(*L);
      for (sv = 0; sv < *nP; sv++){    
        w[sv] = rgamma(shape, scale);
      }
    }
    else{                  // use a Bartlett's decomposition for sampling
      int lw = ((*p)*(*p+1))/2;                    // length of one w
      double* norm = new double[lw];
      double* v = new double[lw];
      double zeta;
      scale = 2.0;
      for (sv = 0; sv < *nP; sv++){

        // V matrix
        for (j = 0; j < *p; j++){    // loop over columns
          shape = 0.5*(*nu - j);                // Note: C++ index i <=> index i + 1 in a theory
          zeta = rgamma(shape, scale);
          v[diagI[j]] = zeta;
          norm[diagI[j]] = sqrt(zeta);
          for (i = j + 1; i < *p; i++){
            norm[diagI[j] + i - j] = rnorm(0, 1);
            v[diagI[j] + i - j] = norm[diagI[j] + i - j] * norm[diagI[j]];
          }
          for (i = j; i < *p; i++){
            for (k = j - 1; k >= 0; k--){
              v[diagI[j] + i - j] += norm[diagI[k] + i - k] * norm[diagI[k] + j - k];
            }
          }
        }

        // LVL'
        LxMxtL(w + sv*lw, L, v, p, diagI);        
      }
      delete [] norm;      
      delete [] v;
    }

    if (*callFromR) PutRNGstate();
    return;
  }
  catch(returnR){
    if (*callFromR){
      PutRNGstate();
      return;
    }
    throw;
  }  
}    // end of the function rwishart


// ================================================================
// Sample from a Wishart distribution (nu, S), version 2
//   * inverse decomposition of an inverse scale matrix is given
//   * parametrization as in Gelman et al. (2004), appendix A
// ================================================================
//
// OUTPUT PARAMETERS:
//
// w ..... lower triangles of the sampled matrices with the Wishart distribution
//
// INPUT PARAMETERS:
//
// p ........... dimension of the Wishart distribution
// nu .......... degrees of freedom of the Wishart distribution
//               * it does not have to be an integer
//               * it must only satisfy: nu > p - 1
// iLi ......... inverse decomposition of an inverse of the scale matrix S,
//               i.e. S^{-1} = L*L', iLi = L^{-1},
//                    S = L^{-1}'*L^{-1}
//               * only its lower triangle stored in an array
//               * it is assumed that the upper triangle contains only zeros
// diagI ....... indeces of diagonal elements in arrays
// nP .......... number of the matrices that should be sampled
// callFromR ... 0/1 indicating whether the function is called from R
//
void
rwishart2(double* w,
          const int* p,      const double* nu,   const double* iLi,
          const int* diagI,  const int* nP,      const int* callFromR)
{
  try{
    if (*nu <= *p - 1) returnR("C++ Error: Incorrect degrees of freedom for a Wishart distribution.", -1);
    if (*callFromR) GetRNGstate();

    int i, j, k, sv;
    double shape, scale;

    if (*p == 1){          // Wishart(nu, S) = gamma(nu/2, 1/2S) = gamma(shape = nu/2, scale = 2S)
      shape = *nu/2;
      scale = 2*(*iLi)*(*iLi);
      for (sv = 0; sv < *nP; sv++){    
        w[sv] = rgamma(shape, scale);
      }
    }
    else{                  // use a Bartlett's decomposition for sampling
      int lw = ((*p)*(*p+1))/2;                    // length of one w
      double* norm = new double[lw];
      double* v = new double[lw];
      double zeta;
      scale = 2.0;
      for (sv = 0; sv < *nP; sv++){

        // V matrix
        for (j = 0; j < *p; j++){    // loop over columns
          shape = 0.5*(*nu - j);                // Note: C++ index i <=> index i + 1 in a theory
          zeta = rgamma(shape, scale);
          v[diagI[j]] = zeta;
          norm[diagI[j]] = sqrt(zeta);
          for (i = j + 1; i < *p; i++){
            norm[diagI[j] + i - j] = rnorm(0, 1);
            v[diagI[j] + i - j] = norm[diagI[j] + i - j] * norm[diagI[j]];
          }
          for (i = j; i < *p; i++){
            for (k = j - 1; k >= 0; k--){
              v[diagI[j] + i - j] += norm[diagI[k] + i - k] * norm[diagI[k] + j - k];
            }
          }
        }

        // LVL'
        tLxMxL(w + sv*lw, iLi, v, p, diagI);        
      }
      delete [] norm;      
      delete [] v;
    }

    if (*callFromR) PutRNGstate();
    return;
  }
  catch(returnR){
    if (*callFromR){
      PutRNGstate();
      return;
    }
    throw;
  }  
}    // end of the function rwishart2


// ==============================================================================
// Sample from an inverse-Wishart distribution (nu, S^{-1}), where S^{-1} = L*L'
//   * parametrization as in Gelman et al. (2004), appendix A,
//   * i.e. S = scale of the inverse Wishart distribution
// ==============================================================================
//
// OUTPUT PARAMETERS:
//
// w ..... lower triangles of the sampled matrices with the inverse-Wishart distribution
//
// INPUT PARAMETERS:
//
// p ........... dimension of the Wishart distribution
// nu .......... degrees of freedom of the inverse-Wishart distribution
//               * it does not have to be an integer
//               * it must only satisfy: nu > p - 1
// L ........... decomposition of the inverse of the scale matrix S such that S^{-1} = L*L',
//               * only its lower triangle stored in an array
//               * it is assumed that the upper triangle contains only zeros
// diagI ....... indeces of diagonal elements in arrays
// nP .......... number of the matrices that should be sampled
// callFromR ... 0/1 indicating whether the function is called from R
//
void
rinvwishart(double* w,
            const int* p,      const double* nu,   const double* L,
            const int* diagI,  const int* nP,      const int* callFromR)
{
  const double tolerChol = 1e-7;
  int sv;

  // Sample from Wishart(nu, L*L')
  rwishart(w, p, nu, L, diagI, nP, callFromR);

  // Compute inversions
  if (*p == 1){
    for (sv = 0; sv < *nP; sv++){
      if (w[sv] <= invFLT_MAX) w[sv] = FLT_MAX;
      else                     w[sv] = 1/w[sv];
    }
  }
  else{
    int lw = ((*p)*(*p+1))/2;                    // length of one w 
    int rank;
    for (sv = 0; sv < *nP; sv++){
      cholesky(w + sv*lw, &rank, p, diagI, &tolerChol);
      chinv(w + sv*lw, p, diagI, &ZERO_INT);
    }
  }

  return;
}    // end of function rinvwishart



}    // end of extern "C"
