// Routines to compute a Cholesky decomposition and related functions of a real symmetric matrix
//   given as a lower triangle stored in an array in column major order.

// 19/01/2004: start working on it

// These programs are based on functions 'cholesky2', 'chsolve2' and 'chinv2' of Terry Therneau
// included in the R library 'survival'.

#include "bayessurvreg.h"

extern "C"{

using namespace std;

// ================================================================================
// Function to compute a Cholesky decomposition of a real symmetric matrix
//   stored as its lower triangle in an array in column major order.
// ================================================================================
//
// Cholesky decompostion on a matrix: C = LL'
//   where L is lower triangular matrix with positive entries on a diagonal for C = positive semidefinite
//
// If C is not positive semidefinite, we return C = LDL'
//   where L is a lower triangular with 1 on a diagonal
//   and D is diagonal
//
// INPUT PARAMETERS:
//    C ....... input matrix (it must be of length 0.5 * nC * (nC+1))
//              !!! C is rewritten by its decomposition
//    nC ...... number of rows (and columns) of matrix C
//    toler ... threshhold value for detecting singularity
//    diagI ... indeces of diagonal elements of C in the array
//              (diagI[j] = (j * (2*(*nC) - j + 1)) / 2)
//              -> array of length nC
//
//  OUTPUT PARAMETERS:
//    C ....... Cholesky decomposition of C.
//              * the factorization is returned in the lower triangle of the original matrix,
//              * D occupies the diagonal 
//    rankC ... the rank of the matrix (non-negative definite), 
//              or -rank if not SPD or NND
//      
void
cholesky(double* C,  int* rankC,  const int* nC,  const int* diagI, const double* toler)
{
    double temp;
    int  i, j, k;
    double eps, pivot;
    int nonneg;

    if (*nC == 1){
      *rankC = 1*(C[0] > *toler) + (-1)*(C[0] < -(*toler));
      if (*rankC) C[0] = sqrt(C[0]);
      return;
    }

    nonneg= 1;                               // At the beginning, I believe that C is positive definite
    eps = 0;
    for (j = 0; j < *nC; j++) {
      if (fabs(C[diagI[j]]) > eps)  eps = fabs(C[diagI[j]]);              // eps = max(C(i,i))
    }
    eps *= (*toler);                                                      // eps = toler * max(C(i,i))

    *rankC = 0;
    for (j = 0; j < *nC; j++) {
	pivot = C[diagI[j]];
        if (pivot < -8*eps) nonneg= -1;
	if (fabs(pivot) < eps) {
	   C[diagI[j]] = 0;
	}
	else  {
          (*rankC)++;
	  for (i = (j+1); i < *nC; i++) {
	    temp = C[diagI[j] + i - j]/pivot;
	    C[diagI[j] + i - j] = temp;
	    C[diagI[i]] -= temp*temp*pivot;
	    for (k = (i+1); k < *nC; k++) C[diagI[i] + k - i] -= temp*C[diagI[j] + k - j];
	  }
	}
    }
    (*rankC) *= nonneg;

    // Compute the square root of the diagonal
    if (nonneg == 1){
      for (j = 0; j < *nC; j++){
        C[diagI[j]] = sqrt(C[diagI[j]]);
      }        
    }

    // Multiply each column under the diagonal by the diagonal element
    for (j = 0; j < *nC - 1; j++){
      for (i = j + 1; i < *nC; i++ ){
        C[diagI[j] + i - j] *= C[diagI[j]];
      }
    }

    return;
}


// ================================================================================
// Function to compute a Cholesky decomposition of a real symmetric matrix
//   stored in an array in column major order.
// ================================================================================
// INPUT PARAMETERS:
//    C ....... input matrix (it must be an array of length nC * nC)
//              the upper triangle is taken to represent the matrix C,
//              on input: the lower triangle is ignored
//              !!! lower triangle of C is rewritten by its decomposition
//    nC ...... number of rows (and columns) of matrix C
//    toler ... threshhold value for detecting singularity
//
//  OUTPUT PARAMETERS:
//    C ....... Cholesky decomposition of C.
//              * the factorization is returned in the lower triangle of the original matrix,
//              * D (if C is not positive semidefinite) occupies the diagonal 
//    rankC ... the rank of the matrix (non-negative definite), 
//              or -rank if not SPD or NND
//      
void
cholesky2(double* C,  int* rankC,  const int* nC,  const double* toler)
{
    double temp;
    int  i, j, k;
    double eps, pivot;
    int nonneg;

    if (*nC == 1){
      *rankC = 1*(C[0] > *toler) + (-1)*(C[0] < -(*toler));
      return;
    }

    int n = *nC;
    nonneg= 1;          // At the beginning, I believe that C is positive definite
    eps = 0;
    for (i = 0; i < n; i++) {
      if (fabs(C[i*n+i]) > eps)  eps = fabs(C[i*n+i]);                // eps = max(C(i,i))
      for (j = (i+1); j < n; j++)  C[i*n+j] = C[j*n+i];
    }
    eps *= (*toler);                                                  // eps = toler * max(C(i,i))

    *rankC = 0;
    for (i = 0; i < n; i++) {
	pivot = C[i*n+i];
        if (pivot < -8*eps) nonneg= -1;
	if (fabs(pivot) < eps) {
	   C[i*n+i] = 0;
	}
	else  {
	    (*rankC)++;
	    for (j = (i+1); j < n; j++) {
	      temp = C[i*n+j]/pivot;
	      C[i*n+j] = temp;
	      C[j*n+j] -= temp*temp*pivot;
	      for (k = (j+1); k < n; k++) C[j*n+k] -= temp*C[i*n+k];
	    }
	}
    }
    (*rankC) *= nonneg;

    // Compute the square root of the diagonal
    if (nonneg == 1){
      for (j = 0; j < *nC; j++){
        C[j*n+j] = sqrt(C[j*n+j]);
      }        
    }

    // Multiply each column under the diagonal by the diagonal element
    for (j = 0; j < *nC - 1; j++){
      for (i = j + 1; i < *nC; i++ ){
        C[j*n+i] *= C[j*n+j];
      }
    }


    return;
}


// ==================================================================================
// Function to compute an inversion of a real symmetric, positive definite matrix
//  from its Cholesky decomposition stored as a lower triangle in an array.
// ==================================================================================
// INPUT PARAMETERS:
//
//    C ............. input matrix (it must be of length 0.5 * nC * (nC+1))
//                    = Cholesky decomposition of the original matrix obtained by the function 'cholesky'
//                    !!! C is rewritten by its inversion
//                    !!! It is assumed that C is positive definite, it is not checked!!!
//                    !!! Check this using the output from 'cholesky' function.
//    nC ............ number of rows (and columns) of matrix C
//    diagI ... indeces of diagonal elements of C in the array
//              (diagI[j] = (j * (2*(*nC) - j + 1)) / 2)
//              -> array of length nC
//    onlyCholInv ... 1 if only inversion of the Cholesky decomposition (L matrix) is desired
//
//  OUTPUT PARAMETERS:
//    C ....... inversion of C (its lower triangle)
//
void 
chinv(double* C , const int* nC,  const int* diagI, const int* onlyCholInv)
{

  double temp;
  int i, j, k;

  // * divide each column of the Cholesky decomposition under the diagonal by the diagonal element
  for (j = 0; j < *nC - 1; j++){
    if (C[diagI[j]] != 0){
      for (i = j + 1; i < *nC; i++ ){
        C[diagI[j] + i - j] /= C[diagI[j]];
      }
    }
  }

  // * (almost) invert the cholesky in the lower triangle
  // * do not multiply the rows in the lower strict triangle by the inverse of the diagonal element on that row
  // * this corresponds to L^{-1} where L had ones on the diagonal,
  for (j = 0; j < *nC; j++){                   // loop over columns
    if (C[diagI[j]] > 0) {
      C[diagI[j]] = 1/(C[diagI[j]]);                            /*this line inverts a diagonal */
      for (i= (j+1); i < *nC; i++){
	C[diagI[j] + i - j] = -C[diagI[j] + i - j];
  	for (k = 0; k < j; k++){                                /*sweep operator */
          C[diagI[k] + i - k] += C[diagI[j] + i - j]*C[diagI[k] + j - k];
        }
      }
    }
  }

  if (*onlyCholInv){
    for (i = 0; i < *nC; i++)        // loop over rows (each of them must be multiplied by the corresponding diagonal element)
      if (C[diagI[i]] == 0) for (j = 0 ; j < i; j++) C[diagI[j] + i - j] = 0;              /* singular row */
      else                  for (j = 0 ; j < i; j++) C[diagI[j] + i - j] *= C[diagI[i]];
    return;
  }


  // * lower triangle now contains inverse of cholesky
  // * calculate (L^{-1})'DL^{-1} (inverse of cholesky decomp process) 
  //   where D is a diagonal matrix with 1/(diag(L)*diag(L))
  //   and L had ones on a diagonal
  // * this will be the inverse of the original matrix

  // Store L^{-1}D matrix in a new array
  double* LiD = new double[((*nC)*(*nC + 1))/2];
  for (j = 0; j < *nC; j++){
    LiD[diagI[j]] = C[diagI[j]] * C[diagI[j]];
    for (i = j + 1; i < *nC; i++)
      LiD[diagI[j] + i - j] = C[diagI[j] + i - j];
  }

  // Compute L^{-1}DL^{-1}
  for (j= 0; j < *nC; j++) {             // loop over columns
    if (C[diagI[j]] == 0){                         /* singular column */
      for (i = j; i < *nC; i++) C[diagI[j] + i - j] = 0;
    }
    else {
      C[diagI[j]] = LiD[diagI[j]];                                           // initialize the jth diagonal element
      if (*nC - 1 > j){
        i = *nC - 1;
        C[diagI[j] + i - j] = LiD[diagI[i]] * LiD[diagI[j] + i - j];         // last row of the jth column
        for (k = i - 1; k > j; k--)                                          // initialize rows j + 1, ..., n - 2 of the jth column
          C[diagI[j] + k - j] = C[diagI[j] + i - j] * LiD[diagI[k] + i - k];
        C[diagI[j]] += C[diagI[j] + i - j] * LiD[diagI[k] + i - j];          // update the jth diagonal element
        
      }
      for (i = *nC - 2; i > j; i--){
        temp = LiD[diagI[i]] * LiD[diagI[j] + i - j];
        C[diagI[j] + i - j] += temp;                                   // update the (i,j)th element
        for (k = i - 1; k >= j; k--)                                   // update rows j , ..., i - 1 of the jth column
          C[diagI[j] + k - j] += temp * LiD[diagI[k] + i - k];
      }
    }
  }

  delete [] LiD;
  return;
}    // end of function chinv


// ==================================================================================
// Function to compute an inversion of a real symmetric, positive definite matrix
//  from its Cholesky decomposition stored as a lower triangle in an array.
// This version returns both the full inversion and the inversion 
//  of the Cholesky decomposition only
// ==================================================================================
// INPUT PARAMETERS:
//
//    C ............. input matrix (it must be of length 0.5 * nC * (nC+1))
//                    = Cholesky decomposition of the original matrix obtained by the function 'cholesky'
//                    !!! C is rewritten by its inversion
//                    !!! It is assumed that C is positive definite, it is not checked!!!
//                    !!! Check this using the output from 'cholesky' function.
//    nC ............ number of rows (and columns) of matrix C
//    diagI ... indeces of diagonal elements of C in the array
//              (diagI[j] = (j * (2*(*nC) - j + 1)) / 2)
//              -> array of length nC
//    onlyCholInv ... 1 if only inversion of the Cholesky decomposition (L matrix) is desired
//
//  OUTPUT PARAMETERS:
//    C ....... inversion of C (its lower triangle)
//    ichol ... inversion of the Cholesky decomposition (its lower triangle)
//
void 
chinv2(double* C , double* ichol, const int* nC,  const int* diagI)
{

  double temp;
  int i, j, k;


  // * divide each column of the Cholesky decomposition under the diagonal by the diagonal element
  for (j = 0; j < *nC - 1; j++){
    if (C[diagI[j]] != 0){
      for (i = j + 1; i < *nC; i++ ){
        C[diagI[j] + i - j] /= C[diagI[j]];
      }
    }
  }

  // * (almost) invert the cholesky in the lower triangle
  // * do not multiply the rows in the lower strict triangle by the inverse of the diagonal element on that row
  // * this corresponds to L^{-1} where L had ones on the diagonal,
  for (j = 0; j < *nC; j++){                   // loop over columns
    if (C[diagI[j]] > 0) {
      C[diagI[j]] = 1/(C[diagI[j]]);                            /*this line inverts a diagonal */
      for (i= (j+1); i < *nC; i++){
	C[diagI[j] + i - j] = -C[diagI[j] + i - j];
  	for (k = 0; k < j; k++){                                /*sweep operator */
          C[diagI[k] + i - k] += C[diagI[j] + i - j]*C[diagI[k] + j - k];
        }
      }
    }
  }

  // * inverse of the Cholesky decomposition
  for (i = 0; i < *nC; i++){        // loop over rows (each element must be multiplied by the corresponding diagonal elements)
    if (C[diagI[i]] == 0) 
      for (j = 0 ; j <= i; j++) ichol[diagI[j] + i - j] = 0;              /* singular row */
    else{                  
      for (j = 0 ; j < i; j++) ichol[diagI[j] + i - j] = C[diagI[j] + i - j] * C[diagI[i]];
      ichol[diagI[i]] = C[diagI[i]];
    }
  }


  // * lower triangle now contains inverse of cholesky
  // * calculate (L^{-1})'DL^{-1} (inverse of cholesky decomp process) 
  //   where D is a diagonal matrix with 1/(diag(L)*diag(L))
  //   and L had ones on a diagonal
  // * this will be the inverse of the original matrix

  // Store L^{-1}D matrix in a new array
  double* LiD = new double[((*nC)*(*nC + 1))/2];
  for (j = 0; j < *nC; j++){
    LiD[diagI[j]] = C[diagI[j]] * C[diagI[j]];
    for (i = j + 1; i < *nC; i++)
      LiD[diagI[j] + i - j] = C[diagI[j] + i - j];
  }

  // Compute L^{-1}DL^{-1}
  for (j= 0; j < *nC; j++) {             // loop over columns
    if (C[diagI[j]] == 0){                         /* singular column */
      for (i = j; i < *nC; i++) C[diagI[j] + i - j] = 0;
    }
    else {
      C[diagI[j]] = LiD[diagI[j]];                                           // initialize the jth diagonal element
      if (*nC - 1 > j){
        i = *nC - 1;
        C[diagI[j] + i - j] = LiD[diagI[i]] * LiD[diagI[j] + i - j];         // last row of the jth column
        for (k = i - 1; k > j; k--)                                          // initialize rows j + 1, ..., n - 2 of the jth column
          C[diagI[j] + k - j] = C[diagI[j] + i - j] * LiD[diagI[k] + i - k];
        C[diagI[j]] += C[diagI[j] + i - j] * LiD[diagI[k] + i - j];          // update the jth diagonal element
        
      }
      for (i = *nC - 2; i > j; i--){
        temp = LiD[diagI[i]] * LiD[diagI[j] + i - j];
        C[diagI[j] + i - j] += temp;                                   // update the (i,j)th element
        for (k = i - 1; k >= j; k--)                                   // update rows j , ..., i - 1 of the jth column
          C[diagI[j] + k - j] += temp * LiD[diagI[k] + i - k];
      }
    }
  }

  delete [] LiD;
  return;
}    // end of function chinv


// =======================================================================
// Function to compute a cholesky decomposition of a matrix,
//  if the matrix is not positive definite, the function will successively try
//  to add a multiple of eps to the diagonal to make this matrix
//  positive definite

// INPUT PARAMETERS:
//    C .......... input matrix (it must be of length 0.5 * nC * (nC+1))
//                 C is NOT rewritten by its decomposition
//    nC ......... number of rows (and columns) of matrix C
//    diagI ...... indeces of diagonal elements of C in the array
//                 (diagI[j] = (j * (2*(*nC) - j + 1)) / 2)
//                 -> array of length nC
//    toler ...... threshhold value for detecting singularity
//    eps ........ coefficient which is added to a diagonal of C to get a positive definite matrix
//    nattempt ... maximal number of attempts when it is tried to add eps to a diagonal of C
//
//  OUTPUT PARAMETERS:
//    cholC ..... last Cholesky decomposition of C.
//                * the factorization is returned in the lower triangle of the original matrix,
//                * D occupies the diagonal 
//    rankC ..... the rank of the last C matrix (non-negative definite), 
//                or -rank if not SPD or NND
//    attempt ... how many times eps was added to a diagonal of C to get a positive definite matrix
//
void
chposDef(const double* C,      double* cholC,      int* rankC,           
         int* attempt,         const int* nC,      const int* diagI,  
         const double* toler,  const double* eps,  const int* nattempt)
{
  int i, j, k;

  for (k = 0; k < (*nC * (*nC + 1))/2; k++) cholC[k] = C[k];    

  *rankC = 0;
  *attempt = 0;
  while ((*rankC) < (*nC) && (*attempt) < (*nattempt)){
    cholesky(cholC, rankC, nC, diagI, toler);
    if (*rankC < *nC){                               // add a multiple of eps to the diagonal
      (*attempt)++;
      for (j = 0; j < *nC; j++){
        cholC[diagI[j]] = C[diagI[j]] + (*attempt)*(*eps);
        for (i = j + 1; i < *nC; i++)
          cholC[diagI[j] + i - j] = C[diagI[j] + i - j];
      }              
    }
  }

}    // end of function chPosDef

}    // end of extern "C"
