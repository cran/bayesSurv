// ======================================================================================
// ***** Header for covMatrix.cpp: a class to store covariance matrices
// ======================================================================================
#ifndef _COV_MATRIX_H_
#define _COV_MATRIX_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "constants.h"
#include "cholesky.h"
#include "templatefun.h"
#include "printArray.h"
#include "qrdecomp.h"

class covMatrix
{
 private:
  int _nrow;          // number of rows
  int _larray;        // length of the lower triangle
  int _rank;          // rank of the matrix  

 public:
  double* covm;       // lower triangle of the matrix stored in an array
                      // This points always to an external array!!!

  double* ichicovm;   // lower triangle of the inverse of the Cholesky decomposition of the inverse of the covariance matrix stored in an array
  double* icovm;      // lower triangle of the inverse of this matrix stored in an array
  int* diagI;         // indeces of diagonal elements in above mentioned arrays

  double* qr;         // QR decomposition as returned by the function dqrdc2CPP
  double* qraux;      // further information required to recover the orthogonal part of the QR decomposition
  int* jpvt;          // pivot information from QR decomposition
  double det;         // determinant of the matrix

  covMatrix();

  covMatrix(const covMatrix& cm);  

  covMatrix(double* covmA, const int* nr, const double* tolChol, const double* tolQR);

  ~covMatrix();

  covMatrix& 
  operator=(const covMatrix& cm);

  inline int
  nrow() const { return _nrow;}

  inline int
  larray() const { return _larray;}

  inline int
  rank() const { return _rank;}

  inline bool
  isNull() const { if (_nrow == 0) return true;
                   else            return false;}

  void
  print() const;

  void
  updateFromInv(const double* tolChol, const double* tolQR);

  void
  update(const double* tolChol, const double* tolQR);

};   // end of the class covMatrix

#endif
