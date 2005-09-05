// Class to store a covariance matrix together with the information
//   on the prior distribution and functions to update it
// * used primarily for covariance matrix of random effects in 'bayessurvreg2'
// * it should replace 'covMatrix' class in all future versions of bayessurvreg
//
#ifndef _CLASS_COV_MATRIX_H_
#define _CLASS_COV_MATRIX_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "constants.h"
#include "cholesky.h"
#include "qrdecomp.h"

const double _toler_chol_CovMatrix = 1e-10;
const double _toler_QR_CovMatrix = 1e-10;

class CovMatrix
{
  protected:
    int _nrow;              // number of rows
    int _larray;            // length of the lower triangle = (_nrow * (_nrow + 1))/2
    int _rank;              // rank of the matrix
    double _det;            // determinant of the matrix

    int _type_prior;        // type of the prior distribution, 0 = Inverse-Wishart, 1 = unifor on standard deviation (only if _nrow==1)
    double _dfD;            // _type_prior = InvWishart: prior degrees of freedoom
                            //               SDUniform:  1/range^2
    double* _scaleD;        // [_larray] _type_prior = InvWishart: prior scale matrix (lower triangle)
                            //                         SDUniform: upper limit for the standard deviation
    
    double* _covm;           // [_larray] lower triangle of the matrix
    double* _ichicovm;       // [_larray] lower triangle of the inverse of the Cholesky decomposition of the inverse of the covariance mat.
    double* _icovm;          // [_larray] lower triangle of the inverse of the matrix
    int* _diagI;             // [_nrow]   indeces of diagonal elements in above arrays

    double* _qr;             // [_nrow*_nrow] QR decomposition as returned by the function dqrdc2CPP
    double* _qraux;          // [_nrow]       further information required to recover the orthogonal part of the QR decomposition
    int* _jpvt;              // [_nrow]       pivot information from QR decomposition

 public:
    CovMatrix();
    CovMatrix(const int* parmI,  const double* parmD);
    CovMatrix(const CovMatrix& cm);
    CovMatrix& 
    operator=(const CovMatrix& cm);
    ~CovMatrix();

    void print() const;

    void CovMatrix2initArray(int* parmI,  double* parmD) const;

    inline int
    nrow() const { return _nrow;}

    inline int
    larray() const { return _larray;}

    inline int
    rank() const { return _rank;}

    inline double
    det() const { return _det;}

    inline bool
    isNull() const { if (_nrow == 0) return true;
                     else            return false;}

    inline void
    new_covm(const int& i, const double& dd)
    {
     if (i >= _larray || i < 0) throw returnR("C++ Error: Incorrect i in CovMatrix::new_covm(i, dd).", 1);
     _covm[i] = dd;
     return;
    }

    inline double
    operator()(const int& i, const int& j) const
    {
      if (i < 0 || i >= _nrow) throw returnR("Index i out of range in CovMatrix operator()", 1);
      if (j < 0 || j >= _nrow) throw returnR("Index j out of range in CovMatrix operator()", 1);
      if (i >= j) return _covm[_diagI[j] + i - j];
      else        return _covm[_diagI[i] + j - i];
    }

    inline double
    dfD() const { return _dfD; }

    inline double
    scaleD(const int& i) const
    {
      if (i < 0 || i >= _larray) throw returnR("Index i out of range in CovMatrix.scaleD(i)", 1);
      return _scaleD[i];
    }

    inline double
    covm(const int& i) const
    {
      if (i < 0 || i >= _larray) throw returnR("Index i out of range in CovMatrix.covm(i)", 1);
      return _covm[i];
    }

    inline double
    ichicovm(const int& i) const
    {
      if (i < 0 || i >= _larray) throw returnR("Index i out of range in CovMatrix.ichicovm(i)", 1);
      return _ichicovm[i];
    }

    inline double
    icovm(const int& i) const
    {
      if (i < 0 || i >= _larray) throw returnR("Index i out of range in CovMatrix.icovm(i)", 1);
      return _icovm[i];
    }

    inline int
    diagI(const int& i) const
    {
      if (i < 0 || i >= _nrow) throw returnR("Index i out of range in CovMatrix.diagI(i)", 1);
      return _diagI[i];
    }

    inline double*
    scaleDP() const{ return _scaleD; }

    inline double*
    covmP() const{ return _covm; }

    inline double*
    ichicovmP() const{ return _ichicovm; }

    inline double*
    icovmP() const{ return _icovm; }

    inline int*
    diagIP() const{ return _diagI; }

    void
    update_after_change_icovm();

    void
    update_after_change_covm();





};    // end of the class CovMatrix

#endif
