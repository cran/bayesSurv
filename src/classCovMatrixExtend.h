// Derived class from CovMatrix
// * added method to update the covariance matrix under the assumption that 
//   the random effects are normal
// * again, a chicken-egg problem is present if we don't use a derived class
//
#ifndef _CLASS_COV_MATRIX_EXTEND_H_
#define _CLASS_COV_MATRIX_EXTEND_H_

#include "random.h"
#include "mvtdist.h"

#include "classCovMatrix.h"
#include "classRandomEff.h"
#include "classBetaGamma.h"

class CovMatrixExtend : public CovMatrix
{
  private:

  public:
    CovMatrixExtend();
    CovMatrixExtend(const int* parmI,  const double* parmD);
    CovMatrixExtend(const CovMatrixExtend& cm);
    CovMatrixExtend(const CovMatrix& cm);
    CovMatrixExtend& operator=(const CovMatrixExtend& cm);

    void
    GIBBSnormalRE(const RandomEff* b_obj,  const BetaGamma* bbg);

};    /*** end of the class CovMatrixExtend ***/


#endif

