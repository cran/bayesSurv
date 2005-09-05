// Derived class from CovMatrix
// * added method to update the covariance matrix under the assumption that 
//   the random effects are normal
// * again, a chicken-egg problem is present if we don't use a derived class
//
// 29/01/2005: 'GIBBSnormalRE'
//
#include "classCovMatrixExtend.h"

/***** Nonparametric constructor *****/
CovMatrixExtend::CovMatrixExtend() 
  : CovMatrix()
{}


/***** Parametric constructor  *****/
CovMatrixExtend::CovMatrixExtend(const int* parmI,  const double* parmD) 
  : CovMatrix(parmI, parmD)
{}


/***** Copy constructor 1 *****/
CovMatrixExtend::CovMatrixExtend(const CovMatrixExtend& cm) 
  : CovMatrix(cm)
{}


/***** Copy constructor 2 *****/
CovMatrixExtend::CovMatrixExtend(const CovMatrix& cm) 
  : CovMatrix(cm)
{}

/***** Assignment operator *****/
CovMatrixExtend&
CovMatrixExtend::operator=(const CovMatrixExtend& cm)
{
  int i;

  if (_nrow){
    free(_scaleD);    free(_covm);    free(_ichicovm);    free(_icovm);
    free(_diagI);     free(_qr);      free(_qraux);       free(_jpvt);
  }

  if (cm._nrow == 0){
    _nrow = _larray = _rank = _type_prior = 0;
    _det = 0.0;
    _scaleD = _covm = _ichicovm = _icovm = _qr = _qraux = NULL;
    _diagI = _jpvt = NULL;
  }
  else{
    _nrow   = cm._nrow;
    _larray = cm._larray;
    _rank   = cm._rank;
    _det    = cm._det;

    _type_prior = cm._type_prior;
    _dfD        = cm._dfD;
    _scaleD = (double*) calloc(_larray, sizeof(double));
    if (!_scaleD) throw returnR("Not enough memory available in CovMatrixExtend assignment operator (_scaleD)", 1);

    _covm     = (double*) calloc(_larray, sizeof(double));
    _ichicovm = (double*) calloc(_larray, sizeof(double));
    _icovm    = (double*) calloc(_larray, sizeof(double));
    if (!_covm || !_ichicovm || !_icovm) throw returnR("Not enough memory avail. in CovMatrixExtend ass. op. (_covm/_ichicovm/_icovm)", 1);
    for (i = 0; i < _larray; i++){
      _scaleD[i] = cm._scaleD[i];
      _covm[i]   = cm._covm[i];
      _ichicovm[i]   = cm._ichicovm[i];
      _icovm[i]   = cm._icovm[i];
    }

    _diagI = (int*) calloc(_nrow, sizeof(int));
    if (!_diagI) throw returnR("Not enough memory available in CovMatrixExtend assignment operator (_diagI)", 1);

    _qr    = (double*) calloc(_nrow*_nrow, sizeof(double));
    _qraux = (double*) calloc(_nrow, sizeof(double));
    _jpvt  = (int*) calloc(_nrow, sizeof(int));
    if (!_qr || !_qraux || !_jpvt) throw returnR("Not enough memory avail. in CovMatrixExtend assignment operator (_qr/_qraux/_jpvt)", 1);
    for (i = 0; i < _nrow; i++){
      _diagI[i] = cm._diagI[i];
      _qraux[i] = cm._qraux[i];
      _jpvt[i]  = cm._jpvt[i];
    }
    for (i = 0; i < _nrow*_nrow; i++) _qr[i] = cm._qr[i];
  }

  return *this;
}    /** end of the assignment operator **/


//
// ***** GIBBSnormalRE *****
// 
// * update the covariance matrix of normal random effects using the Gibbs move
//
void
CovMatrixExtend::GIBBSnormalRE(const RandomEff* b_obj,  const BetaGamma* bbg)
{
  if (!_nrow) return;

  static int i, rank;
  static double df, shape, rate;

    /** Compute sum of squares of random effects, store it in Dcm->covm **/
  b_obj->sumSquare(_covm, bbg);

  switch (_type_prior){
  case InvWishart:

      /** Degrees of freedom for the full conditional **/
    df = _dfD + b_obj->nCluster();

      /** Compute a scale matrix of the inverse-Wishart distribution which is a full conditional,  **/
      /**   store it in _covm                                                                      **/
    for (i = 0; i < _larray; i++) _covm[i] += _scaleD[i];

      /** Cholesky decomposition of a scale matrix                                                 **/
      /**   and an inverse of this Cholesky decomposition                                          **/
    cholesky(_covm, &rank, &_nrow, _diagI, &_toler_chol_CovMatrix);
    if (rank < _nrow){
      REprintf("nrow = %d,  rank = %d\n", _nrow, rank);
      throw returnR("Trap: Scale matrix for update of D is not positive definite.", 1);
    }
    chinv(_covm, &_nrow, _diagI, &ONE_INT);

      /** Sample new inverse D matrix                                                              **/
    rwishart2(_icovm, &_nrow, &df, _covm, _diagI, &ONE_INT, &ZERO_INT);
    break;

  case SDUniform:
      /** Shape and rate of the full conditional distribution   **/
    shape = 0.5 * (b_obj->nCluster() - 1);  
    rate = 0.5 * _covm[0];
    
      /** Sample new inverse var(b), i.e. new inverse D matrix  **/
    rltruncGamma(_icovm, &shape, &rate, &_dfD, &ONE_INT, &ZERO_INT);     
    break;
  }

    /** Compute D, and iLi where iLi = L^{-1} and L*L' = D^{-1} **/
  this->update_after_change_icovm();
  return;
}    /*** end of the function CovMatrixExtend::GIBBSnormalRE ***/

