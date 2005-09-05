// Class to store a covariance matrix together with the information
//   on the prior distribution and functions to update it
// * used primarily for covariance matrix of normal random effects in 'bayessurvreg2'
// * it should replace 'covMatrix' class in all future versions of bayessurvreg
//
// 28/01/2005: basics
//             'update_after_change_icovm'
//             'update_after_change_covm'
//

#include "classCovMatrix.h"

/***** Nonparametric constructor *****/
CovMatrix::CovMatrix()
  : _nrow(0), _larray(0), _rank(0), _det(0),
    _type_prior(0), _dfD(0), _scaleD(NULL),
    _covm(NULL), _ichicovm(NULL), _icovm(NULL), _diagI(NULL),
    _qr(NULL), _qraux(NULL), _jpvt(NULL)
{
}


/***** Parametrix constructor *****/
CovMatrix::CovMatrix(const int* parmI,  const double* parmD)
{
  int i;

  /** Set up integer pointers **/
  int pnrow       = 0;
  int ptype_prior = pnrow + 1;
  // int pnext = ptype_prior + 1;

  if (parmI[pnrow] <= 0){
    _nrow = _larray = _rank = _type_prior = 0;
    _det = 0.0;
    _scaleD = _covm = _ichicovm = _icovm = _qr = _qraux = NULL;
    _diagI = _jpvt = NULL;
  }
  else{
    _nrow   = parmI[pnrow];
    _larray = (_nrow * (_nrow + 1))/2;

    /** Set up double pointers **/
    int pcovm   = 0;
    int pdfD    = pcovm + _larray;
    int pscaleD = pdfD + 1;
    // int pnext = pscaleD + _larray;

    /** Prior information **/
    _scaleD = (double*) calloc(_larray, sizeof(double));
    if (!_scaleD) throw returnR("Not enough memory available in CovMatrix constructor (_scaleD)", 1);
    switch (parmI[ptype_prior]){
    case InvWishart:
      _type_prior = InvWishart;
      if (parmD[pdfD] <= _nrow - 1){
        REprintf("_nrow = %d, df = %f\n", _nrow, parmD[pdfD]);
        throw returnR("Error: df for the covariance matrix must be > nrow - 1", 1);
      }
      _dfD = parmD[pdfD];
      for (i = 0; i < _larray; i++) _scaleD[i] = parmD[pscaleD + i];
      break;
    case SDUniform:
      if (_nrow > 1) throw returnR("Error: SDUniform prior not implemented for matrices, only for scalars", 1);
      _type_prior = SDUniform;
      if (parmD[pscaleD] <= 0){
        REprintf("upper limit = %f\n", parmD[pscaleD]);
        throw returnR("Error: Upper limit of the uniform prior for std. dev.(b) <= 0.", 1);
      }
      _scaleD[0] = parmD[pscaleD];
      _dfD       = 1/(_scaleD[0]*_scaleD[0]);
      break;
    default:
      throw returnR("Error: Unimplemented prior appeared in CovMatrix constructor", 1);
    }

    /** Components of the matrix **/
    _covm     = (double*) calloc(_larray, sizeof(double));
    _ichicovm = (double*) calloc(_larray, sizeof(double));
    _icovm    = (double*) calloc(_larray, sizeof(double));
    if (!_covm || !_ichicovm || !_icovm) throw returnR("Not enough memory available in CovMatrix constructor (_covm/_ichicovm/_icovm)", 1);

    _diagI = (int*) calloc(_nrow, sizeof(int));
    if (!_diagI) throw returnR("Not enough memory available in CovMatrix constructor (_diagI)", 1);
    for (i = 0; i < _nrow; i++) _diagI[i] =  (i * (2*_nrow - i + 1)) / 2;

    _qr    = (double*) calloc(_nrow*_nrow, sizeof(double));
    _qraux = (double*) calloc(_nrow, sizeof(double));
    _jpvt  = (int*) calloc(_nrow, sizeof(int));
    if (!_qr || !_qraux || !_jpvt) throw returnR("Not enough memory available in CovMatrix constructor (_qr/_qraux/_jpvt)", 1);

    for (i = 0; i < _larray; i++) _covm[i] = parmD[pcovm + i];
    update_after_change_covm();
  }
}    /** end of the parametric constructor **/


/***** Copy constructor *****/
CovMatrix::CovMatrix(const CovMatrix& cm)
{
  int i;
  
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
    if (!_scaleD) throw returnR("Not enough memory available in CovMatrix copy constructor (_scaleD)", 1);

    _covm     = (double*) calloc(_larray, sizeof(double));
    _ichicovm = (double*) calloc(_larray, sizeof(double));
    _icovm    = (double*) calloc(_larray, sizeof(double));
    if (!_covm || !_ichicovm || !_icovm) throw returnR("Not enough memory available in CovMatrix copy const. (_covm/_ichicovm/_icovm)", 1);
    for (i = 0; i < _larray; i++){
      _scaleD[i] = cm._scaleD[i];
      _covm[i]   = cm._covm[i];
      _ichicovm[i]   = cm._ichicovm[i];
      _icovm[i]   = cm._icovm[i];
    }

    _diagI = (int*) calloc(_nrow, sizeof(int));
    if (!_diagI) throw returnR("Not enough memory available in CovMatrix copy constructor (_diagI)", 1);

    _qr    = (double*) calloc(_nrow*_nrow, sizeof(double));
    _qraux = (double*) calloc(_nrow, sizeof(double));
    _jpvt  = (int*) calloc(_nrow, sizeof(int));
    if (!_qr || !_qraux || !_jpvt) throw returnR("Not enough memory available in CovMatrix copy constructor (_qr/_qraux/_jpvt)", 1);
    for (i = 0; i < _nrow; i++){
      _diagI[i] = cm._diagI[i];
      _qraux[i] = cm._qraux[i];
      _jpvt[i]  = cm._jpvt[i];
    }
    for (i = 0; i < _nrow*_nrow; i++) _qr[i] = cm._qr[i];
  }
}    /** end of the copy constructor **/


/***** Assignment operator *****/
CovMatrix&
CovMatrix::operator=(const CovMatrix& cm)
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
    if (!_scaleD) throw returnR("Not enough memory available in CovMatrix assignment operator (_scaleD)", 1);

    _covm     = (double*) calloc(_larray, sizeof(double));
    _ichicovm = (double*) calloc(_larray, sizeof(double));
    _icovm    = (double*) calloc(_larray, sizeof(double));
    if (!_covm || !_ichicovm || !_icovm) throw returnR("Not enough memory available in CovMatrix assign. op. (_covm/_ichicovm/_icovm)", 1);
    for (i = 0; i < _larray; i++){
      _scaleD[i] = cm._scaleD[i];
      _covm[i]   = cm._covm[i];
      _ichicovm[i]   = cm._ichicovm[i];
      _icovm[i]   = cm._icovm[i];
    }

    _diagI = (int*) calloc(_nrow, sizeof(int));
    if (!_diagI) throw returnR("Not enough memory available in CovMatrix assignment operator (_diagI)", 1);

    _qr    = (double*) calloc(_nrow*_nrow, sizeof(double));
    _qraux = (double*) calloc(_nrow, sizeof(double));
    _jpvt  = (int*) calloc(_nrow, sizeof(int));
    if (!_qr || !_qraux || !_jpvt) throw returnR("Not enough memory available in CovMatrix assignment operator (_qr/_qraux/_jpvt)", 1);
    for (i = 0; i < _nrow; i++){
      _diagI[i] = cm._diagI[i];
      _qraux[i] = cm._qraux[i];
      _jpvt[i]  = cm._jpvt[i];
    }
    for (i = 0; i < _nrow*_nrow; i++) _qr[i] = cm._qr[i];
  }

  return *this;
}    /** end of the assignment operator **/


/***** Destructor *****/
CovMatrix::~CovMatrix()
{
  if (_nrow){
    free(_scaleD);    free(_covm);    free(_ichicovm);    free(_icovm);
    free(_diagI);     free(_qr);      free(_qraux);       free(_jpvt);
  }
}


/***** print *****/
void
CovMatrix::print() const
{
  Rprintf("\nCovariance matrix object:\n");

  if (!_nrow){
    Rprintf("   Empty CovMatrix.\n");
    return;
  }

  int j;
  Rprintf("   nrow = %d,  larray = %d\n", _nrow, _larray);
  Rprintf("   rank = %d,  det = %g\n", _rank, _det);
  char zz[15];
  switch (_type_prior){
  case InvWishart : strcpy(zz, "Inverse-Wishart"); break;
  case SDUniform  : strcpy(zz, "SD Uniform"); break;
  default         : strcpy(zz, "unimplemented"); break;  
  }
  Rprintf("   Prior distribution = %s\n", zz);
  Rprintf("   df(prior) = %g\n", _dfD);
  Rprintf("   scale(prior) = "); for (j = 0; j < _larray; j++) Rprintf("%g,  ", _scaleD[j]);
  Rprintf("\n   covm = ");     for (j = 0; j < _larray; j++) Rprintf("%g,  ", _covm[j]);  
  Rprintf("\n   ichicovm = "); for (j = 0; j < _larray; j++) Rprintf("%g,  ", _ichicovm[j]);  
  Rprintf("\n   icovm = ");    for (j = 0; j < _larray; j++) Rprintf("%g,  ", _icovm[j]);  
  Rprintf("\n   diagI = ");    for (j = 0; j < _nrow; j++) Rprintf("%d,  ", _diagI[j]);

  return;
}


/***** CovMatrix2initArray *****/
void
CovMatrix::CovMatrix2initArray(int* parmI,  double* parmD) const
{
  int i;

  /** Set up integer pointers **/
  int pnrow       = 0;
  int ptype_prior = pnrow + 1;
  // int pnext = ptype_prior + 1;

  /** Set up double pointers **/
  int pcovm   = 0;
  int pdfD    = pcovm + _larray;
  int pscaleD = pdfD + 1;
  // int pnext = pscaleD + _larray;

  parmI[pnrow]       = _nrow;
  parmI[ptype_prior] = _type_prior;

  for (i = 0; i < _larray; i++){
    parmD[pcovm + i]   = _covm[i];
    parmD[pscaleD + i] = _scaleD[i];
  }
  parmD[pdfD] = _dfD;
 
  return;
}


/***** update_after_change_icovm  *****/
//
// * update the additional information on the covariance matrix
//   when the elements of the array icovm change
// * member function 'updateFromInv' of the original class 'covMatrix'
//
void
CovMatrix::update_after_change_icovm()
{
  int i, j;

  if (_nrow > 0){

    // Cholesky decomposition of the inverse, store it in covm
    for (i = 0; i < _larray; i++) _covm[i] = _icovm[i];
    cholesky(_covm, &_rank, &_nrow, _diagI, &_toler_chol_CovMatrix);
    if (_rank < 0) 
      throw returnR("Error: Non positive-semidefinite covariance matrix in CovMatrix::update_after_change_icovm", 1);

    // Inverse of the Cholesky decomposition of the inverse and the full inverse of the inverse    
    chinv2(_covm, _ichicovm, &_nrow, _diagI);

    // QR decomposition and determinant
    for (i = 0; i < _nrow; i++){
      _jpvt[i] = i;
      _qr[_nrow*i + i] = _covm[_diagI[i]];
      for (j = i+1; j < _nrow; j++){
        _qr[_nrow*i + j] = _covm[_diagI[i] + j - i];
        _qr[_nrow*j + i] = _covm[_diagI[i] + j - i];
      }
    }
    dqrdc2CPP(_qr, &_nrow, &_nrow, &_toler_QR_CovMatrix, &_rank, _qraux, _jpvt);
    if (_rank < _nrow) 
      _det = 0.0; 
    else{
      _det = _qr[0];
      for (i = 1; i < _nrow; i++) _det *= _qr[_nrow*i + i];
      if (_nrow % 2 == 0) _det *= (-1);
    }
  }

  return;
}


/*****  update_after_change_covm  *****/
//
// * update the additional information on the covariance matrix
//   when the elements of the array _covm change
// * member function 'update' of the original class 'covMatrix'
//
void
CovMatrix::update_after_change_covm()
{
  int i, j;
 
  if (_nrow > 0){

    // Cholesky decomposition 
    for (i = 0; i < _larray; i++) _ichicovm[i] = _covm[i];
    cholesky(_ichicovm, &_rank, &_nrow, _diagI, &_toler_chol_CovMatrix);
    if (_rank < 0) 
      throw returnR("Error: Non positive-semidefinite covariance matrix in CovMatrix::update_after_change_covm", 1);

    // Inverse
    for (i = 0; i < _larray; i++) _icovm[i] = _ichicovm[i];  
    if (_rank == _nrow)  
      chinv(_icovm, &_nrow, _diagI, &ZERO_INT);
    else
      throw returnR("Error: Singular covariance matrix in CovMatrix::update_after_change_covm", 1);

    // Cholesky decomposition of the inverse and its inverse
    for (i = 0; i < _larray; i++) _ichicovm[i] = _icovm[i];  
    cholesky(_ichicovm, &_rank, &_nrow, _diagI, &_toler_chol_CovMatrix);
    chinv(_ichicovm, &_nrow, _diagI, &ONE_INT);

    // QR decomposition and determinant
    for (i = 0; i < _nrow; i++){
      _jpvt[i] = i;
      _qr[_nrow*i + i] = _covm[_diagI[i]];
      for (j = i+1; j < _nrow; j++){
        _qr[_nrow*i + j] = _covm[_diagI[i] + j - i];
        _qr[_nrow*j + i] = _covm[_diagI[i] + j - i];
      }
    }
    dqrdc2CPP(_qr, &_nrow, &_nrow, &_toler_QR_CovMatrix, &_rank, _qraux, _jpvt);
    if (_rank < _nrow) 
      _det = 0.0; 
    else{
      _det = _qr[0];
      for (i = 1; i < _nrow; i++) _det *= _qr[_nrow*i + i];
      if (_nrow % 2 == 0) _det *= (-1);
    }
  }

  return;
}


