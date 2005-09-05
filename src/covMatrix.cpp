// Class to store covariance matrices (or other symmetric positive definite matrices)
//
// 30/01/2004: start working on it
//
#include "covMatrix.h"

using namespace std;

// Non-parametric constructor
// ===========================
covMatrix::covMatrix()
  : _nrow(0), _larray(0), _rank(0), det(0)
{
  covm = NULL;
  ichicovm = new double[_larray];
  icovm = new double[_larray];
  diagI = new int[_nrow];

  qr = new double[_nrow*_nrow];
  qraux = new double[_nrow];

  if (ichicovm == NULL || icovm == NULL || diagI == NULL || qr == NULL || qraux == NULL)
    throw returnR("C++ Error: 'covMatrix::covMatrix()' could not allocate a memory for a working space.", 1);
}


// Copy constructor
// =================
covMatrix::covMatrix(const covMatrix& cm)
  : _nrow(cm._nrow), _larray(cm._larray), _rank(cm._rank), det(cm.det)
{
  int i, j;

  covm = cm.covm;
  ichicovm = new double[_larray];
  icovm = new double[_larray];
  diagI = new int[_nrow];

  qr = new double[_nrow*_nrow];
  qraux = new double[_nrow];
  jpvt = new int[_nrow];

  if (ichicovm == NULL || icovm == NULL || diagI == NULL || qr == NULL || qraux == NULL || jpvt == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  
  for (i = 0; i < _larray; i++){
    ichicovm[i] = cm.ichicovm[i];
    icovm[i] = cm.icovm[i];
  }
  for (i = 0; i < _nrow; i++){
    diagI[i] = cm.diagI[i];
    qraux[i] = cm.qraux[i];
    for (j = 0; j < _nrow; j++) qr[i*_nrow + j] = cm.qr[i*_nrow + j];
  }
}
  

// Assignment operator
// ===================
covMatrix&
covMatrix::operator=(const covMatrix& cm)
{
  int i, j;

  delete [] ichicovm;
  delete [] icovm;
  delete [] diagI;

  _nrow = cm._nrow;
  _larray = cm._larray;
  _rank = cm._rank;
  det = cm.det;

  covm = cm.covm;
  ichicovm = new double[_larray];
  icovm = new double[_larray];
  diagI = new int[_nrow];

  qr = new double[_nrow*_nrow];
  qraux = new double[_nrow];
  jpvt = new int[_nrow];

  if (ichicovm == NULL || icovm == NULL || diagI == NULL || qr == NULL || qraux == NULL || jpvt == NULL)
    throw returnR("C++ Error: 'covMatrix' operator= could not allocate a memory for a working space.", 1);
  
  for (i = 0; i < _larray; i++){
    ichicovm[i] = cm.ichicovm[i];
    icovm[i] = cm.icovm[i];
  }
  for (i = 0; i < _nrow; i++){
    diagI[i] = cm.diagI[i];
    qraux[i] = cm.qraux[i];
    jpvt[i] = cm.jpvt[i];
    for (j = 0; j < _nrow; j++) qr[i*_nrow + j] = cm.qr[i*_nrow + j];
  }

  return *this;
}


// Destructor
// ==========
covMatrix::~covMatrix()
{
  delete [] ichicovm;
  delete [] icovm;
  delete [] diagI;
  delete [] qr;
  delete [] qraux;
  delete [] jpvt;
}


// Parametric constructor
// =======================
covMatrix::covMatrix(double* covmA, const int* nr, const double* tolChol, const double* tolQR)
  : _nrow(*nr), _rank(0), det(0)
{
  int i, j;
  _larray = (_nrow * (_nrow + 1))/2;
 
  if (_nrow == 0) covm = NULL;
  else            covm = covmA;
  ichicovm = new double[_larray];
  icovm = new double[_larray];
  diagI = new int[_nrow];

  qr = new double[_nrow*_nrow];
  qraux = new double[_nrow];
  jpvt = new int[_nrow];

  if (ichicovm == NULL || icovm == NULL || diagI == NULL || qr == NULL || qraux == NULL || jpvt == NULL)
    throw returnR("C++ Error: 'covMatrix' constructor could not allocate a memory for a working space.", 1);

  // Indeces of diagonal elements in an array
  for (i = 0; i < _nrow; i++) diagI[i] =  (i * (2*_nrow - i + 1)) / 2;

  if (_nrow > 0){
    // Cholesky decomposition 
    for (i = 0; i < _larray; i++) ichicovm[i] = covm[i];
    cholesky(ichicovm, &_rank, &_nrow, diagI, tolChol);
    if (_rank < 0) 
      throw returnR("C++ Error: Non positive-semidefinite covariance matrix appeared somewhere", 1);    

    // Inverse
    for (i = 0; i < _larray; i++) icovm[i] = ichicovm[i];  
    if (_rank == _nrow)  
      chinv(icovm, &_nrow, diagI, &ZERO_INT);
    else
      throw returnR("C++ Error: Singular covariance matrix appeared somewhere", 1);    

    // Cholesky decomposition of the inverse and its inverse
    for (i = 0; i < _larray; i++) ichicovm[i] = icovm[i];  
    cholesky(ichicovm, &_rank, &_nrow, diagI, tolChol);
    chinv(ichicovm, &_nrow, diagI, &ONE_INT);

    // QR decomposition and determinant
    for (i = 0; i < _nrow; i++){
      jpvt[i] = i;
      qr[_nrow*i + i] = covm[diagI[i]];
      for (j = i+1; j < _nrow; j++){
        qr[_nrow*i + j] = covm[diagI[i] + j - i];
        qr[_nrow*j + i] = covm[diagI[i] + j - i];
      }
    }
    dqrdc2CPP(qr, &_nrow, &_nrow, tolQR, &_rank, qraux, jpvt);
    if (_rank < _nrow) 
      det = 0.0; 
    else{
      det = qr[0];
      for (i = 1; i < _nrow; i++) det *= qr[_nrow*i + i];
      if (_nrow % 2 == 0) det *= (-1);
    }
  }
}


// print
// =====
void
covMatrix::print() const
{
  cout << "nrow = " << _nrow;
  cout << ",  larray = " << _larray;
  cout << ",  rank = " << _rank;
  cout << ",  det = " << det << endl;;
  cout << "covm = "; printArray(covm, &_larray);
  cout << "ichicovm = "; printArray(ichicovm, &_larray);
  cout << "diagI = "; printArray(diagI, &_nrow);
  int r_r = _nrow*_nrow;
  cout << "qr = "; printArray(qr, &r_r);
  cout << "qraux = "; printArray(qraux, &_nrow);
  cout << "jpvt = "; printArray(jpvt, &_nrow);

  return;
}

// updateFromInv
// =============
// * update the additional information on the covariance matrix
//   when the elements of the array icovm change
//
void
covMatrix::updateFromInv(const double* tolChol, const double* tolQR)
{
  int i, j;

  if (_nrow > 0){
    // Cholesky decomposition of the inverse, store it in covm
    for (i = 0; i < _larray; i++) covm[i] = icovm[i];
    cholesky(covm, &_rank, &_nrow, diagI, tolChol);
    if (_rank < 0) 
      throw returnR("C++ Error: Non positive-semidefinite covariance matrix appeared somewhere", 1);    

    // Inverse of the Cholesky decomposition of the inverse and the full inverse of the inverse    
    chinv2(covm, ichicovm, &_nrow, diagI);

    // QR decomposition and determinant
    for (i = 0; i < _nrow; i++){
      jpvt[i] = i;
      qr[_nrow*i + i] = covm[diagI[i]];
      for (j = i+1; j < _nrow; j++){
        qr[_nrow*i + j] = covm[diagI[i] + j - i];
        qr[_nrow*j + i] = covm[diagI[i] + j - i];
      }
    }
    dqrdc2CPP(qr, &_nrow, &_nrow, tolQR, &_rank, qraux, jpvt);
    if (_rank < _nrow) 
      det = 0.0; 
    else{
      det = qr[0];
      for (i = 1; i < _nrow; i++) det *= qr[_nrow*i + i];
      if (_nrow % 2 == 0) det *= (-1);
    }
  }

  return;
}


// update
// =======
// * update the additional information on the covariance matrix
//   when the elements of the array covm change
//
void
covMatrix::update(const double* tolChol, const double* tolQR)
{
  int i, j;
 
  if (_nrow > 0){
    // Cholesky decomposition 
    for (i = 0; i < _larray; i++) ichicovm[i] = covm[i];
    cholesky(ichicovm, &_rank, &_nrow, diagI, tolChol);
    if (_rank < 0) 
      throw returnR("C++ Error: Non positive-semidefinite covariance matrix appeared somewhere", 1);    

    // Inverse
    for (i = 0; i < _larray; i++) icovm[i] = ichicovm[i];  
    if (_rank == _nrow)  
      chinv(icovm, &_nrow, diagI, &ZERO_INT);
    else
      throw returnR("C++ Error: Singular covariance matrix appeared somewhere", 1);    

    // Cholesky decomposition of the inverse and its inverse
    for (i = 0; i < _larray; i++) ichicovm[i] = icovm[i];  
    cholesky(ichicovm, &_rank, &_nrow, diagI, tolChol);
    chinv(ichicovm, &_nrow, diagI, &ONE_INT);

    // QR decomposition and determinant
    for (i = 0; i < _nrow; i++){
      jpvt[i] = i;
      qr[_nrow*i + i] = covm[diagI[i]];
      for (j = i+1; j < _nrow; j++){
        qr[_nrow*i + j] = covm[diagI[i] + j - i];
        qr[_nrow*j + i] = covm[diagI[i] + j - i];
      }
    }
    dqrdc2CPP(qr, &_nrow, &_nrow, tolQR, &_rank, qraux, jpvt);
    if (_rank < _nrow) 
      det = 0.0; 
    else{
      det = qr[0];
      for (i = 1; i < _nrow; i++) det *= qr[_nrow*i + i];
      if (_nrow % 2 == 0) det *= (-1);
    }
  }

  return;
}
