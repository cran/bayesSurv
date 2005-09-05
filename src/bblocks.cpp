// Implementation of the class bblocks
//
// 05/02/2004: start working on it
// 30/07/2004: added possibility to have a uniform prior for sd(b) if the random effect is univariate
//
#include "bblocks.h"

using namespace std;

// Non-parametric constructor
// ==========================
bblocks::bblocks()
  : _nRandom(0), _nBlocks(0), _nCluster(0)
{
    bM = NULL;
    priorD = NULL;
    dfD = NULL;
    scaleD = NULL;
    typeUpd = NULL;
    nInBlock = NULL;
    indBlock = new int*[0];
    diagI = new int*[0];
    covpar = new double*[0];
    chcovpar = new double*[0];
    weightUnif = NULL;
    halfRangeUnif = NULL;
    sumAccept = NULL;
}


// Copy constructor
// ================
bblocks::bblocks(const bblocks& bb)
 :  _nRandom(bb._nRandom), _nBlocks(bb._nBlocks), _nCluster(bb._nCluster)
{
  int i, j, dimi;

  bM = bb.bM;

  priorD = bb.priorD;
  dfD = bb.dfD;
  scaleD = bb.scaleD;
  typeUpd = bb.typeUpd;

  nInBlock = bb.nInBlock;

  indBlock = new int*[_nBlocks];
  if (indBlock == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    indBlock[i] = bb.indBlock[i];
  }

  diagI = new int*[_nBlocks];
  if (diagI == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    diagI[i] = new int[nInBlock[i]];
    if (diagI[i] == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (j = 0; j < nInBlock[i]; j++)
      diagI[i][j] = bb.diagI[i][j];
  }

  covpar = new double*[_nBlocks];
  chcovpar = new double*[_nBlocks];
  if (covpar == NULL || chcovpar == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    covpar[i] = bb.covpar[i];
    dimi = (nInBlock[i] * (nInBlock[i] + 1))/2;
    chcovpar[i] = new double[dimi];
    if (chcovpar[i] == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (j = 0; j < dimi; j++) chcovpar[i][j] = bb.chcovpar[i][j];
  }

  weightUnif = bb.weightUnif;

  halfRangeUnif = bb.halfRangeUnif;

  sumAccept = bb.sumAccept;
}


// Assignment operator
// ===================
bblocks&
bblocks::operator=(const bblocks& bb)
{
  int i, j, dimi;

  delete [] indBlock;
  for (i = 0; i < _nBlocks; i++) delete [] diagI[i];
  delete [] diagI;
  delete [] covpar;
  for (i = 0; i < _nBlocks; i++) delete [] chcovpar[i];
  delete [] chcovpar;
  
  _nRandom = bb._nRandom; 
  _nBlocks = bb._nBlocks;
  _nCluster = bb._nCluster;  

  bM = bb.bM;

  priorD = bb.priorD;
  dfD = bb.dfD;
  scaleD = bb.scaleD;

  typeUpd = bb.typeUpd;

  nInBlock = bb.nInBlock;

  indBlock = new int*[_nBlocks];
  if (indBlock == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    indBlock[i] = bb.indBlock[i];
  }

  diagI = new int*[_nBlocks];
  if (diagI == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    diagI[i] = new int[nInBlock[i]];
    if (diagI[i] == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (j = 0; j < nInBlock[i]; j++)
      diagI[i][j] = bb.diagI[i][j];
  }

  covpar = new double*[_nBlocks];
  chcovpar = new double*[_nBlocks];
  if (covpar == NULL || chcovpar == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    covpar[i] = bb.covpar[i];
    dimi = (nInBlock[i] * (nInBlock[i] + 1))/2;
    chcovpar[i] = new double[dimi];
    if (chcovpar[i] == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (j = 0; j < dimi; j++) chcovpar[i][j] = bb.chcovpar[i][j];
  }

  weightUnif = bb.weightUnif;

  halfRangeUnif = bb.halfRangeUnif;

  sumAccept = bb.sumAccept;
  return *this;
}


// Parametric constructor
// =======================
bblocks::bblocks(double* bbM,              int* parmI,   double* parmD,
                 int* sumAcceptP, 
                 const double* tolerChol)
{
  int i, j, helpi, dimi;

  _nRandom = parmI[0];

  if (_nRandom == 0){
    _nBlocks = 0;
    _nCluster = 0;

    bM = NULL;
    priorD = NULL;
    typeUpd = NULL;
    dfD = NULL;
    scaleD = NULL;
    nInBlock = NULL;
    indBlock = new int*[0];
    diagI = new int*[0];
    covpar = new double*[0];
    chcovpar = new double*[0];
    weightUnif = NULL;
    halfRangeUnif = NULL;
  }
  else{
    
    // Integer components
    // ==================
    _nCluster = parmI[1];    
    if (_nCluster <= 0)
      throw returnR("C++ Error: Number of clusters must be positive in bblocks constructor", 1);

    priorD = parmI + 2; 
    if (*priorD != InvWishart && _nRandom > 1)
      throw returnR("C++ Error: Prior for D must be inverse-wishart when dim(b) > 1", 1);    

    typeUpd = parmI + 3;
 
    _nBlocks = parmI[4];
 
    nInBlock = parmI + 5;
    int* cumnInBlock = new int[_nBlocks + 1];
    if (cumnInBlock == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    cumnInBlock[0] = 0;
    for (i = 0; i < _nBlocks; i++){
      if (nInBlock[i] <= 0) throw returnR("C++ Error: Incorrect nInBlock parameter supplied", 1); 
      cumnInBlock[i + 1] = cumnInBlock[i] + nInBlock[i];
    }
    if (cumnInBlock[_nBlocks] != _nRandom) throw returnR("C++ Error: Incorrect nInBlock parameter supplied", 1);  

    int lcovparLV = parmI[5 + _nBlocks];      // length of the array with covariance matrices for proposals

    int* indBlockLV = parmI + 5 + _nBlocks + 1;
    for (i = 0; i < _nRandom; i++){
      if (indBlockLV[i] < 0 || indBlockLV[i] >= _nRandom)
        throw returnR("C++ Error: Incorrect indBlockLV parameter supplied", 1);
    }
    indBlock = new int*[_nBlocks];
    if (indBlock == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < _nBlocks; i++){
      indBlock[i] = indBlockLV + cumnInBlock[i];
    }
    delete [] cumnInBlock;

    diagI = new int*[_nBlocks];
    if (diagI == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < _nBlocks; i++){
      diagI[i] = new int[nInBlock[i]];
      if (diagI[i] == NULL)
        throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
      for (j = 0; j < nInBlock[i]; j++)
        diagI[i][j] = (j * (2*nInBlock[i] - j + 1)) / 2;
    }

    sumAccept = sumAcceptP;
    for (i = 0; i < _nBlocks*_nCluster; i++) sumAccept[i] = 0;  

    // Double components
    // =================
    bM = bbM;

    switch (*priorD){
    case InvWishart:
      dfD = parmD + 0;
      scaleD = parmD + 1;      
      break;
    case SDUniform:
      dfD = parmD + 0;
      scaleD = parmD + 1;
      if (*scaleD <= 0) throw returnR("C++ Error: Upper limit of the uniform prior for std. dev.(b) <= 0.", 1);
      *dfD = 1/((*scaleD)*(*scaleD));
      break;
    }

    covpar = new double*[_nBlocks];
    chcovpar = new double*[_nBlocks];
    if (covpar == NULL || chcovpar == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    double* covparLV = parmD + 1 + (_nRandom * (1 + _nRandom))/2;
    int cumdimi = 0;
    for (i = 0; i < _nBlocks; i++){
      covpar[i] = covparLV + cumdimi;
      dimi = (nInBlock[i] * (nInBlock[i] + 1))/2;
      cumdimi += dimi;
      chcovpar[i] = new double[dimi];
      if (chcovpar[i] == NULL)
        throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);

      switch (*typeUpd){
        case RandomWalk:                // standard MH step 
          for (j = 0; j < dimi; j++) chcovpar[i][j] = covpar[i][j];
          cholesky(chcovpar[i], &helpi, nInBlock + i, diagI[i], tolerChol);
          if (helpi < nInBlock[i]) 
            throw returnR("C++ Error: The covariance matrix for MH proposal is not positive definite", 1);
          break;
   
        case Gibbs:                     // Gibbs step, covpar and chcovpar are just working spaces
          for (j = 0; j < dimi; j++) covpar[i][j] = 0.0;
          for (j = 0; j < dimi; j++) chcovpar[i][j] = 0.0;
          break;

        default:
          throw returnR("C++ Error: Unknown or unimplemented update type for random effects.", 1);
      }
    }

    halfRangeUnif = parmD + 1 + (_nRandom * (1 + _nRandom))/2 + lcovparLV;
    for (i = 0; i < _nRandom; i++)
      if (halfRangeUnif[i] < 0) 
        throw returnR("C++ Error: Invalid 'halfRangeUnif", 1);

    weightUnif = parmD  + 1 + (_nRandom * (1 + _nRandom))/2 + lcovparLV + _nRandom;
    for (i = 0; i < _nBlocks; i++){
      if (weightUnif[i] < 0) weightUnif[i] = 0;
      if (weightUnif[i] > 1) weightUnif[i] = 1;
    }    
  }
}


// Destructor
// ==========
bblocks::~bblocks()
{
  int i;

  delete [] indBlock;
  for (i = 0; i < _nBlocks; i++) delete [] diagI[i];
  delete [] diagI;
  delete [] covpar;
  for (i = 0; i < _nBlocks; i++) delete [] chcovpar[i];
  delete [] chcovpar;
}


// Print the whole object on a screen
// ================================== 
void
bblocks::print() const
{
  int i, dimi;

  cout << "nBlocks = " << _nBlocks;
  cout << ",   nRandom = " << _nRandom;
  cout << ",   nCluster = " << _nCluster;

  if (_nRandom > 0){
    switch (*priorD){
    case InvWishart:
      cout << ",   priorD = Inverse-Wishart";
      break;
    case SDUniform:
      cout << ",   priorD = SD-Uniform";
      break;
    default:
      cout << ",   priorD = ERROR";
    }

    switch (*typeUpd){
    case RandomWalk:
      cout << ",   typeUpd = Random Walk" << endl;
      break;
    case AdaptiveM:
      cout << ",   typeUpd = Adaptive Metropolis" << endl;
      break;
    case Gibbs:
      cout << ",   typeUpd = Gibbs" << endl;      
      break;
    default:
      cout << ",   priorD = ERROR" << endl;
    }

    i = _nRandom * _nCluster;
    cout << "bM = "; printArray(bM, &i);
    cout << "halfRangeUnif = "; printArray(halfRangeUnif, &_nRandom);

    cout << "nInBlock = "; printArray(nInBlock, &_nBlocks);
    for (i = 0; i < _nBlocks; i++){
      dimi = (nInBlock[i] * (nInBlock[i] + 1))/2;
      cout << "Block " << i << ":";
      cout << "  indBlock = "; printArray(indBlock[i], nInBlock + i);
      cout << "          diagI = "; printArray(diagI[i], nInBlock + i);    
      cout << "          covpar = "; printArray(covpar[i], &dimi);    
      cout << "          chcovpar = "; printArray(chcovpar[i], &dimi);        
    }

    cout << "weightUnif = "; printArray(weightUnif, &_nBlocks);
    i = _nBlocks * _nCluster;
    cout << "sumAccept = "; printArray(sumAccept, &i);
    cout << endl;
  }

  return;
}
