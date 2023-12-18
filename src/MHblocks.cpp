// Implementation of the class MHblocks.
//
//  MHblocks: a class to store information concerning the parameters that
//            will be possibly updated in blocks using a Metropolis-Hastings step
//
//   * this class is primarily designed for beta parameters
//     (i.e. regression parameters and means of random effects, except the random intercept)

// 28/01/2004: start working on it
// 19/06/2012: std::cout replaced by Rprintf()
//
#include "MHblocks.h"

using namespace std;

// Non-parametric constructor
// ===========================
MHblocks::MHblocks()
  : _nBlocks(0), _nParams(0), _maxnInBlock(0), _isdprior(0)
{
  par = NULL;
  proppar = new double[0];
  priorInvVar = new double[0];
  meanpar = NULL;
  halfRangeUnif = NULL;
  priorMean = NULL;
  priorSD = NULL;
  typeUpd = NULL;
  nInBlock = NULL;
  nRandomB = new int[0];
  nFixedB = new int[0];
  indBlock = new int*[0];
  diagI = new int*[0];
  covpar = new double*[0];
  chcovpar = new double*[0];
  logdprior = new double[0];
  weightUnif = NULL;
  eps = NULL;
  sdNum = new double[0];
  sumAccept = NULL;
}


// Copy constructor
// =================
MHblocks::MHblocks(const MHblocks& mh)
  : _nBlocks(mh._nBlocks), _nParams(mh._nParams), _maxnInBlock(mh._maxnInBlock), _isdprior(mh._isdprior)
{
  int i, j, dimi;

  par = mh.par;

  proppar = new double[_nParams];
  if (proppar == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nParams; i++) proppar[i] = mh.proppar[i];

  meanpar = mh.meanpar;

  halfRangeUnif = mh.halfRangeUnif;

  priorMean = mh.priorMean;
  priorSD = mh.priorSD;
  priorInvVar = new double[_nParams];
  if (priorInvVar == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nParams; i++) priorInvVar[i] = mh.priorInvVar[i];

  typeUpd = mh.typeUpd;

  nInBlock = mh.nInBlock;

  nRandomB = new int[_nBlocks];
  nFixedB = new int[_nBlocks];
  if (nRandomB == NULL || nFixedB == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    nRandomB[i] = mh.nRandomB[i];
    nFixedB[i] = mh.nFixedB[i];
  }

  indBlock = new int*[_nBlocks];
  if (indBlock == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++) indBlock[i] = mh.indBlock[i];

  diagI = new int*[_nBlocks];
  if (diagI == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    diagI[i] = new int[nInBlock[i]];
    if (diagI[i] == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (j = 0; j < nInBlock[i]; j++) diagI[i][j] = mh.diagI[i][j];
  }


  covpar = new double*[_nBlocks];
  chcovpar = new double*[_nBlocks];
  if (covpar == NULL || chcovpar == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    covpar[i] = mh.covpar[i];
    dimi = (nInBlock[i] * (nInBlock[i] + 1))/2;
    chcovpar[i] = new double[dimi];
    if (chcovpar[i] == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (j = 0; j < dimi; j++) chcovpar[i][j] = mh.chcovpar[i][j];
  }

  if (_isdprior){
    logdprior = new double[_nBlocks];
    if (logdprior == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < _nBlocks; i++) logdprior[i] = mh.logdprior[i];
  }
  else{
    logdprior = new double[0];
  }

  weightUnif = mh.weightUnif;

  eps = mh.eps;

  sdNum = new double[_nBlocks];
  if (sdNum == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++) sdNum[i] = mh.sdNum[i];
  
  sumAccept = mh.sumAccept;
}

// Assignment operator
// ===================
MHblocks&
MHblocks::operator=(const MHblocks& mh)
{
  int i, j, dimi;

  delete [] proppar; 
  delete [] priorInvVar;
  delete [] nRandomB;
  delete [] nFixedB;
  delete [] indBlock;
  for (i = 0; i < _nBlocks; i++) delete [] diagI[i];
  delete [] diagI;
  delete [] covpar;
  for (i = 0; i < _nBlocks; i++) delete [] chcovpar[i];
  delete [] chcovpar;
  delete [] logdprior;
  delete [] sdNum;

  _nBlocks = mh._nBlocks; 
  _nParams = mh._nParams;
  _maxnInBlock = mh._maxnInBlock;
  _isdprior = mh._isdprior;

  par = mh.par;

  proppar = new double[_nParams];
  if (proppar == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nParams; i++) proppar[i] = mh.proppar[i];

  meanpar = mh.meanpar;

  halfRangeUnif = mh.halfRangeUnif;

  priorMean = mh.priorMean;
  priorSD = mh.priorSD;
  priorInvVar = new double[_nParams];
  if (priorInvVar == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nParams; i++) priorInvVar[i] = mh.priorInvVar[i];

  typeUpd = mh.typeUpd;

  nInBlock = mh.nInBlock;

  nRandomB = new int[_nBlocks];
  nFixedB = new int[_nBlocks];
  if (nRandomB == NULL || nFixedB == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    nRandomB[i] = mh.nRandomB[i];
    nFixedB[i] = mh.nFixedB[i];
  }

  indBlock = new int*[_nBlocks];
  if (indBlock == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++) indBlock[i] = mh.indBlock[i];

  diagI = new int*[_nBlocks];
  if (diagI == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    diagI[i] = new int[nInBlock[i]];
    if (diagI[i] == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (j = 0; j < nInBlock[i]; j++) diagI[i][j] = mh.diagI[i][j];
  }


  covpar = new double*[_nBlocks];
  chcovpar = new double*[_nBlocks];
  if (covpar == NULL || chcovpar == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++){
    covpar[i] = mh.covpar[i];
    dimi = (nInBlock[i] * (nInBlock[i] + 1))/2;
    chcovpar[i] = new double[dimi];
    if (chcovpar[i] == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (j = 0; j < dimi; j++) chcovpar[i][j] = mh.chcovpar[i][j];
  }

  if (_isdprior){
    logdprior = new double[_nBlocks];
    if (logdprior == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < _nBlocks; i++) logdprior[i] = mh.logdprior[i];
  }
  else{
    logdprior = new double[0];
  }

  weightUnif = mh.weightUnif;

  eps = mh.eps;

  sdNum = new double[_nBlocks];
  if (sdNum == NULL)
    throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
  for (i = 0; i < _nBlocks; i++) sdNum[i] = mh.sdNum[i];
  
  sumAccept = mh.sumAccept;

  return *this;
}


// Parametric constructor
// ======================
MHblocks::MHblocks(double* parP,             int* parmI,        double* parmD,                   
                   int* sumAcceptP,          const int* indbA,
                   const double* tolerChol,  const int* iter,   const int* wantdprior)
  : _isdprior(*wantdprior)
{
  // All priors are assumed to be normal. Include here this part of code to allow for possible future changes
  // ========================================================================================================
  double (*dprior) (double, double, double, int);
  dprior = dnorm;

  int i, j, helpi, dimi;

  _nBlocks = parmI[0];

  if (_nBlocks == 0){
    _nParams = 0;
    _maxnInBlock = 0;
    _isdprior = 0;

    par = NULL;
    proppar = new double[0];
    meanpar = NULL;
    halfRangeUnif = NULL;
    priorMean = NULL;
    priorSD = NULL;
    priorInvVar = new double[0];
    typeUpd = NULL;
    nInBlock = NULL;
    nRandomB = new int[0];
    nFixedB = new int[0];
    indBlock = new int*[0];
    diagI = new int*[0];
    covpar = new double*[0];
    chcovpar = new double*[0];
    logdprior = new double[0];
    weightUnif = NULL;
    eps = NULL;
    sdNum = new double[0];
    sumAccept = NULL;
  }
  else{

    // Integer components
    // ===================
    _nParams = parmI[1];

    nInBlock = parmI + 2;
    int* cumnInBlock = new int[_nBlocks + 1];
    if (cumnInBlock == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    cumnInBlock[0] = 0;
    for (i = 0; i < _nBlocks; i++){
      if (nInBlock[i] <= 0) throw returnR("C++ Error: Incorrect nInBlock parameter supplied", 1); 
      cumnInBlock[i + 1] = cumnInBlock[i] + nInBlock[i];
    }
    if (cumnInBlock[_nBlocks] != _nParams) throw returnR("C++ Error: Incorrect nInBlock parameter supplied", 1);  

    _maxnInBlock = parmI[2 + _nBlocks];

    int lcovparLV = parmI[2 + _nBlocks + 1];      // length of the array with covariance matrices

    int* indBlockLV = parmI + 2 + _nBlocks + 2;
    for (i = 0; i < _nParams; i++){
      if (indBlockLV[i] < 0 || indBlockLV[i] >= _nParams)
        throw returnR("C++ Error: Incorrect indBlockLV parameter supplied", 1);
    }
    indBlock = new int*[_nBlocks];
    if (indBlock == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < _nBlocks; i++){
      indBlock[i] = indBlockLV + cumnInBlock[i];      
    }

    nRandomB = new int[_nBlocks];
    nFixedB = new int[_nBlocks];
    if (nRandomB == NULL || nFixedB == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < _nBlocks; i++){
      nRandomB[i] = 0;
      nFixedB[i] = 0;
      for (j = 0; j < nInBlock[i]; j++){
        if (indbA[indBlock[i][j]] < 0) nFixedB[i]++;
        else                           nRandomB[i]++;
      }      
    }

    typeUpd = parmI + 2 + _nBlocks + 2 + _nParams;

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
    for (i = 0; i < _nBlocks; i++) sumAccept[i] = 0;  


    // Double components
    // ===================
    par = parP;
    proppar = new double[_nParams];
    if (proppar == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < _nParams; i++) proppar[i] = par[i];

    priorMean = parmD;
    priorSD = parmD + _nParams;
    priorInvVar = new double[_nParams];
    if (priorInvVar == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < _nParams; i++){
      if (priorSD[i] <= 0) throw returnR("C++ Error: Negative and zero prior variances are not allowed.", 1);
      priorInvVar[i] = 1/priorSD[i];
      priorSD[i] = sqrt(priorSD[i]);
    }

    meanpar = parmD + 2*_nParams;
    if (*iter == 0)
      for (i = 0; i < _nParams; i++) meanpar[i] = par[i];

    halfRangeUnif = parmD + 3*_nParams;
    for (i = 0; i < _nParams; i++)
      if (halfRangeUnif[i] < 0) 
        throw returnR("C++ Error:  Invalid 'halfRangeUnif", 1);

    covpar = new double*[_nBlocks];
    chcovpar = new double*[_nBlocks];
    if (covpar == NULL || chcovpar == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    double* covparLV = parmD + 4*_nParams;
    int cumdimi = 0;
    for (i = 0; i < _nBlocks; i++){
      covpar[i] = covparLV + cumdimi;
      dimi = (nInBlock[i] * (nInBlock[i] + 1))/2;
      cumdimi += dimi;
      chcovpar[i] = new double[dimi];
      if (chcovpar[i] == NULL)
        throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);

      // Cholesky decomposition of the covariance matrix for the proposal
      switch (typeUpd[i]){
        case RandomWalk:                // standard MH step 
    	  for (j = 0; j < dimi; j++) chcovpar[i][j] = covpar[i][j];        
          cholesky(chcovpar[i], &helpi, nInBlock + i, diagI[i], tolerChol);
          if (helpi < nInBlock[i]) 
            throw returnR("C++ Error: The covariance matrix for MH proposal is not positive definite.", 1);
          break;
 
        case AdaptiveM:               // adaptive Metropolis step
	  for (j = 0; j < dimi; j++) chcovpar[i][j] = covpar[i][j];
          cholesky(chcovpar[i], &helpi, nInBlock + i, diagI[i], tolerChol);
          if (helpi < nInBlock[i]) 
            throw returnR("C++ Error: The initial covariance matrix for AM proposal is not positive definite.", 1);
          break;

        case Gibbs:                     // Gibbs step, covpar and chcovpar are just working spaces
          for (j = 0; j < dimi; j++) covpar[i][j] = 0.0;
          for (j = 0; j < dimi; j++) chcovpar[i][j] = 0.0;
          if (nRandomB[i] > 0 && nFixedB[i] > 0)
            throw returnR("C++ Error: Fixed effects and means of random effects may not be updated in one block when Gibbs used.", 1);
          break;

        default:
          throw returnR("C++ Error: Unknown or unimplemented update type for fixed effects.", 1);
      }
    }

    if (_isdprior){
      logdprior = new double[_nBlocks];
      if (logdprior == NULL)
        throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
      for (i = 0; i < _nBlocks; i++){
        logdprior[i] = 0.0;
        for (j = 0; j < nInBlock[i]; j++){
          logdprior[i] += dprior(par[indBlock[i][j]], priorMean[indBlock[i][j]], priorSD[indBlock[i][j]], 1);
        }
      }
    }
    else{
      logdprior = new double[0];
    }


    weightUnif = parmD + 4*_nParams + lcovparLV;
    for (i = 0; i < _nBlocks; i++){
      if (weightUnif[i] < 0) weightUnif[i] = 0;
      if (weightUnif[i] > 1) weightUnif[i] = 1;
    }    
  
    eps = parmD + 4*_nParams + lcovparLV + _nBlocks;
    for (i = 0; i < _nBlocks; i++){
      if (eps[i] < 0) throw returnR("Incorrect epsilon for (AM) algorithm supplied", 1);
    }

    double* sdNumP = parmD + 4*_nParams + lcovparLV + 2*_nBlocks; 
    sdNum = new double[_nBlocks];
    if (sdNum == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < _nBlocks; i++){
      sdNum[i] = sdNumP[nInBlock[i] - 1];
      if (sdNum[i] <= 0) throw returnR("C++ Error: Incorrect s(d) numbers for (AM) algorithm supplied", 1);
    }    

    delete [] cumnInBlock;
  }

}


// Destructor
// ==========
MHblocks::~MHblocks()
{
  int i;

  delete [] proppar; 
  delete [] priorInvVar;
  delete [] nRandomB;
  delete [] nFixedB;
  delete [] indBlock;
  for (i = 0; i < _nBlocks; i++) delete [] diagI[i];
  delete [] diagI;
  delete [] covpar;
  for (i = 0; i < _nBlocks; i++) delete [] chcovpar[i];
  delete [] chcovpar;
  delete [] logdprior;
  delete [] sdNum;
}


// Print the whole object on a screen
// ================================== 
void
MHblocks::print() const
{
  int i, dimi;

  //std::cout << "nBlocks = " << _nBlocks;
  //std::cout << ",   nParams = " << _nParams;
  //std::cout << ",   maxnInBlock = " << _maxnInBlock << std::endl;
  //std::cout << ",   isdprior = " << _isdprior << std::endl;
  Rprintf("nBlocks = %d,   nParams = %d,   maxnInBlock = %d\n,   isdprior = %d\n", _nBlocks, _nParams, _maxnInBlock, _isdprior);

  if (_nBlocks > 0){
    //std::cout << "nFixedB = "; 
    Rprintf("nFixedB = ");
    printArrayI(nFixedB, &_nBlocks);
    //std::cout << "nRandomB = "; 
    Rprintf("nRandomB = ");
    printArrayI(nRandomB, &_nBlocks);    

    //std::cout << "par = "; 
    Rprintf("par = ");
    printArrayD(par, &_nParams);
    //std::cout << "proppar = "; 
    Rprintf("proppar = ");
    printArrayD(proppar, &_nParams);
    //std::cout << "meanpar = "; 
    Rprintf("meanpar = ");
    printArrayD(meanpar, &_nParams);
    //std::cout << "halfRangeUnif = "; 
    Rprintf("halfRangeUnif = ");
    printArrayD(halfRangeUnif, &_nParams);

    //std::cout << "priorMean = "; 
    Rprintf("priorMean = ");
    printArrayD(priorMean, &_nParams);
    //std::cout << "priorSD = "; 
    Rprintf("priorSD = ");
    printArrayD(priorSD, &_nParams);
    //std::cout << "priorInvVar = "; 
    Rprintf("priorInvVar = ");
    printArrayD(priorInvVar, &_nParams);

    //std::cout << "logdprior = "; 
    Rprintf("logdprior = ");
    printArrayD(logdprior, &_nBlocks);
    //std::cout << "typeUpd = ";  
    Rprintf("typeUpd = ");
    printArrayI(typeUpd, &_nBlocks);
    //std::cout << "nInBlock = "; 
    Rprintf("nInBlock = ");
    printArrayI(nInBlock, &_nBlocks);
    for (i = 0; i < _nBlocks; i++){
      dimi = (nInBlock[i] * (nInBlock[i] + 1))/2;
      //std::cout << "Block " << i << ":";
      //std::cout << "  indBlock = "; 
      Rprintf("Block %d:  indBlock = ", i);
      printArrayI(indBlock[i], nInBlock + i);
      //std::cout << "          diagI = "; 
      Rprintf("          diagI = ");
      printArrayI(diagI[i], nInBlock + i);    
      //std::cout << "          covpar = "; 
      Rprintf("          covpar = ");
      printArrayD(covpar[i], &dimi);    
      //std::cout << "          chcovpar = "; 
      Rprintf("          chcovpar = ");
      printArrayD(chcovpar[i], &dimi);        
    }

    //std::cout << "weightUnif = "; 
    Rprintf("weightUnif = ");
    printArrayD(weightUnif, &_nBlocks);
    //std::cout << "eps = "; 
    Rprintf("eps = ");
    printArrayD(eps, &_nBlocks);
    //std::cout << "sdNum = "; 
    Rprintf("sdNum = ");
    printArrayD(sdNum, &_nBlocks);
    //std::cout << "sumAccept = "; 
    Rprintf("sumAccept = ");
    printArrayI(sumAccept, &_nBlocks);
    //std::cout << std::endl;
    Rprintf("\n");
  }

  return;
}
