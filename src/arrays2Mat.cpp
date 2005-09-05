// Function to create covariance matrices, response vectors etc.
//   obtained as pointers.

// 03/11/2003: start working on it
// 19/03/2004: createData splitted into two parts such that the shorter version 
//             can be called by predictive functions
// 23/03/2004: sub-parts of predictive added

#include "arrays2Mat.h"

using namespace std;

// ======================================================================================================
// create2pointer: Allocate space for a second layer of double** and initialize it to 'init' (if wanted)
// ======================================================================================================
//
// PARAMETERS:
//
// point ...... array of arrays
// lpoint ..... length of point
// llayer2 .... desired length of each array that is stored in point
// init ....... initial value
//
void
create2pointer(double**& point,    const int lpoint,  const int llayer2,  
               const double init,  const bool doinit)
{
  int i, j;  

  for (i = 0; i <= lpoint; i++){
    point[i] = new double[llayer2];
    if (point[i] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
    if (doinit) for (j = 0; j < llayer2; j++) point[i][j] = init;
  }
  return;
}


// ======================================================================================================
// create3pointer: PROTOTYPE 1
//                 Allocate space for the second and third layer of double*** and initialize it to 'init'
//                 (if wanted))
//                 * the SECOND layer may have variable length
// ======================================================================================================
//
// PARAMETERS:
//
// point ............... array of arrays of arrays
// lpoint .............. length of point
// llayer2[lpoint] ..... desired length of each array  of arrays that is stored in point
// llayer3 ............. desired length of each array in the 3rd layer
// init ................ initial value
//
void
create3pointer(double***& point,   const int lpoint,  const int* llayer2,  const int llayer3,  
               const double init,  const bool doinit)
{
  int obs, grid, j;

  for (obs = 0; obs <= lpoint; obs++){
    point[obs] = new double*[llayer2[obs]];
    if (point[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
    for (grid = 0; grid < llayer2[obs]; grid++){
      point[obs][grid] = new double[llayer3];
      if (point[obs][grid] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
      if (doinit) for (j = 0; j < llayer3; j++) point[obs][grid][j] = init;
    }
  }

  return;
}


// ======================================================================================================
// create3pointer: PROTOTYPE 2
//                 Allocate space for the second and third layer of double*** and initialize it to 'init'
//                 (if wanted))
//                 * the THIRD layer may have variable length
// ======================================================================================================
//
// PARAMETERS:
//
// point ............... array of arrays of arrays
// lpoint .............. length of point
// llayer2 ............. desired length of each array  of arrays that is stored in point
// llayer3[llayer2] .... desired length of each array in the 3rd layer
// init ................ initial value
//
void
create3pointer(double***& point,   const int lpoint,  const int llayer2,  const int* llayer3,  
               const double init,  const bool doinit)
{
  int i, obs, j;

  for (i = 0; i <= lpoint; i++){
    point[i] = new double*[llayer2];
    if (point[i] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
    for (obs = 0; obs < llayer2; obs++){
      point[i][obs] = new double[llayer3[obs]];
      if (point[i][obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
      if (doinit) for (j = 0; j < llayer3[obs]; j++) point[i][obs][j] = init;
    }
  }

  return;
}


// ==================================================================================================
// ******* createDataShort: subset of overall createData
// ==================================================================================================
//
// PARAMETERS:
//   read the list of parameters in 'createData' function
//
void
createDataShort(int* nwithinA,          int* clusteriA,            List<int>* invclusteriA,
                const double* XA,
                double** ZZt,           int* diagIZZt,             int* indbinXA,     
                const int* nP,          const int* nclusterP,
                const int* nXP,         const int* nfixedP,        const int* nrandomP,
                const int* randomIntP,  const int* indbA)
{
  int i, j, k;

  // Numbers of observations in clusters (assign all to one or check for consistency)
  // ================================================================================
  if (*nP == *nclusterP) 
    for (i = 0; i < *nclusterP; i++) nwithinA[i] = 1;
  else{
    j = 0;
    for (i = 0; i < *nclusterP; i++) j += nwithinA[i];
    if (j != *nP){
      returnR error("C++ Error: Incorrect number of observations supplied.", 99);
      throw error; 
    }
  }


  // Indicators of clusters for observations
  // =======================================
  if (*nP == *nclusterP){
    for (i = 0; i < *nP; i++){
      clusteriA[i] = i;
    }
  }
  else{
    k = 0;
    for (i = 0; i < *nclusterP; i++){
      for (j = 0; j < nwithinA[i]; j++){
        clusteriA[k] = i;
        k++;
      }          
    }
  }


  // Inverse-indicators of clusters for observations
  // =======================================
  if (*nP == *nclusterP){
    for (i = 0; i < *nP; i++){
      invclusteriA[i].addNode(i);
    }
  }
  else{
    k = 0;
    for (i = 0; i < *nclusterP; i++){
      for (j = 0; j < nwithinA[i]; j++){
        invclusteriA[i].addNode(k);
        k++;
      }          
    }
  }


  // Check indb, nrandom and nfixed for consistencies
  // ================================================
  if ((*nfixedP + (*nrandomP) - (*randomIntP)) != (*nXP)){
    returnR error("C++ Error: Incorrect dimensions supplied.", 99);
    throw error; 
  }  

  for (i = 0; i < *nXP; i++){
    if (!(indbA[i] == -1 || (indbA[i] >= *randomIntP && indbA[i] < *nrandomP))){
      returnR error("C++ Error: Incorrect indb vector supplied.", 99);
      throw error; 
    }
  }


  // Indeces of columns in X corresponding to random effects
  // ========================================================
  if (*nrandomP){
    if (*randomIntP) indbinXA[0] = -1;
    for (i = 0; i < *nXP; i++) if (indbA[i] != -1) indbinXA[indbA[i]] = i;
  }


  // z*z' for random effects, random intercept included
  // ====================================================
  for (i = 0; i < *nrandomP; i++) diagIZZt[i] = (i * (2*(*nrandomP) - i + 1)) / 2;

  if (*nrandomP && (*nP)){
    if (*randomIntP){
      for (i = 0; i < *nP; i++){
        ZZt[i][0] = 1.0;
        for (k = 1; k < *nrandomP; k++){
          ZZt[i][k] = XA[indbinXA[k]*(*nP) + i];
        }
        for (j = 1; j < *nrandomP; j++){
          for (k = j; k < *nrandomP; k++){      
            ZZt[i][diagIZZt[j] + k - j] = XA[indbinXA[j]*(*nP) + i] * XA[indbinXA[k]*(*nP) + i];
          }
        }
      }
    }
    else{
      for (i = 0; i < *nP; i++){
        for (j = 0; j < *nrandomP; j++){
          for (k = j; k < *nrandomP; k++){
            ZZt[i][diagIZZt[j] + k - j] = XA[indbinXA[j]*(*nP) + i] * XA[indbinXA[k]*(*nP) + i];
          }
        }
      }
    }
  }

  return;
}


// ==================================================================================================
// ****** createData: create matrices with data (response, covariance matrices, status vector) ******
// ==================================================================================================
//
// INPUT PARAMETERS:
//
// XA ............. design matrix                       
// YA ............. response matrix array
// nP ............. number of observations
// nclusterP ...... number of clusters
// nYP ............ number of columns in the response matrix (2 or 3)
// nXP ............ number of columns in the design matrix
// nfixedP ........
// nrandomP .......
// randomIntP .....
// indbA ..........
//
// OUTPUT PARAMETERS:
//
// nwithinA ....... number of observations within clusters
// clusteriA ...... indeces of clusters for observations
// statusA ........ censoring status
// Y1A ............ lower limit of the response
// Y2A ............ upper limit of the response
// ZZt ............ z*z' for all random effects (lower triangle only)
// diagIZZt ....... indeces of diagonal elemets of above mentioned matrices
// indbinXA .......    
//
// indBlockbeta ... it is not 'const' but it is not changed by this function
//
void 
createData(int* nwithinA,          int* clusteriA,            List<int>* invclusteriA, 
           int* statusA,           double*& Y1A,              double*& Y2A,           
           double** ZZt,           int* diagIZZt,             int* indbinXA,     
           double*** XXt,          int** diagIXXt,      
           const double* XA,       double* YA,
           const int* nP,          const int* nclusterP,      const int* nYP,
           const int* nXP,         const int* nfixedP,        const int* nrandomP,
           const int* randomIntP,  const int* indbA,      
           const int& nBlockbeta,  const int* nInBlockbeta,   int** indBlockbeta,        const int* typeUpdbeta)
{
  int i, j, k, l, block;

  createDataShort(nwithinA, clusteriA, invclusteriA, XA, ZZt, diagIZZt, indbinXA, nP, nclusterP, nXP, nfixedP, nrandomP, randomIntP, indbA);


  // Status vector and response vector/matrix
  // =========================================
    // statusA ..... vector (n x 1) with status 
    //               (0 = right censored, 1 = event, 2 = left censored, 3 = interval censored)
    // Y1A ......... vector (n x 1) with response (exact or right censored 
    //               or lower limit of int. cens.)
    // Y2A ......... vector (n x 1) with response (upper limit of interval censored)
  if (*nYP == 2){                                    // only exact/right censored observations
     for (i = 0; i < *nP; i++) statusA[i] = int(YA[*nP + i]);
     Y1A = YA;
     Y2A = NULL;
  }
  else{             // *nyP is 3                     // also interval censored observations
    for(i = 0; i < *nP; i++) statusA[i] = int(YA[2*(*nP) + i]);
    Y1A = YA;
    Y2A = YA + (*nP);
  }


  // x*x' matrices for fixed effects
  // =================================
  for (block = 0; block < nBlockbeta; block++){
    if (typeUpdbeta[block] == Gibbs){
      for (i = 0; i < nInBlockbeta[block]; i++) diagIXXt[block][i] = (i * (2*(nInBlockbeta[block]) - i + 1)) / 2;
      for (j = 0; j < *nP; j++){
        for (k = 0; k < nInBlockbeta[block]; k++){
          for (l = k; l < nInBlockbeta[block]; l++){
            XXt[block][j][diagIXXt[block][k] + l - k] = XA[indBlockbeta[block][k]*(*nP) + j] * XA[indBlockbeta[block][l]*(*nP) + j];
          }
        }
      }
    }
  }    
  

  return;
}   // end of function createData


// ==================================================================================================
// ****** createParam: create matrices with initial values of the parameters
// ==================================================================================================
//
// INPUT PARAMETERS:
//
// nP ............... number of observations
// kmaxP ............ maximal number of mixture components
// mixtureA ......... mixture parameters in a condensed format
// rM ............... observation allocations (in R indeces, i.e.  1, 2, ...)
//
// OUTPUT PARAMETERS:
//
// wM ............... mixture weights
// muM .............. mixture means             
// invsigma2M ....... mixture inverse-variances 
// rM ............... observation allocations (in C++ indeces, i.e. 0, 1, ...)
// invrM ............ indeces of observations belonging to the mixture components 
// mixtureNM ........ numbers of observations belonging to each mixture component
//
// propwM
// propmuM
// propinvsigma2M ... same as above (used to store proposed values when split-combine or birth death moves are carried out)
// proprM
// propinvrM   
// propmixtureNM
//
void 
createParam(const int* nP,                 const int* kmaxP,             const double* mixtureA,        
            double* wM,                    double* muM,                  double* invsigma2M,
            int* rM,                       List<int>* invrM,             int* mixtureNM,  
            double* propwM,                double* propmuM,              double* propinvsigma2M,
            int* proprM,                   List<int>* propinvrM,         int* propmixtureNM 
           )
{
  int j, k, obs;


  // Mixture parameters
  // ==================
  k = int(mixtureA[0]);
  for (j = 0; j < k; j++){
    wM[j] = mixtureA[1+j];                             
    propwM[j] = wM[j];
    muM[j] = mixtureA[1 + *kmaxP + j];                 
    propmuM[j] = muM[j];
    invsigma2M[j] = 1/mixtureA[1 + 2*(*kmaxP) + j];    
    propinvsigma2M[j] = invsigma2M[j];
  }
  for (j = k; j < *kmaxP; j++){
    wM[j] = 0.0;
    propwM[j] = 0.0; 
    muM[j] = 0.0;
    propmuM[j] = 0.0;
    invsigma2M[j] = 0.0;
    propinvsigma2M[j] = 0.0;
  }


  // Component pertinence
  // ====================
  for (obs = 0; obs < *nP; obs++){
    rM[obs] -= 1;                       // R index --> C++ index
    proprM[obs] = rM[obs];
  }

  
  // Inverse component pertinence
  // ============================
  for (j = 0; j < *kmaxP; j++){      
    invrM[j] = List<int>();
    propinvrM[j] = List<int>();
  }
  for (obs = 0; obs < *nP; obs++){
    invrM[rM[obs]].addNode(obs);
    propinvrM[rM[obs]].addNode(obs);
  }

  // Numbers of observations belonging to each mixture component
  // ===========================================================
  giveMixtureN(mixtureNM, kmaxP, invrM);
  giveMixtureN(propmixtureNM, kmaxP, invrM);
  
  return;

}   // end of function createParam


// ==================================================================================================
// ****** createPriors: create matrices with initial values of some prior parameters
// ==================================================================================================
//
// INPUT PARAMETERS:
//
// kmaxP ............
// priorParD ........
//
// OUTPUT PARAMETERS:
//
// piSplitM ......... probabilities of the split move depending on k	  
// logpiSplitM ...... log-probabilities of the split move depending on k  
// logpiCombineM .... log-probabilities of the combine move depending on k
// piBirthM ......... probabilities of the birth move depending on k	  
// logpiBirthM ...... log-probabilities of the birth move depending on k  
// logpiDeathM ....   log-probabilities of the death move depending on k
//
void
createPriors(const int *kmaxP,     const double* priorParD,
             double* piSplitM,       
             double* logpiSplitM,  double* logpiCombineM,
             double* piBirthM,       
             double* logpiBirthM,  double* logpiDeathM
            ) 
{
  int j;

  piSplitM[0] = 0.0;
  logpiSplitM[0] = 0.0;
  logpiCombineM[0] = 0.0;

  piBirthM[0] = 0.0;
  logpiBirthM[0] = 0.0;
  logpiDeathM[0] = 0.0;

  for (j = 1; j <= *kmaxP; j++){
     piSplitM[j] = priorParD[j - 1];
     piBirthM[j] = priorParD[*kmaxP + j - 1];

     if (piSplitM[j] <= 0.0){
       logpiSplitM[j] = -FLT_MAX;
       logpiCombineM[j] = 0.0;
     }
     else{
       if (piSplitM[j] >= 1.0){
         logpiSplitM[j] = 0.0;
         logpiCombineM[j] = -FLT_MAX;
       }
       else{
         logpiSplitM[j] = log(piSplitM[j]);
         logpiCombineM[j] = log(1 - piSplitM[j]);
       } 
     } 

     if (piBirthM[j] <= 0.0){
       logpiBirthM[j] = -FLT_MAX;
       logpiDeathM[j] = 0.0;
     }
     else{
       if (piBirthM[j] >= 1.0){
         logpiBirthM[j] = 0.0;
         logpiDeathM[j] = -FLT_MAX;
       }
       else{
         logpiBirthM[j] = log(piBirthM[j]);
         logpiDeathM[j] = log(1 - piBirthM[j]);
       } 
     } 

  }


  return;

}   // end of function createPriors

