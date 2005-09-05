// Function to compute predictive quantities based on a Bayesian survival regression
//  from 'bayessurvreg' functions
//
// * predictive survivor, hazard and cum. hazard functions are stored only when asked
// * predictive survivor times are stored always
// * for all quantities, always mean and desired quantiles are computed

//
// Remark: predictive values of random effects are alewys sampled using a Gibbs move

// 07/03/2004: start working on it
// 23/03/2004: version without computation of quantiles and means finished
// 18/08/2004: time0 != 0 allowed
//

#include "predictive.h"

extern "C" {

// PARAMETERS:
//
// errorTypeP[1] ....... type of the error distribution (currently only Mixture or Spline are allowed with this function)
// dirP ................ directory where to store files with simulated values
//

// dimsP[nclusterP + 7] ..... dimensionality parameters same as in bayessurvreg1
//         nP ...................... number of observations
//         nclusterP ............... number of clusters
//         nwithinA[nclustersP] .... numbers of observations within clusters
//         nYP ..................... number of columns in the response matrix          (ignored here)
//         nXP ..................... number of columns in the design matrix
//         nfixedP ................. number of fixed effects
//         nrandomP ................ number of random effects (including the random intercept)
//         randomIntP .............. 0/1, is a random intercept in the model or not
//         nToStore ................ (currently ignored, always changed to nwrite)
//                                   
//
// dims2P[1+nP] ............... additional dimensionality parameters
//        nquantileP[1]
//        ngridM[nP]
//

// XA[nP x nXP] ......... design matrix for both fixed and random effects
// indbA[nXP]..... ...... an array of length *nXP, indicating where to find random effects
//         indbA[j] = -1 if the jth column of X corresponds to the fixed effect
//         indbA[j] = k if the jth column of X corresponds to the kth random effect
//         i.e. values of -1, 1, ..., nrandom - 1 appear in this array if *randomIntP
//         and values of -1, 0, ..., nrandom - 1 appear in this array if !(*randomIntP)  
//

// gridA[sum_{i=0}^{nP-1} dims2P[i]] .... grid to evaluate predictive survivor, hazard and cumulative hazard

// priorParI[3]
//    kmaxP .......... maximal number of mixture components 
//    xxxxx .......... not used
//    Eb0dependmix ... 0/1

// priorParD[1]
//    time0 .......... starting time of the survival model (zero usually)

// nsimulP[3] ......... parameters of either done simulation or simulation to predict other wuantities
//         niter ..... number of done iterations (number of rows in files with sampled values)
//         nthin ..... ignored here
//         nwrite .... how often write information on number of done iterations
//

// predictP[5] ........ 0/1 indicating which quantities are to be predicted
//         predictETP
//         predictTP
//         predictSurvP
//         predicthazardP
//         predictcumhazardP

// storeP[5] ........ 0/1 indicating which functions are to be stored
//         storeETP
//         storeTP
//         storeSurvP
//         storehazardP
//         storecumhazardP

// tolersP[2] ......... tolerances for decompositions
//         tolerCholP . tolerance for the cholesky decomposition
//         tolerQRP ... tolerance for the QR decomposition

// errP[1] ....... error flag
//
void
predictive(const int* errorTypeP,  const char** dirP,        int* dimsP,          int* dims2P, 
           const double* XA,       const int* indbA,         double* quantileA,   double* gridA,
           const int* priorParI,   const double* priorParD,  const int* nsimulP,  const int* skipP,          
           const int* byP,         int* onlyAver,            int* predictP,       int* storeP,
           const double* tolersP,  int* errP)
{
  try{/******/
    double dtemp;
    int itemp;
    int obs, i, j, cl;
    double help;

    if (!(*errorTypeP == Mixture || *errorTypeP == Spline))
      throw returnR("C++ Error: unimplemented 'errorType' in a call to function 'predictive'.", 1);

    GetRNGstate();

    *errP = 0;
    const string dir = *dirP;

    // Transformation and inverse transformation of the response
    // ==========================================================
    double (*itrans)(const double);                itrans = exp;
    double (*trans)(const double);                 trans = log;

    // Starting time for the survival model
    // ====================================
    const double* time0P = priorParD + 0;

    // Dimensionality parameters
    // =========================
    int* nP = dimsP;
    int* nclusterP = dimsP + 1;
    int* nwithinA = dimsP + 2;
    //    int* nYP = dimsP + 2 + (*nclusterP);   /* nowhere used in this function */
    int* nXP = dimsP + 2 + (*nclusterP) + 1;
    int* nfixedP = dimsP + 2 + (*nclusterP) + 2;
    int* nrandomP = dimsP + 2 + (*nclusterP) + 3;
    int* randomIntP = dimsP + 2 + (*nclusterP) + 4;

    // Additional dimensionality parameters
    // =====================================
    int* nquantileP = dims2P + 0;
    const int* ngridM = dims2P + 1;
    if (*nquantileP <= 0) *onlyAver = 1;
    if (*onlyAver) *nquantileP = 0;

    // Parameters of done or future simulation
    // ========================================
    int niter = nsimulP[0]; 
    //    int nthin = nsimulP[1];              /* nowhere used in this function */
    int nwrite = nsimulP[2];

    // Tolerances for decompositions
    // =============================
    const double* tolerCholP = tolersP;
    const double* tolerQRP = tolersP + 1;

    // What to predict?
    // =================
    int* predictETP = predictP + 0;
    int* predictTP = predictP + 1;
    int* predictSurvP = predictP + 2;
    int* predicthazardP = predictP + 3;
    int* predictcumhazardP = predictP + 4;
    if (*predicthazardP || *predictcumhazardP) *predictSurvP = 1;
    const bool predictShch = (*predictSurvP) || (*predicthazardP) || (*predictcumhazardP);

    // What to store?
    // ================
    int* storeETP = storeP + 0;             if (!(*predictETP)) *storeETP = 0;
    int* storeTP = storeP + 1;              if (!(*predictTP)) *storeTP = 0;
    int* storeSurvP = storeP + 2;           if (!(*predictSurvP)) *storeSurvP = 0;
    int* storehazardP = storeP + 3;         if (!(*predicthazardP)) *storehazardP = 0;
    int* storecumhazardP = storeP + 4;      if (!(*predictcumhazardP)) *storecumhazardP = 0;
    const bool store = (*storeETP) || (*storeTP) || (*storeSurvP) || (*storehazardP) || (*storecumhazardP);

    // Indeces of quantile values in sampled chain (indexing starting from 0)
    // ======================================================================
    // nread .................... number of iterations that will be considered now (after discarding skipped and thinned values)
    // indquant1, indquant2 ..... quantile = q*sample[indquant1] + (1-q)sample[indquant2]
    // help ..................... helping number
    //
    if (*skipP >= niter) throw returnR("C++ Error: Too many iterations are to be skipped", 1);
    if (*byP <= 0) throw returnR("C++ Error: by parameter must be positive", 1);
    if (*skipP < 0) throw returnR("C++ Error: skip parameter must not be negative", 1);
    const int nread = 1 + (niter - (*skipP) - 1) / (*byP);

    int *indquant1 = &itemp;
    int *indquant2 = &itemp;
    if (!(*onlyAver)){
      indquant1  = (int*) calloc(*nquantileP, sizeof(int));
      indquant2  = (int*) calloc(*nquantileP, sizeof(int));
      if (!indquant1 || !indquant2) throw returnR("C++ Error: Could not allocate working memory.", 1);
      for (i = 0; i < *nquantileP; i++){
        if (quantileA[i] < 0 || quantileA[i] > 1) throw returnR("C++ Error: Incorrect quantile values supplied.", 1);
        if (quantileA[i] <= 0) indquant1[i] = indquant2[i] = 0;
        else if (quantileA[i] >= 1) indquant1[i] = indquant2[i] = nread - 1;
             else{
               help = quantileA[i] * double(nread);
               if (fabs(help - floor(help + 1e-8)) < 1e-8){   // help is integer
                 indquant1[i] = int(floor(help)) - 1;
                 indquant2[i] = int(floor(help));
               }
               else{
                 indquant1[i] = indquant2[i] = int(floor(help));
               }
             }
//      cout << "indquant1 = "; printArray(indquant1, nquantileP);
//      cout << "indquant2 = "; printArray(indquant2, nquantileP);
      }
    }


    // Create data matrices and variables
    // ==================================
    // *clusteriA ........... indicators of clusters for observations
    // *invclusteriA ........
    // *indbinXA
    // **ZZt ................ this is finally not needed
    // *diagIZZt ............ this is finally not needed
    //
    int* clusteriA          = (int*) calloc(*nP, sizeof(int));                          if (!clusteriA)    *errP = 1;
    List<int>* invclusteriA = new List<int>[*nclusterP];                                if (!invclusteriA) *errP = 1;
    int* indbinXA           = (int*) calloc(*nrandomP ? *nrandomP : 1, sizeof(int));    if (!indbinXA)     *errP = 1;
    double** ZZt            = (double**) calloc(*nP, sizeof(double*));                  if (!ZZt)          *errP = 1;
    int* diagIZZt           = (int*) calloc(*nrandomP, sizeof(int));                    if (!diagIZZt)     *errP = 1;
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for a working space", 1);
    for (j = 0; j < *nP && !(*errP); j++){
      ZZt[j] = (double*) calloc(*nrandomP * (*nrandomP), sizeof(double));               if (!ZZt[j]) *errP = 1;
    }
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for a working space", 1);

    createDataShort(nwithinA, clusteriA, invclusteriA, XA, ZZt, diagIZZt, indbinXA, 
                    nP, nclusterP, nXP, nfixedP, nrandomP, randomIntP, indbA);


    // Grid and log(grid) of values where to predict
    // ===============================================
    // cumngrid ....... helping variable counting cumulative sums of ngridM
    // **gridM ........ grid as given by the user
    // **loggridM ..... log(grid - time0)
    //
    int cumngrid = 0;
    double** gridM    = (double**) calloc(*nP, sizeof(double*));
    double** loggridM = (double**) calloc(*nP, sizeof(double*));
    if (!gridM || !loggridM) 
      throw returnR("C++ Error: Could not allocate a memory for working space, too many predictive quantities are computed.", 1);
    for (obs = 0; obs < *nP; obs++){
      gridM[obs] = gridA + cumngrid;
      loggridM[obs] = (double*) calloc(ngridM[obs], sizeof(double));
      if (!loggridM[obs])
        throw returnR("C++ Error: Could not allocate a memory for working space, too many predictive quantities are computed.", 1);
      for (i = 0; i < ngridM[obs]; i++){
        if (gridM[obs][i] <= 0.0) throw returnR("C++ Error: Non-positive grid value for predictive survivor supplied.", 1);
        loggridM[obs][i] = (*trans)(gridM[obs][i] - (*time0P));
      }
      cumngrid += ngridM[obs];
    }
    

    // Read simulated values of important parameters
    // ==============================================
    // kmaxP ......... maximal number of mixture components
    // lmixture ...... number of columns in 'mixture.sim' file
    // nD ............ number of unique elements of a D matrix
    // mixtureAwork .. sampled mixtures
    // betaAwork ..... sampled betas
    // DAwork ........ sampled D matrices (variances of random effects)
    // niter_now ..... number of iterations that will be considered now (check)
    //
    const int* kmaxP     = priorParI + 0;
    const int lmixture   = 1 + 3*(*kmaxP);
    double* mixtureAwork = (double*) calloc(niter*lmixture, sizeof(double));                   if (!mixtureAwork) *errP = 1;
    double* betaAwork    = (double*) calloc(*nXP ? niter*(*nXP) : 1, sizeof(double));          if (!betaAwork)    *errP = 1;
    const int nD         = ((*nrandomP) * (*nrandomP + 1)) / 2;
    double* DAwork       = (double*) calloc(*nrandomP ? niter*(1 + nD) : 1, sizeof(double));   if (!DAwork)       *errP = 1;
    if (*errP) 
      throw returnR("C++ Error: Could not allocate a memory for sampled values.", *errP);

    int niter_now;
    readMixtureFromFiles(mixtureAwork, &niter_now, niter, *skipP, *byP, *kmaxP, dir, 
                         "/mixmoment.sim", "/mweight.sim", "/mmean.sim", "/mvariance.sim");
    if (niter_now != nread) throw returnR("Different MCMC sample sizes indicated by mixture files", 1);
    if (*nXP){
      readFromFile(betaAwork, &niter_now, niter, *nXP, 1, *skipP, *byP, dir, "/beta.sim", 0);
      if (niter_now != nread) throw returnR("Different MCMC sample sizes indicated 'beta.sim'", 1);
    }
    if (*nrandomP){
      readFromFile(DAwork, &niter_now, niter, nD, 1, *skipP, *byP, dir, "/D.sim", 1);              // skip first column
      if (niter_now != nread) throw returnR("Different MCMC sample sizes indicated by 'D.sim'", 1);
    }

    // Space for simulated quantities before they are written to files or used in other way
    // =====================================================================================
    // ETsA[obs][iter] .................... predictive means of survivor times
    // TsA[obs][iter] ..................... predictive survivor times
    // SurvA[obs][grid.value][iter] ....... predictive survivor function
    // hazardA[obs][grid.value][iter] ..... predictive hazard function
    // cumhazardA[obs][grid.value][iter] .. predictive cumulative hazard
    //
    double** ETsA        = (double**) calloc((*predictETP)*(*nP), sizeof(double*));             if (!ETsA)       *errP = 1;
    double** TsA         = (double**) calloc((*predictTP)*(*nP), sizeof(double*));              if (!TsA)        *errP = 1;
    double*** SurvA      = (double***) calloc((*predictSurvP)*(*nP), sizeof(double**));         if (!SurvA)      *errP = 1;
    double*** hazardA    = (double***) calloc((*predicthazardP)*(*nP), sizeof(double**));       if (!hazardA)    *errP = 1;
    double*** cumhazardA = (double***) calloc((*predictcumhazardP)*(*nP), sizeof(double**));    if (!cumhazardA) *errP = 1;
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for working space.", *errP);

    int nToAlloc = (*onlyAver) ? 1 : nread;
    if (*predictETP){
      for (obs = 0; obs < *nP; obs++){
        ETsA[obs] = (double*) calloc(nToAlloc, sizeof(double));
        if (!ETsA[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
      }
    }
    if (*predictTP){
      for (obs = 0; obs < *nP; obs++){
        TsA[obs] = (double*) calloc(nToAlloc, sizeof(double));
        if (!TsA[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
      }
    }
    if (*predictSurvP){
      for (obs = 0; obs < *nP; obs++){
        SurvA[obs] = (double**) calloc(ngridM[obs], sizeof(double*));
        if (!SurvA[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
	  SurvA[obs][j] = (double*) calloc(nToAlloc, sizeof(double));
          if (!SurvA[obs][j]) throw returnR("C++ Error: Could not allocate a memory for SurvA[obs][j].", 1);
        }
      }
    }
    if (*predicthazardP){
      for (obs = 0; obs < *nP; obs++){
        hazardA[obs] = (double**) calloc(ngridM[obs], sizeof(double*));
        if (!hazardA[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
	  hazardA[obs][j] = (double*) calloc(nToAlloc, sizeof(double));
          if (!hazardA[obs][j]) throw returnR("C++ Error: Could not allocate a memory for working space for hazardA[obs][j].", 1);
        }
      }
    }
    if (*predictcumhazardP){
      for (obs = 0; obs < *nP; obs++){
        cumhazardA[obs] = (double**) calloc(ngridM[obs], sizeof(double*));
        if (!cumhazardA[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          cumhazardA[obs][j] = (double*) calloc(nToAlloc, sizeof(double));
          if (!cumhazardA[obs][j]) throw returnR("C++ Error: Could not allocate a memory for working space for cumhazardA[obs][j].", 1);
        }
      }
    }
    //    if (*predictETP) create2pointer(ETsA, *nP, nToAlloc, -99, true);
    //    if (*predictTP) create2pointer(TsA, *nP, nToAlloc, -99, true);
    //    if (*predictSurvP) create3pointer(SurvA, *nP, ngridM, nToAlloc, -99, true);
    //    if (*predicthazardP) create3pointer(hazardA, *nP, ngridM, nToAlloc, -99, true);
    //    if (*predictcumhazardP) create3pointer(cumhazardA, *nP, ngridM, nToAlloc, -99, true);


    // Space for quantiles and means
    // ==============================
    // ETquant[obs][quantile] .......................
    // Tquant[obs][quantile] ........................
    // Survquant[obs][grid.value][quantile] .........
    // hazardquant[obs][grid.value][quantile] .......
    // cumhazardquant[obs][grid.value][quantile] ....
    //
    double** ETquant         = (double**) calloc((*predictETP)*(*nP), sizeof(double*));             if (!ETquant)        *errP = 1;   
    double** Tquant          = (double**) calloc((*predictTP)*(*nP), sizeof(double*));              if (!Tquant)         *errP = 1;   
    double*** Survquant      = (double***) calloc((*predictSurvP)*(*nP), sizeof(double**));         if (!Survquant)      *errP = 1;
    double*** hazardquant    = (double***) calloc((*predicthazardP)*(*nP), sizeof(double**));       if (!hazardquant)    *errP = 1;
    double*** cumhazardquant = (double***) calloc((*predictcumhazardP)*(*nP), sizeof(double**));    if (!cumhazardquant) *errP = 1;
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for working space.", *errP);

    if (*predictETP){
      for (obs = 0; obs < *nP; obs++){
        ETquant[obs] = (double*) calloc(*nquantileP + 1, sizeof(double));
        if (!ETquant[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (i = 0; i < *nquantileP + 1; i++) ETquant[obs][i] = 0.0;
      }
    }
    if (*predictTP){
      for (obs = 0; obs < *nP; obs++){
        Tquant[obs] = (double*) calloc(*nquantileP + 1, sizeof(double));
        if (!Tquant[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (i = 0; i < *nquantileP + 1; i++) Tquant[obs][i] = 0.0;
      }
    }
    if (*predictSurvP){
      for (obs = 0; obs < *nP; obs++){
        Survquant[obs] = (double**) calloc(ngridM[obs], sizeof(double*));
        if (!Survquant[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          Survquant[obs][j] = (double*) calloc(*nquantileP + 1, sizeof(double));
          if (!Survquant[obs][j]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
          for (i = 0; i < *nquantileP + 1; i++) Survquant[obs][j][i] = 0.0;
        }
      }
    }
    if (*predicthazardP){
      for (obs = 0; obs < *nP; obs++){
        hazardquant[obs] = (double**) calloc(ngridM[obs], sizeof(double*));
        if (!hazardquant[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          hazardquant[obs][j] = (double*) calloc(*nquantileP + 1, sizeof(double));
          if (!hazardquant[obs][j]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
          for (i = 0; i < *nquantileP + 1; i++) hazardquant[obs][j][i] = 0.0;
        }
      }
    }
    if (*predictcumhazardP){
      for (obs = 0; obs < *nP; obs++){
        cumhazardquant[obs] = (double**) calloc(ngridM[obs], sizeof(double*));
        if (!cumhazardquant[obs]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          cumhazardquant[obs][j] = (double*) calloc(*nquantileP + 1, sizeof(double));
          if (!cumhazardquant[obs][j]) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
          for (i = 0; i < *nquantileP + 1; i++) cumhazardquant[obs][j][i] = 0.0;
        }
      }
    }
    //    if (*predictETP) create2pointer(ETquant, *nP, *nquantileP + 1, 0.0, true);
    //    if (*predictTP) create2pointer(Tquant, *nP, *nquantileP + 1, 0.0, true);
    //    if (*predictSurvP) create3pointer(Survquant, *nP, ngridM, *nquantileP + 1, 0.0, true);
    //    if (*predicthazardP) create3pointer(hazardquant, *nP, ngridM, *nquantileP + 1, 0.0, true);
    //    if (*predictcumhazardP) create3pointer(cumhazardquant, *nP, ngridM, *nquantileP + 1, 0.0, true);


    // Space for values that are sampled but not stored
    //  and pointers to current values of parameters
    // =================================================
    // *mixtureM ....... current mixture
    // *kP ............. current mixture components
    // *wM ............. current mixture weights
    // *cumwM .......... cumulative current mixture weights
    // *muM ............ current mixture means
    // *invsigma2M ..... current mixture inv-variances
    // *sigma2M ........ current mixture variances
    // *sigmaM ......... current mixture standard deviations
    // *betaM .......... current betas
    // *Dcm ............ current covariance matrix object 
    // *rM ............. sampled component pertinences
    // *YsM ............ sampled predictive log(times)
    // *regresPredM .... sampled regression predictors (x'beta + z'b)
    // nMHinfo2 ........ here: only working integer
    // *MHinfo2M ....... here: only working space
    // *bM ............. sampled predictive values of random effects
    // *bb ............. object to handle random effects
    // *priorbI ........ helping array to create bb
    // *priorbD ........ helping array to create bb
    // *Eb0dependMix ... 0/1 whether a mean of the random intercept depends on the mixture
    // *Eb0 ............ pointer to a mean of random intercept
    //
    double* mixtureM;
    int kP[1];
    double* wM;
    double* cumwM      = (double*) calloc(*kmaxP, sizeof(double));       if (!cumwM)      *errP = 1;
    double* muM;
    double* invsigma2M = (double*) calloc(*kmaxP, sizeof(double));       if (!invsigma2M) *errP = 1;
    double* sigma2M;
    double* sigmaM     = (double*) calloc(*kmaxP, sizeof(double));       if (!sigmaM)     *errP = 1;
    double* betaM = &dtemp;

    covMatrix* Dcm = new covMatrix;
    if (!Dcm) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    if (*nrandomP) *Dcm = covMatrix(DAwork, nrandomP, tolerCholP, tolerQRP);

    int* rM             = (int*) calloc(*nP, sizeof(int));             if (!rM)          *errP = 1;
    double* YsM         = (double*) calloc(*nP, sizeof(double));       if (!YsM)         *errP = 1;
    double* regresPredM = (double*) calloc(*nP, sizeof(double));       if (!regresPredM) *errP = 1;
    for (i = 0; i < *nP; i++) YsM[i] = 0.0;

    const int nMHinfo2 = 1*(*nclusterP);
    int* MHinfo2M      = (int*) calloc(nMHinfo2, sizeof(int));                        if (!MHinfo2M) *errP = 1;
    double* bM         = (double*) calloc((*nclusterP)*(*nrandomP), sizeof(double));  if (!bM)       *errP = 1;
    bblocks* bb        = new bblocks;                                                 if (!bb)       *errP = 1;
    int* priorbI       = (int*) calloc(7+(*nrandomP), sizeof(int));                   if (!priorbI)  *errP = 1;
    double* priorbD    = (double*) calloc(1+2*nD+(*nrandomP)+1, sizeof(double));      if (!priorbD)  *errP = 1;
    if (*errP)
      throw returnR("C++ Error: Could not allocate a memory for working space.", *errP);
    if (*nrandomP){
      priorbI[0] = *nrandomP;  priorbI[1] = *nclusterP;  priorbI[2] = InvWishart;  priorbI[3] = Gibbs;
      priorbI[4] = 1;          priorbI[5] = *nrandomP;   priorbI[6] = nD;
      for (i = 0; i < *nrandomP; i++) priorbI[7+i] = i;
      priorbD[0] = *nrandomP + 2;
      for (i = 0; i < nD; i++) priorbD[1+i] = 1.0;
      for (i = 0; i < nD+(*nrandomP)+1; i++) priorbD[1+nD+i] = 0.0;
      *bb = bblocks(bM, priorbI, priorbD, MHinfo2M, tolerCholP);
    }
 
    const int Eb0dependMix[1] = {priorParI[2]};
    double Eb0[1] = {0};

 
    // THIS SECTION IS ACTUALLY NOT NECESSARY...
    // Give initial values for unknown parameters:
    // Component pertinence: independently sample according to weights of the first mixture
    // Random effects: their means
    // Y = log(T): sample according to mean given by the first regression vector and comp. pertinences
    // ==============================================================================
    mixtureM = mixtureAwork;
    wM       = mixtureAwork + 1;
    muM      = wM + (*kmaxP);
    sigma2M  = muM + (*kmaxP);
    *kP      = int(mixtureM[0]);
    cumsum(cumwM, wM, kP);
    giveSigmaAndInvsigma2(sigmaM, invsigma2M, sigma2M, kP);

    if (*nXP) betaM = betaAwork;
    if (*nrandomP){
      Dcm->covm = DAwork;
      if (*Eb0dependMix) mixMean(Eb0, kP, wM, muM);
      for (cl = 0; cl < *nclusterP; cl++){
        for (j = 0; j < *nrandomP; j++){
          bM[cl*(*nrandomP) + j] = ((indbinXA[j] == -1) ? (*Eb0) : betaM[indbinXA[j]]);
        }
      }
    }
    for (obs = 0; obs < *nP; obs++) YsM[obs] = 0.0;
    regresPredictor(regresPredM, betaM, bM, XA, clusteriA, randomIntP, indbA, nP, nXP, nrandomP);
    discreteSampler(rM, cumwM, kP, nP, &ONE_INT, &ZERO_INT);
    predictData(YsM, regresPredM, rM, cumwM, muM, sigmaM, Eb0, kP, nP, errorTypeP, randomIntP);


    // ========== Main loop over iterations =========
    int iter, index;
    int nwithoutinfo = 0;
    char write_flag = 'a';
    int backs = 0;

    Rprintf("Iteration ");
    for (iter = 1; iter <= nread; iter++){
      index = (*onlyAver) ? 0 : (iter - 1);

      *kP = int(mixtureM[0]);
      cumsum(cumwM, wM, kP);
      giveSigmaAndInvsigma2(sigmaM, invsigma2M, sigma2M, kP);

      if (*Eb0dependMix) mixMean(Eb0, kP, wM, muM);
      if (*nrandomP)     Dcm->update(tolerCholP, tolerQRP);

      // Sample and update cumulative sums
      if (*nrandomP){
        predictRandom(bM, betaM, Eb0, Dcm, nrandomP, nclusterP, indbinXA, bb->indBlock[0]);
	//        GIBBSproposalRandom(bM, regresResM, bb->covpar[0], bb->chcovpar[0], betaM, Eb0, Dcm, XA, ZZt,
	//                            wM, muM, invsigma2M, rM, randomIntP, nrandomP, nXP, nclusterP, nP, errorTypeP,
	//                            invclusteriA, indbinXA, bb->indBlock[0], tolerCholP);
      }
      regresPredictor(regresPredM, betaM, bM, XA, clusteriA, randomIntP, indbA, nP, nXP, nrandomP);
//writeToFile(regresPredM, 1, *nP, "/home/arnost/temp/", "regresPred.sim", 'a');
      if (*predictETP){
        predictET(ETsA, time0P, index, betaM, wM, muM, sigma2M, Dcm, XA, kP, nP, nXP, indbinXA, randomIntP, nrandomP, errorTypeP);
        cumsumQuantile1(ETquant, ETsA, *nquantileP, *nP, index);
      }

      if (*predictTP){
        predictData(YsM, regresPredM, rM, cumwM, muM, sigmaM, Eb0, kP, nP, errorTypeP, randomIntP);
        Y2T(TsA, YsM, time0P, index, nP, itrans);
        cumsumQuantile1(Tquant, TsA, *nquantileP, *nP, index);
      }
      if (predictShch){
  	 predictSurv(SurvA, hazardA, cumhazardA, index, gridM, loggridM, time0P, regresPredM, rM, wM, muM, sigmaM,
                     Eb0, kP, nP, ngridM, errorTypeP, randomIntP, predicthazardP, predictcumhazardP);
         if (*predictSurvP) cumsumQuantile2(Survquant, SurvA, *nquantileP, *nP, ngridM, index);
         if (*predicthazardP) cumsumQuantile2(hazardquant, hazardA, *nquantileP, *nP, ngridM, index);
         if (*predictcumhazardP) cumsumQuantile2(cumhazardquant, cumhazardA, *nquantileP, *nP, ngridM, index);
      }
      nwithoutinfo++;

      if (nwithoutinfo == nwrite){
        for (i = 0; i < backs; i++) Rprintf("\b");
        Rprintf("%d", iter);
        backs = int(log10(double(iter))) + 1;
        nwithoutinfo = 0;
      } 

      // Shift pointers for previously sampled values
      mixtureM += lmixture;
      wM       = mixtureM + 1;
      muM      = wM + (*kmaxP);
      sigma2M  = muM + (*kmaxP);
      if (*nrandomP) Dcm->covm += nD;
      if (*nXP)      betaM     += (*nXP);
    }    // end of the simulation

    if (nwithoutinfo > 0){
      for (i = 0; i < backs; i++) Rprintf("\b");
      Rprintf("%d", iter);
    } 

    // Write to files simulated values
    if (store){
      Rprintf("\nStoring simulated values.");
      writeToFiles2(ETsA, TsA, SurvA, hazardA, cumhazardA, nread, dir, write_flag, nP, ngridM,
                    storeETP, storeTP, storeSurvP, storehazardP, storecumhazardP);
    }

    // Compute means and quantiles
    Rprintf("\nComputing quantiles.");
    if (*predictETP) meanQuantile1(ETquant, ETsA, quantileA, indquant1, indquant2, *nP, *nquantileP, nread);
    if (*predictTP) meanQuantile1(Tquant, TsA, quantileA, indquant1, indquant2, *nP, *nquantileP, nread);
    if (*predictSurvP) meanQuantile2(Survquant, SurvA, quantileA, indquant1, indquant2, *nP, ngridM, *nquantileP, nread);
    if (*predicthazardP) meanQuantile2(hazardquant, hazardA, quantileA, indquant1, indquant2, *nP, ngridM, *nquantileP, nread);
    if (*predictcumhazardP) meanQuantile2(cumhazardquant, cumhazardA, quantileA, indquant1, indquant2, *nP, ngridM, *nquantileP, nread);

    // Write to files quantiles and means
    Rprintf("\nStoring quantiles.");
    writeToFiles3(ETquant, Tquant, Survquant, hazardquant, cumhazardquant, *nquantileP, dir, write_flag, nP, ngridM,
                  predictETP, predictTP, predictSurvP, predicthazardP, predictcumhazardP);
    Rprintf("\n");

    PutRNGstate();

    // Cleaning
    if (!(*onlyAver)){
      free(indquant1);         
      free(indquant2);
    }

    free(clusteriA);         
    delete [] invclusteriA;     
    free(indbinXA);

    for (j = 0; j < *nP; j++) free(ZZt[j]);
    free(ZZt);               
    free(diagIZZt);

    for (obs = 0; obs < *nP; obs++) free(loggridM[obs]);
    free(loggridM);            
    free(gridM);

    free(mixtureAwork);        
    free(betaAwork);        
    free(DAwork);

    if (*predictETP) for (obs = 0; obs < *nP; obs++){
                       free(ETsA[obs]);
                       free(ETquant[obs]);
                     }
    if (*predictTP) for (obs = 0; obs < *nP; obs++){
                       free(TsA[obs]);
                       free(Tquant[obs]);
                     }
    if (*predictSurvP) for (obs = 0; obs < *nP; obs++){ 
                         for (i = 0; i < ngridM[obs]; i++){
                           free(SurvA[obs][i]);
                           free(Survquant[obs][i]);
                         }
                         free(SurvA[obs]);
                         free(Survquant[obs]); 
                       }
    if (*predicthazardP) for (obs = 0; obs < *nP; obs++){
                           for (i = 0; i < ngridM[obs]; i++){
                             free(hazardA[obs][i]);
                             free(hazardquant[obs][i]);
                           }
                           free(hazardA[obs]);
                           free(hazardquant[obs]); 
                         }
    if (*predictcumhazardP) for (obs = 0; obs < *nP; obs++){
                              for (i = 0; i < ngridM[obs]; i++){
                                free(cumhazardA[obs][i]);
                                free(cumhazardquant[obs][i]);
                              }
                              free(cumhazardA[obs]);
                              free(cumhazardquant[obs]);
                            }

    free(ETsA);     free(TsA);     free(SurvA);      free(hazardA);      free(cumhazardA);
    free(ETquant);  free(Tquant);  free(Survquant);  free(hazardquant);  free(cumhazardquant);

    free(cumwM);    free(invsigma2M);    free(sigmaM);
    free(rM);       free(YsM);           free(regresPredM);
    free(bM);       free(MHinfo2M);      free(priorbI);       free(priorbD);
    delete bb;      delete Dcm;

    return;
  }                     // end of try
  catch(returnR rr){
    *errP = rr.errflag();
    PutRNGstate();
    return;
  }

}    // end of the function predictive

}    // end of extern "C"
