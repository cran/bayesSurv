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

#include "bayessurvreg.h"

extern "C" {

using namespace std;

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
//         nToStore ................ number of values that will be stored before written to files
//                                   (currently ignored, always changed to nwrite)
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
predictive(int* errorTypeP,     char** dirP,
           int* dimsP,          int* dims2P, 
           double* XA,          int* indbA,
           double* quantileA,   double* gridA,
           int* priorParI,      double* priorParD,
           int* nsimulP,        int* predictP,    int* storeP,
           double* tolersP,     int* errP)
{
  try{
    int obs, i, j, cl;

    if (!(*errorTypeP == Mixture || *errorTypeP == Spline))
      throw returnR("C++ Error: unimplemented 'errorType' in a call to function 'predictive'.", 1);

    GetRNGstate();

    *errP = 0;
    string dir = *dirP;

    // Transformation and inverse transformation of the response
    // ==========================================================
    double (*itrans)(const double);                itrans = exp;
    double (*trans)(const double);                 trans = log;

    // Starting time for the survival model
    // ====================================
    double* time0P = priorParD + 0;

    // Dimensionality parameters
    // =========================
    int* nP = dimsP;
    int* nclusterP = dimsP + 1;
    int* nwithinA = dimsP + 2;
    int* nYP = dimsP + 2 + (*nclusterP);
    int* nXP = dimsP + 2 + (*nclusterP) + 1;
    int* nfixedP = dimsP + 2 + (*nclusterP) + 2;
    int* nrandomP = dimsP + 2 + (*nclusterP) + 3;
    int* randomIntP = dimsP + 2 + (*nclusterP) + 4;
    int nToStore = dimsP[2 + (*nclusterP) + 5];

    // Additional dimensionality parameters
    // =====================================
    int* nquantileP = dims2P + 0;
    int* ngridM = dims2P + 1;

    // Parameters of done or future simulation
    // ========================================
    int niter = nsimulP[0]; 
    int nthin = nsimulP[1];
    int nwrite = nsimulP[2];         nToStore = nwrite;

    // Tolerances for decompositions
    // =============================
    double* tolerCholP = tolersP;
    double* tolerQRP = tolersP + 1;

    // What to predict?
    // =================
    int* predictETP = predictP + 0;
    int* predictTP = predictP + 1;
    int* predictSurvP = predictP + 2;
    int* predicthazardP = predictP + 3;
    int* predictcumhazardP = predictP + 4;
    if (*predicthazardP || *predictcumhazardP) *predictSurvP = 1;
    bool predictShch = (*predictSurvP) || (*predicthazardP) || (*predictcumhazardP);

    // What to store?
    // ================
    int* storeETP = storeP + 0;             if (!(*predictETP)) *storeETP = 0;
    int* storeTP = storeP + 1;              if (!(*predictTP)) *storeTP = 0;
    int* storeSurvP = storeP + 2;           if (!(*predictSurvP)) *storeSurvP = 0;
    int* storehazardP = storeP + 3;         if (!(*predicthazardP)) *storehazardP = 0;
    int* storecumhazardP = storeP + 4;      if (!(*predictcumhazardP)) *storecumhazardP = 0;
    bool store = (*storeETP) || (*storeTP) || (*storeSurvP) || (*storehazardP) || (*storecumhazardP);

    // Indeces of quantile values in sampled chain (indexing starting from 0)
    // ======================================================================
    // indquant1, indquant2 ..... quantile = q*sample[indquant1] + (1-q)sample[indquant2]
    // help ..................... helping number
    //
    int* indquant1 = new int[*nquantileP];
    int* indquant2 = new int[*nquantileP];
    double help;
    for (i = 0; i < *nquantileP; i++){
      if (quantileA[i] < 0 || quantileA[i] > 1) throw returnR("C++ Error: Incorrect quantile values supplied.", 1);
      if (quantileA[i] <= 0) indquant1[i] = indquant2[i] = 0;
      else if (quantileA[i] >= 1) indquant1[i] = indquant2[i] = niter - 1;
           else{
             help = quantileA[i] * double(niter);
             if (fabs(help - floor(help + 1e-8)) < 1e-8){   // help is integer
               indquant1[i] = int(floor(help)) - 1;
               indquant2[i] = int(floor(help));
             }
             else{
               indquant1[i] = indquant2[i] = int(floor(help));
             }
           }
    }

//    cout << "indquant1 = "; printArray(indquant1, nquantileP);
//    cout << "indquant2 = "; printArray(indquant2, nquantileP);


    // Create data matrices and variables
    // ==================================
    // *clusteriA ........... indicators of clusters for observations
    // *invclusteriA ........
    // *indbinXA
    // **ZZt ................ this is finally not needed
    // *diagIZZt ............ this is finally not needed
    //
    int* clusteriA = new int[*nP];                           if (clusteriA == NULL) *errP = 1;
    List<int>* invclusteriA = new List<int>[*nclusterP];     if (invclusteriA == NULL) *errP = 1;
    int* indbinXA = new int[(*nrandomP ? *nrandomP : 1)];    if (indbinXA == NULL) *errP = 1;
    double** ZZt = new double*[*nP];                         if (ZZt == NULL) *errP = 1;
    int* diagIZZt = new int[*nrandomP];                      if (diagIZZt == NULL) *errP = 1;
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (j = 0; j < *nP && !(*errP); j++){
      ZZt[j] = new double[*nrandomP * (*nrandomP)];          if (ZZt[j] == NULL) *errP = 1;
    }
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);

    createDataShort(nwithinA, clusteriA, invclusteriA, XA, ZZt, diagIZZt, indbinXA, 
                    nP, nclusterP, nXP, nfixedP, nrandomP, randomIntP, indbA);


    // Grid and log(grid) of values where to predict
    // ===============================================
    // cumngrid ....... helping variable counting cumulative sums of ngridM
    // **gridM ........ grid as given by the user
    // **loggridM ..... log(grid - time0)
    //
    int cumngrid = 0;
    double** gridM = new double*[*nP];
    double** loggridM = new double*[*nP];
    if (gridM == NULL || loggridM == NULL) 
      throw returnR("C++ Error: Could not allocate a memory for working space, too many predictive quantities are computed.", 1);
    for (obs = 0; obs < *nP; obs++){
      gridM[obs] = gridA + cumngrid;
      loggridM[obs] = new double[ngridM[obs]];
      if (loggridM[obs] == NULL)
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
    //
    int* kmaxP = priorParI + 0;
    int lmixture = 1 + 3*(*kmaxP);
    double* mixtureAwork = new double[niter*lmixture];                 if (mixtureAwork == NULL) *errP = 1;
    double* betaAwork = new double[(*nXP ? niter*(*nXP) : 1)];         if (betaAwork == NULL) *errP = 1;
    int nD = ((*nrandomP) * (*nrandomP + 1)) / 2;
    double* DAwork = new double[(*nrandomP ? niter*(1 + nD) : 1)];     if (DAwork == NULL) *errP = 1;
    if (*errP) 
      throw returnR("C++ Error: Could not allocate a memory for sampled values.", *errP);

    readMixtureFromFiles(mixtureAwork, niter, *kmaxP, dir, "/mixmoment.sim", "/mweight.sim", "/mmean.sim", "/mvariance.sim", 1);    
    if (*nXP) readFromFile(betaAwork, niter, *nXP, dir, "/beta.sim", 0, 1);
    if (*nrandomP) readFromFile(DAwork, niter, nD, dir, "/D.sim", 1, 1);              // skip first column


    // Space for simulated quantities before they are written to files or used in other way
    // =====================================================================================
    // ETsA[obs][iter] .................... predictive means of survivor times
    // TsA[obs][iter] ..................... predictive survivor times
    // SurvA[obs][grid.value][iter] ....... predictive survivor function
    // hazardA[obs][grid.value][iter] ..... predictive hazard function
    // cumhazardA[obs][grid.value][iter] .. predictive cumulative hazard
    //
    double** ETsA = new double*[(*predictETP)*(*nP)];                      if (ETsA == NULL) *errP = 1;
    double** TsA = new double*[(*predictTP)*(*nP)];                        if (TsA == NULL) *errP = 1;
    double*** SurvA = new double**[(*predictSurvP)*(*nP)];                 if (SurvA == NULL) *errP = 1;
    double*** hazardA = new double**[(*predicthazardP)*(*nP)];             if (hazardA == NULL) *errP = 1;
    double*** cumhazardA = new double**[(*predictcumhazardP)*(*nP)];       if (cumhazardA == NULL) *errP = 1;
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for working space.", *errP);
    if (*predictETP){
      for (obs = 0; obs < *nP; obs++){
        ETsA[obs] = new double[niter];
        if (ETsA[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
      }
    }
    if (*predictTP){
      for (obs = 0; obs < *nP; obs++){
        TsA[obs] = new double[niter];
        if (TsA[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
      }
    }
    if (*predictSurvP){
      for (obs = 0; obs < *nP; obs++){
        SurvA[obs] = new double*[ngridM[obs]];
        if (SurvA[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          SurvA[obs][j] = new double[niter];
          if (SurvA[obs][j] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        }
      }
    }
    if (*predicthazardP){
      for (obs = 0; obs < *nP; obs++){
        hazardA[obs] = new double*[ngridM[obs]];
        if (hazardA[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          hazardA[obs][j] = new double[niter];
          if (hazardA[obs][j] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        }
      }
    }
    if (*predictcumhazardP){
      for (obs = 0; obs < *nP; obs++){
        cumhazardA[obs] = new double*[ngridM[obs]];
        if (cumhazardA[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          cumhazardA[obs][j] = new double[niter];
          if (cumhazardA[obs][j] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        }
      }
    }
    //    if (*predictETP) create2pointer(ETsA, *nP, niter, -99, true);
    //    if (*predictTP) create2pointer(TsA, *nP, niter, -99, true);
    //    if (*predictSurvP) create3pointer(SurvA, *nP, ngridM, niter, -99, true);
    //    if (*predicthazardP) create3pointer(hazardA, *nP, ngridM, niter, -99, true);
    //    if (*predictcumhazardP) create3pointer(cumhazardA, *nP, ngridM, niter, -99, true);


    // Space for quantiles and means
    // ==============================
    // ETquant[obs][quantile] .......................
    // Tquant[obs][quantile] ........................
    // Survquant[obs][grid.value][quantile] .........
    // hazardquant[obs][grid.value][quantile] .......
    // cumhazardquant[obs][grid.value][quantile] ....
    //
    double** ETquant = new double*[(*predictETP)*(*nP)];                    if (ETquant == NULL) *errP = 1;   
    double** Tquant = new double*[(*predictTP)*(*nP)];                      if (Tquant == NULL) *errP = 1;   
    double*** Survquant = new double**[(*predictSurvP)*(*nP)];              if (Survquant == NULL) *errP = 1;
    double*** hazardquant = new double**[(*predicthazardP)*(*nP)];          if (hazardquant == NULL) *errP = 1;
    double*** cumhazardquant = new double**[(*predictcumhazardP)*(*nP)];    if (cumhazardquant == NULL) *errP = 1;
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for working space.", *errP);

    if (*predictETP){
      for (obs = 0; obs < *nP; obs++){
        ETquant[obs] = new double[*nquantileP + 1];
        if (ETquant[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (i = 0; i < *nquantileP + 1; i++) ETquant[obs][i] = 0.0;
      }
    }
    if (*predictTP){
      for (obs = 0; obs < *nP; obs++){
        Tquant[obs] = new double[*nquantileP + 1];
        if (Tquant[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (i = 0; i < *nquantileP + 1; i++) Tquant[obs][i] = 0.0;
      }
    }
    if (*predictSurvP){
      for (obs = 0; obs < *nP; obs++){
        Survquant[obs] = new double*[ngridM[obs]];
        if (Survquant[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          Survquant[obs][j] = new double[*nquantileP + 1];
          if (Survquant[obs][j] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
          for (i = 0; i < *nquantileP + 1; i++) Survquant[obs][j][i] = 0.0;
        }
      }
    }
    if (*predicthazardP){
      for (obs = 0; obs < *nP; obs++){
        hazardquant[obs] = new double*[ngridM[obs]];
        if (hazardquant[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          hazardquant[obs][j] = new double[*nquantileP + 1];
          if (hazardquant[obs][j] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
          for (i = 0; i < *nquantileP + 1; i++) hazardquant[obs][j][i] = 0.0;
        }
      }
    }
    if (*predictcumhazardP){
      for (obs = 0; obs < *nP; obs++){
        cumhazardquant[obs] = new double*[ngridM[obs]];
        if (cumhazardquant[obs] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
        for (j = 0; j < ngridM[obs]; j++){
          cumhazardquant[obs][j] = new double[*nquantileP + 1];
          if (cumhazardquant[obs][j] == NULL) throw returnR("C++ Error: Could not allocate a memory for working space.", 1);
          for (i = 0; i < *nquantileP + 1; i++) cumhazardquant[obs][j][i] = 0.0;
        }
      }
    }
    //    if (*predictETP) create2pointer(ETquant, *nP, *nquantileP, 0.0, true);
    //    if (*predictTP) create2pointer(Tquant, *nP, *nquantileP, 0.0, true);
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
    double* cumwM = new double[*kmaxP];                 if (cumwM == NULL) *errP = 1;
    double* muM;
    double* invsigma2M = new double[*kmaxP];            if (invsigma2M == NULL) *errP = 1;
    double* sigma2M;
    double* sigmaM = new double[*kmaxP];                if (sigmaM == NULL) *errP = 1;
    double* betaM;

    covMatrix* Dcm = new covMatrix;
    if (Dcm == NULL) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    if (*nrandomP) *Dcm = covMatrix(DAwork, nrandomP, tolerCholP, tolerQRP);

    int* rM = new int[*nP];                             if (rM == NULL) *errP = 1;
    double* YsM = new double[*nP];                      if (YsM == NULL) *errP = 1;
    for (i = 0; i < *nP; i++) YsM[i] = 0.0;
    double* regresPredM = new double[*nP];               if (regresPredM == NULL) *errP = 1;

    const int nMHinfo2 = 1*(*nclusterP);
    int* MHinfo2M = new int[nMHinfo2];                  if (MHinfo2M == NULL) *errP = 1;
    double* bM = new double[(*nclusterP)*(*nrandomP)];  if (bM == NULL) *errP = 1;
    bblocks* bb = new bblocks;                          if (bb == NULL) *errP = 1;
    int* priorbI = new int[7+(*nrandomP)];              if (priorbI == NULL) *errP = 1;
    double* priorbD = new double[1+2*nD+(*nrandomP)+1]; if (priorbD == NULL) *errP = 1;
    if (*errP)
      throw returnR("C++ Error: Could not allocate a memory for working space, buy more memory.", *errP);
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
    *kP = int(mixtureAwork[0]);
    wM = mixtureAwork + 1;
    muM = mixtureAwork + 1 + (*kmaxP);
    sigma2M = mixtureAwork + 1 + 2*(*kmaxP);
    giveSigmaAndInvsigma2(sigmaM, invsigma2M, sigma2M, kP);
    cumsum(cumwM, wM, kP);

    if (*nXP) betaM = betaAwork + 0;
    if (*nrandomP){
      if (*Eb0dependMix){
        *kP = int(mixtureAwork[0]);
        mixMean(Eb0, kP, wM, muM);
      }
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
    int iter;
    int nwithoutinfo = 0;
    char write_flag = 'a';
    int backs = 0;

    Rprintf("Iteration ");
    for (iter = 1; iter <= niter; iter++){

//cout << "iter = " << iter << endl;
      // Shift pointers for previously sampled values
      mixtureM = mixtureAwork + (iter - 1)*lmixture;
      *kP = int(mixtureM[0]); 
      wM = mixtureM + 1;
      cumsum(cumwM, wM, kP);
      muM = mixtureM + 1 + (*kmaxP);
      sigma2M = mixtureM + 1 + 2*(*kmaxP);
      giveSigmaAndInvsigma2(sigmaM, invsigma2M, sigma2M, kP);
      if (*Eb0dependMix) mixMean(Eb0, kP, wM, muM);
      if (*nrandomP){
        Dcm->covm = DAwork + (iter - 1)*nD;
        Dcm->update(tolerCholP, tolerQRP);
      }
      if (*nXP) betaM = betaAwork + (iter - 1)*(*nXP);


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
        predictET(ETsA, time0P, iter - 1, betaM, wM, muM, sigma2M, Dcm, XA, kP, nP, nXP, indbinXA, randomIntP, nrandomP, errorTypeP);
        cumsumQuantile1(ETquant, ETsA, *nquantileP, *nP, iter - 1);
      }

      if (*predictTP){
        predictData(YsM, regresPredM, rM, cumwM, muM, sigmaM, Eb0, kP, nP, errorTypeP, randomIntP);
        Y2T(TsA, YsM, time0P, iter - 1, nP, itrans);
        cumsumQuantile1(Tquant, TsA, *nquantileP, *nP, iter - 1);
      }
      if (predictShch){
  	 predictSurv(SurvA, hazardA, cumhazardA, iter - 1, gridM, loggridM, time0P, regresPredM, rM, wM, muM, sigmaM,
                     Eb0, kP, nP, ngridM, errorTypeP, randomIntP, predicthazardP, predictcumhazardP);
         if (*predictSurvP) cumsumQuantile2(Survquant, SurvA, *nquantileP, *nP, ngridM, iter - 1);
         if (*predicthazardP) cumsumQuantile2(hazardquant, hazardA, *nquantileP, *nP, ngridM, iter - 1);
         if (*predictcumhazardP) cumsumQuantile2(cumhazardquant, cumhazardA, *nquantileP, *nP, ngridM, iter - 1);
      }
      nwithoutinfo++;

      if (nwithoutinfo == nwrite){
        for (i = 0; i < backs; i++) Rprintf("\b");
        Rprintf("%d", iter);
        backs = int(log10(double(iter))) + 1;
        nwithoutinfo = 0;
      } 
    }    // end of the simulation

    // Write to files simulated values
    if (store){
      Rprintf("\nStoring simulated values.");
      writeToFiles2(ETsA, TsA, SurvA, hazardA, cumhazardA, niter, dir, write_flag, nP, ngridM,
                    storeETP, storeTP, storeSurvP, storehazardP, storecumhazardP);
    }

    // Compute means and quantiles
    Rprintf("\nComputing quantiles.");
    if (*predictETP) meanQuantile1(ETquant, ETsA, quantileA, indquant1, indquant2, *nP, *nquantileP, niter);
    if (*predictTP) meanQuantile1(Tquant, TsA, quantileA, indquant1, indquant2, *nP, *nquantileP, niter);
    if (*predictSurvP) meanQuantile2(Survquant, SurvA, quantileA, indquant1, indquant2, *nP, ngridM, *nquantileP, niter);
    if (*predicthazardP) meanQuantile2(hazardquant, hazardA, quantileA, indquant1, indquant2, *nP, ngridM, *nquantileP, niter);
    if (*predictcumhazardP) meanQuantile2(cumhazardquant, cumhazardA, quantileA, indquant1, indquant2, *nP, ngridM, *nquantileP, niter);

    // Write to files quantiles and means
    Rprintf("\nStoring quantiles.");
    writeToFiles3(ETquant, Tquant, Survquant, hazardquant, cumhazardquant, *nquantileP, dir, write_flag, nP, ngridM,
                  predictETP, predictTP, predictSurvP, predicthazardP, predictcumhazardP);
    Rprintf("\n");

    PutRNGstate();

    // Cleaning
    for (obs = 0; obs < *nP; obs++) delete [] loggridM[obs];
    delete [] loggridM;         delete [] gridM;

    delete [] mixtureAwork;     delete [] betaAwork;    delete [] DAwork;

    delete bb;                  delete [] bM;           delete [] MHinfo2M;
    delete Dcm;

    delete [] clusteriA;        delete [] invclusteriA; delete [] indbinXA;

    for (j = 0; j < *nP; j++) delete [] ZZt[j];
    delete [] ZZt;                                       delete [] diagIZZt;

    delete [] cumwM;            delete [] invsigma2M;       delete [] sigmaM;
    delete [] rM;               delete [] YsM;              delete [] regresPredM;
 
    delete [] indquant1;        delete [] indquant2;

    if (*predictETP) for (obs = 0; obs < *nP; obs++){
                       delete [] ETsA[obs];
                       delete [] ETquant[obs];
                    }
    if (*predictTP) for (obs = 0; obs < *nP; obs++){
                       delete [] TsA[obs];
                       delete [] Tquant[obs];
                    }
    if (*predictSurvP) for (obs = 0; obs < *nP; obs++){ 
                         for (i = 0; i < ngridM[obs]; i++){
                           delete [] SurvA[obs][i];
                           delete [] Survquant[obs][i];
                         }
                         delete [] SurvA[obs];
                         delete [] Survquant[obs]; 
                       }
    if (*predicthazardP) for (obs = 0; obs < *nP; obs++){
                           for (i = 0; i < ngridM[obs]; i++){
                             delete [] hazardA[obs][i];
                             delete [] hazardquant[obs][i];
                           }
                           delete [] hazardA[obs];
                           delete [] hazardquant[obs]; 
                         }
    if (*predictcumhazardP) for (obs = 0; obs < *nP; obs++){
                              for (i = 0; i < ngridM[obs]; i++){
                                delete [] cumhazardA[obs][i];
                                delete [] cumhazardquant[obs][i];
                              }
                              delete [] cumhazardA[obs];
                              delete [] cumhazardquant[obs];
                            }
    delete [] TsA;          delete [] Tquant;
    delete [] SurvA;        delete [] Survquant;
    delete [] hazardA;      delete [] hazardquant;
    delete [] cumhazardA;   delete [] cumhazardquant;

    return;
  }  // end of try
  catch(returnR rr){
    *errP = rr.errflag();
    PutRNGstate();
    return;
  }

}    // end of the function predictive

}    // end of extern "C"
