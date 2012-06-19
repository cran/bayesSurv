// Bayesian survival regression, type 1 (classical mixture as an error)
// ====================================================================

//   * mixture of normals with unknown number of components is used 
//     as en error distribution

//   * main function to run the McMC simulation

// 03/11/2003: start working on it:
// 13/02/2004: first working version with all McMC steps:
// 18/02/2004: possibility for auxiliary variables to propose reversible moves:
// 13/03/2004: change handling of a random intercept (allow that its mean depends on a mixture)
// 09/12/2004: some checks for possible NaN in mixture specification added after encountering 
//             such problem with EORTC breast cancer data
//
// ==================================================================

#include "bayessurvreg1.h"

extern "C" {

using namespace std;

// PARAMETERS:
//
// dirP ........ directory where to store files with simulated values
//
// Dimension parameters
// ====================
// dimsP ....... dimensionality parameters     (dim = nclusters + 7)
//         nP .......... number of observations                                       (dim = 1)
//         nclusterP ... number of clusters                                           (dim = 1) 
//         nwithinA .... numbers of observations within clusters                      (dim = *nclustersP)
//         nYP ......... number of columns in the response matrix (2 or 3)            (dim = 1)
//         nXP ......... number of columns in the design matrix                       (dim = 1)
//         nfixedP ..... number of fixed effects                                      (dim = 1)
//         nrandomP .... number of random effects (including the random intercept)    (dim = 1)               
//         randomIntP .. 0/1, is a random intercept in the model or not               (dim = 1)
//         nToStore .... number of values that will be stored before written to files (dim = 1)
//                       (this is here from historical reasons, it is curently ignored and always changed to nwrite)
//
// Data and design matrices
// ========================
// YA ......... response matrix                                    (dim = *nP x *nyP)
// XA ......... design matrix for both fixed and random effects    (dim = *nP x *nXP)
// indbA ...... an array of length *nXP, indicating where to find random effects
//         indbA[j] = -1 if the jth column of X corresponds to the fixed effect
//         indbA[j] = k if the jth column of X corresponds to the kth random effect
//         i.e. values of -1, 1, ..., nrandom - 1 appear in this array if *randomIntP
//         and values of -1, 0, ..., nrandom - 1 appear in this array if !(*randomIntP)  
//
// Initial values and arrays to store the current simulated results
// =================================================================             
// iterM ...... index of the initial iteration                                                             (1)
// loglikM .... used to store a log-likelihood and randomloglik                                            (2)
// mixtureM ... mixture in the form k, w1, ..., w(kmax), mu1, ... mu[kmax], sigma2[1], ... sigma2[kmax] 
//         at rows                                                                                         (3*kmax+1) 
// mixMomentM . used to store mixture mean and standard deviation                                          (2)
// betaM ...... initial values of the regression parameters                                                (*nXP x 1)
// bM ......... initial values of the random effects                                                       (*nclusterP times *nrandomP)
// DM ......... initial lower triangle of the D matrix                                      (0.5 times (*nrandomP) times (*nrandomP+1))
// rM ......... initial values of components pertinence                                                    (*nP)
//          R indeces (1, 2, ...) on input and on return
//          C++ indeces (0, 1, ...) during the computation
// YsM ........ initial sampled latent responses                                                           (*nP)
// otherpM .... initial values of the scale of mixture inv-variances                                       (1)
// uM ......... initial value of the canonical proposal vector for the split-combine move                  (3)
//
// Specification of priors and parameters for Metropolis-Hastings updates
// =======================================================================
// priorParI[3] ..... integer prior parameters
//          kmax = maximal number of mixture components 
//          priorFork = indicator of the prior for k (0 = Poisson, 1 = uniform)  
//          Eb0dependMix = indicator whether Eb0 should depend on a mixture or not
// priorParD[2*kmax + 7] ..... double prior parameters
//          piSplit = probabilities of the split move
//          piBirth = probabilities of the birth move
//          lambda = prior hyperparameter for the number of components 
//          delta = prior parameter for component weights   
//          xi = prior mean of component means   
//          kappa = prior variance of component means
//          zeta = prior shape of component inverse-variances    
//          g = prior shape of the scale of comp. inverse-variances 
//          h = prior rate of the scale of comp. inverse-variances
//
// revJumpParI[ ] ... integer parameters identifying the way the proposals of reversible jumps
//                    will be generated
//          revJumpAlg = type of the algorithm to produce canonical variables for reversible jumps
//          revJumpTransSC = type of the transformation to produce proposed weight, mean and variance from canonical variables
//                 during the split-combine move
//          revJumpTransBD = type of the transformation to produce proposed weight, mean and variance from canonical variables
//                 during the birth-death move
// revJumpParD ... double parameters identifying the way the proposals of reversible jumps
//                 will be generated
//          mrParm = epsilon and delta for a moody ring                            (2)
//          transParmuSC = transformation parameters for a split-combine move 
//                 (typically parameters of beta distributions)                    (6)          
//          transParmuBD = transformsation parameters for a birth-death move      
//                 (with this version, give only space to store 6 doubles)         (6)
//
// priorBetaI .... integer prior and simulation parameters for beta
//          nBlocks = number of blocks in which beta will be updated                       (1)
//          nBeta = nX = number of beta parameters                                         (1)
//          nInBlock = numbers of betas within blocks                                      (nBlocks)
//          max(nInBlock) = length of the biggest block                                    (1)
//          lcovparLV = lengths of lower triangles of covariance matrices 
//                  for each block                                                         (nBlocks)
//          indBlockLV = indeces of betas within blocks (R indeces)                        (nX)
//          typeUpd = update type for each block                                           (nBlocks)
// priorBetaD .... double prior and simulation parameters for beta
//           mean.prior                                                                    (nX)
//           var.prior                                                                     (nX)
//           mean.sampled = mean of up to now sampled values
//                  this is updated for components with adaptive metropolis algorithm      (nX)
//           halfRangeUnif = specification of a uniform distribution for proposal          (nX)
//           covparLV = lower triangles of (initial) covariance matrices for
//                  normal proposals,
//                  this is updated if adaptive Metropolis algorithm is used               (sum(lcovparLV))
//           weightUnif = weights of the uniform component for proposal
//                  for each block                                                         (nBlocks)
//           eps = epsilon numbers for the adaptive Metropolis algorithm
//                  for each block                                                         (nBlocks)
//           sdNum = s_d constants for the adaptive Metropolis algorithm
//                  for each dimension up to max(nInBlock)                                 (maxnInBlock)
//
// priorbI ....... integer prior and simulation parameters for random effects
//           nrandom
//           ncluster
//           priorD
//           typeUpd
//           nBlocks
//           nInBlock
//           lcovparLV
//           indBlockLV
// priorbD ....... double prior and simulation parameters for random effects 
//          tau = prior degrees of freedom of the inverse Wishart distr. of a covariance matrix D
//          S = prior scale matrix of the inverse Wishart distr. of a covariance matrix D (lower triangle)
//           covparLV
//           halfRangeUnif
//           weightUnif
//
// Parameters of the simulation and the rest
// =========================================
// nsimulP .... vector with niter, nthin, nwrite where                                                 (3)
//           niter = number of simulated values (after discadring thinned ones)
//           nthin = thinning interval
//           nwrite = how many simulated values at once are to be written to files 
// storeP ..... flags whether to store Ys, r, b, propVec, MH info for random effects, regresRes        (6)
// tolersP .... tolerances for decompositions                                                          (2)       
//           tolerCholP . tolerance for the cholesky decomposition                                     (1)
//           tolerQRP ... tolerance for the QR decomposition                                           (1)
// errP ....... error flag                                                                             (1)
//
void 
bayessurvreg1(char** dirP,
              int* dimsP,
              double* YA,           double* XA,            int* indbA,
              int* iterM,           double* loglikM,       double* mixtureM,   double* mixMomentM,
              double* betaM,        double* bM,            double* DM,
              int* rM,              double* YsM,           double* otherpM,
              double* uM,
              int* priorParI,       double* priorParD,
              int* revJumpParI,     double* revJumpParD,
              int* priorBetaI,      double* priorBetaD,
              int* priorbI,         double* priorbD,
              int* nsimulP,         int* storeP,           double* tolersP,
              int* errP 
             )
{  
  try{
    int i, j, l;

    GetRNGstate();

    *errP = 0;
    string dir = *dirP;
    int errorTypeP[1] = {Mixture};

    // Logarithm of the Jacobian of the transformation of the response
    // ===============================================================
    // * this is actually not necessary, I do not add log(Jacobian) to log-likelihood
    double (*logdtrans) (const double);
    logdtrans = logdidenttrans;          


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
    int nToStore = dimsP[2 + (*nclusterP) + 5];          // later it's changed to nwrite

    // Parameters of the simulation
    // =============================
    // storeyP ......... whether to store simulated augmented log(event times)
    // storerP ......... whether to store simulated component pertinences
    // storebP ......... whether to store simulated values of random effects
    // storeuP ......... whether to store simulated canonical proposal vector for the split/combine move
    // storeMHbP ....... whether to store information on acceptance for random effects
    // storeregresResP . whether to store regression residuals
    // niter ........... number of simulated values (after discarding thinned ones)
    // nthin ........... thinning interval
    // nwrite .......... interval for writting to files
    // tolerCholP ...... tolerance for a Cholesky decomposition
    // tolerQRP ........ tolerance for a QR decomposition
    //
    int* storeyP = storeP;
    int* storerP = storeP + 1;
    int* storebP = storeP + 2;
    int* storeuP = storeP + 3;
    int* storeMHbP = storeP + 4;
    int* storeregresResP = storeP + 5;
    int niter = nsimulP[0]; 
    int nthin = nsimulP[1];
    int nwrite = nsimulP[2];                    
    nToStore = nwrite;
    double* tolerCholP = tolersP;
    double* tolerQRP = tolersP + 1;

    // Additional monitoring variables and helping constants
    // ======================================================
    // nMHinfo ...... number of "columns" in MHinfoAwork 
    // nMHinfo2 ..... number of "columns" in MHinfoAwork2
    // ncmSM ........ number of "columns" in mixtureAwork
    // MHinfoM ...... an array to store information on cumulative acceptance/rejectance for some parameters (reversible jumps, regression parameters)
    // MHinfo2M ..... an array to store information on cumulative acceptance/rejectance for additional parameters (random effects)
    //    MHinfoM[0] = number of accepted split/combine moves
    //    MHinfoM[1] = number of proposed split moves
    //    MHinfoM[2] = number of accepted birth/death moves
    //    MHinfoM[3] = number of proposed birth moves
    //    MHinfoM[4] to MHinfoM[4 + betaMH->nBlock() - 1] = number of accepted moves for blocks of beta parameters
    //    MHinfo2M[0] to MHinfo2M[bb->nBlock()*bb->nCluster() - 1] = number of accepted
    //          moves for blocks of random effects, sorted first according to clusters and then according to blocks 
    //
    const int nMHinfo = 4 + priorBetaI[0];
    const int nMHinfo2 = priorbI[4]*(*nclusterP);
    const int ncmSM = 1 + 3*priorParI[0];
    int* MHinfoM = new int[nMHinfo];
    int* MHinfo2M = new int[nMHinfo2];
    if (MHinfoM == NULL || MHinfo2M == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (i = 0; i < nMHinfo; i++) MHinfoM[i] = 0;
    for (i = 0; i < nMHinfo2; i++) MHinfo2M[i] = 0;

    // Allocate a working space to store simulated values before writting them to files
    // ==================================================================================
    // kmaxP .......... maximal number of mixture components
    // iterAwork ...... to store indeces of stored iterations
    // loglikAwork ....
    // mixtureAwork ...
    // betaAwork ......
    // bAwork .........
    // DAwork ......... determinant of D and a lower triangle of D
    // rAwork .........
    // YsAwork ........
    // otherpAwork ....
    // uAwork .........
    // mixMomentAwork .
    // MHinfoAwork .... acceptance ratios for the Metropolis-Hasting 
    //                  and ratio of split and birth moves,
    //                  MHinfoM is copied after each iteration to one row of MHinfoAwork
    // MHinfo2Awork ... similar 
    //
    int* kmaxP = priorParI;
    int* iterAwork = new int[nToStore];                                                              if (iterAwork == NULL) *errP = 1;
    double* loglikAwork = new double[nToStore*2];                                                    if (loglikAwork == NULL) *errP = 1;
    double* mixtureAwork = new double[nToStore*(1+3*(*kmaxP))];                                      if (mixtureAwork == NULL) *errP = 1;
    double* betaAwork = new double[(*nXP ? nToStore*(*nXP) : 1)];                                    if (betaAwork == NULL) *errP = 1;
    double* DAwork = new double[(*nrandomP ? nToStore*(1 + ((*nrandomP)*(*nrandomP + 1))/2) : 1)];   if (DAwork == NULL) *errP = 1;

    double* bAwork;
    if (*storebP){ bAwork = new double[nToStore*(*nrandomP)*(*nclusterP)];          if (bAwork == NULL) *errP = 1; }
    else           bAwork = bM;
    int* rAwork;
    if (*storerP){ rAwork = new int[nToStore*(*nP)];                                if (rAwork == NULL) *errP = 1; }
    else           rAwork = rM;
    double* YsAwork;
    if (*storeyP){ YsAwork = new double[nToStore*(*nP)];                            if (YsAwork == NULL) *errP = 1; }
    else           YsAwork = YsM;
    double* uAwork;
    if (*storeuP){ uAwork = new double[nToStore*3*(*kmaxP)];                        if (uAwork == NULL) *errP = 1; }
    else           uAwork = uM;

    double* otherpAwork = new double[nToStore*1];                                                    if (otherpAwork == NULL) *errP = 1;
    double* mixMomentAwork = new double[nToStore*2];                                                 if (mixMomentAwork == NULL) *errP = 1;
    double* MHinfoAwork = new double[nToStore*nMHinfo];                                              if (MHinfoAwork == NULL) *errP = 1;
    double* MHinfo2Awork = new double[(*storeMHbP ? nToStore*nMHinfo2 : 1)];                         if (MHinfoAwork == NULL) *errP = 1;
    double* regresResAwork = new double[(*storeregresResP ? nToStore*(*nP) : 1)];                    if (regresResAwork == NULL) *errP = 1;

    if (*errP) 
      throw returnR("C++ Error: Could not allocate a memory for a working space, try to decrease 'nwrite'", *errP);


    // Create objects to work with beta parameters (fixed effects and means of random effects)
    // ========================================================================================
    // *betaMH ..... object to store fixed effects and simulation related information
    //
    MHblocks* betaMH = new MHblocks; 
    if (betaMH == NULL) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    *betaMH =  MHblocks(betaM, priorBetaI, priorBetaD, MHinfoM + 4, indbA, tolerCholP, iterM, &ONE_INT);

//    betaMH->print();
//    cout << endl;


    // Create data matrices and variables
    // ==================================
    // *clusteriA ........... indicators of clusters for observations
    // *invclusteriA ........ for each cluster: list of observations belonging to that cluster
    // *statusA ............. censoring status vector
    // *Y1A ................. lower limit of interval censored responses, value of left, right, exact responses
    // *Y2A ................. upper limit of interval censored responses
    // *indbinXA ............ columns in X that corresponds to a given random effect
    //                        = -1 for the random intercept
    // **ZZt ................ ZZt[j] = an array where a lower triangle of z_j*z_j' is stored, j = 1, ..., *nP
    // *diagIZZt ............ indeces of their diagonal elements
    // ***XXt ...............
    // **diagIXXt ...........
    //
    int* clusteriA = new int[*nP];                           if (clusteriA == NULL) *errP = 1;       
    List<int>* invclusteriA = new List<int>[*nclusterP];     if (invclusteriA == NULL) *errP = 1;
    int* statusA = new int[*nP];                             if (statusA == NULL) *errP = 1;
    double* Y1A;
    double* Y2A;
    int* indbinXA = new int[(*nrandomP ? *nrandomP : 1)];    if (indbinXA == NULL) *errP = 1;
    double** ZZt = new double*[*nP];                         if (ZZt == NULL) *errP = 1;            
    int* diagIZZt = new int[*nrandomP];                      if (diagIZZt == NULL) *errP = 1;
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);

    for (j = 0; j < *nP && !(*errP); j++){
      ZZt[j] = new double[*nrandomP * (*nrandomP)];          if (ZZt[j] == NULL) *errP = 1;
    }
    double*** XXt = new double**[betaMH->nBlock()];          if (XXt == NULL) *errP = 1;
    int** diagIXXt = new int*[betaMH->nBlock()];             if (diagIXXt == NULL) *errP = 1;
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);

    for (i = 0; i < betaMH->nBlock() && !(*errP); i++){
      l = betaMH->nInBlock[i];
      if (betaMH->typeUpd[i] == Gibbs){      
        XXt[i] = new double*[*nP];                             if (XXt[i] == NULL) *errP = 1;
        diagIXXt[i] = new int[l];                              if (diagIXXt[i] == NULL) *errP = 1;
        if (*errP) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
        for (j = 0; j < *nP && !(*errP); j++){
          XXt[i][j] = new double[l * l];                       if (XXt[i][j] == NULL) *errP = 1;
        }
      } 
    } 
    if (*errP) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);

    createData(nwithinA, clusteriA, invclusteriA, statusA, Y1A, Y2A, ZZt, diagIZZt, indbinXA, XXt, diagIXXt, 
               XA, YA, nP, nclusterP, nYP, nXP, nfixedP, nrandomP, randomIntP, indbA, 
               betaMH->nBlock(), betaMH->nInBlock, betaMH->indBlock, betaMH->typeUpd);


    // Create an object for a covariance matrix of random effects
    // ===========================================================
    // *Dcm ...... object to store a covariance matrix of random effects
    // nD ........ number of unique elements of the matrix D plus one (space for a determinant)
    covMatrix* Dcm = new covMatrix;  
    if (Dcm == NULL) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    *Dcm = covMatrix(DM, nrandomP, tolerCholP, tolerQRP);
    int nD = Dcm->larray() + 1;                                     

//    Dcm->print();
//    cout << endl;

    // Create an object for random effects and related quantities
    // ===========================================================
    // *bb .............. object to store random effects and some more information
    // *randomllcl ...... cluster contributions to the random log-likelihood
    // *proprandomllcl .. cluster contributions to the random log-likelihood used in the proposal
    //
    // *Eb0 ............. mean of the random intercept
    // *Eb0dependMix .... 1/0 indicating whether Eb0 depend on a mixture or not
    // NULA ............. helping to be pointed to
    //
    bblocks* bb = new bblocks;  
    double* randomllcl = new double[*nclusterP * (*nrandomP)];
    double* proprandomllcl = new double[*nclusterP * (*nrandomP)];
    const int Eb0dependMix[1] = {priorParI[2]};
    double NULA = 0.0;                      
    double* Eb0 = &NULA;                                       
    if (bb == NULL || randomllcl == NULL || proprandomllcl == NULL) 
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    *bb = bblocks(bM, priorbI, priorbD, MHinfo2M, tolerCholP);
    if (*Eb0dependMix) Eb0 = mixMomentM;
  
    //    bb->print();

    // Create matrices with initial estimates for mixture
    //  and matrices where proposals during reversible jumps will be stored
    // =====================================================================
    // *regresResM ....... regression residuals (y - x'beta - z'b)
    // *propregresResM ... 
    // *wM ............... current mixture weights 
    // *propwM ........... mixture weights used for the proposal during split-combine move
    // *muM .............. current mixture means 
    // *propmuM .......... mixture means used for the proposal during split-combine move
    // *invsigma2M ....... current mixture inverse-variances
    // *propinvsigma2M ... mixture inverse-variances used for the proposal during split-combine move
    // *proprM ........... component pertinences for the proposal during split-combine move
    // *invrM ............ matrices with indeces of observations belonging to each mixture component
    // *propinvrM ........
    // *mixtureNM ........ numbers of observations belonging to each mixture comp.
    // *propmixtureNM ....
    // *kP ............... current number of mixture components
    // *propkP ...........
    //
    double* regresResM = new double[*nP];                                                   
    double* propregresResM = new double[*nP];
    if (regresResM == NULL || propregresResM == NULL)
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    regresResidual(regresResM, YsM, betaM, bM, XA, clusteriA, randomIntP, indbA, nP, nXP, nrandomP);
    for (i = 0; i < *nP; i++) propregresResM[i] = regresResM[i];

    double* wM = new double[*kmaxP];
    double* muM = new double[*kmaxP];
    double* invsigma2M = new double[*kmaxP];
    double* propwM = new double[*kmaxP];
    double* propmuM = new double[*kmaxP];
    double* propinvsigma2M = new double[*kmaxP];
    int* proprM = new int[*nP];
    List<int>* invrM = new List<int>[*kmaxP];
    List<int>* propinvrM = new List<int>[*kmaxP];
    int* mixtureNM = new int[*kmaxP];
    int* propmixtureNM = new int[*kmaxP];
    int kP[1] = {int(mixtureM[0])};
    int propkP[1] = {*kP};
    if (wM == NULL || muM == NULL || invsigma2M == NULL || propwM == NULL || propmuM == NULL || propinvsigma2M == NULL ||
        proprM == NULL || invrM == NULL || propinvrM == NULL || mixtureNM == NULL || propmixtureNM == NULL)
           throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    createParam(nP, kmaxP, mixtureM, wM, muM, invsigma2M, rM, invrM, mixtureNM, 
 		                     propwM, propmuM, propinvsigma2M, proprM, propinvrM, propmixtureNM);
    mixMoments(mixMomentM, kP, wM, muM, invsigma2M);
    if (*nrandomP){
      randomLogLikelihood(loglikM + 1, randomllcl, &ZERO_INT, nclusterP, nclusterP, bM, betaM, Dcm, Eb0, indbinXA);
      for (j = 0; j < *nclusterP; j++) proprandomllcl[j] = randomllcl[j];
    }
//    if (*nrandomP) printArray(randomllcl, nclusterP);
//    cout << endl;


    // Objects to store log-likelihood contributions
    // ==============================================
    double* loglikobs = new double[*nP];
    double* proploglikobs = new double[*nP];
    if (loglikobs == NULL || proploglikobs == NULL) 
      throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    logLikelihood(loglikM + 0, loglikobs, nP, regresResM, YsM, kP, rM, wM, muM, invsigma2M, Eb0, errorTypeP, randomIntP, logdtrans);
    for (i = 0; i < *nP; i++) proploglikobs[i] = loglikobs[i];


    // Create matrices with some prior parameters and transform some prior hyperparameters   
    // ===================================================================================
    double* piSplitM = new double[*kmaxP + 1];
    double* logpiSplitM = new double[*kmaxP + 1];
    double* logpiCombineM = new double[*kmaxP + 1];

    double* piBirthM = new double[*kmaxP + 1];
    double* logpiBirthM = new double[*kmaxP + 1];
    double* logpiDeathM = new double[*kmaxP + 1];

    if (piSplitM == NULL || logpiSplitM == NULL || logpiCombineM == NULL || piBirthM == NULL || logpiBirthM == NULL || logpiDeathM == NULL)
       throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);

    createPriors(kmaxP, priorParD, piSplitM, logpiSplitM, logpiCombineM, piBirthM, logpiBirthM, logpiDeathM);


    // Further hyperparameters and their transformations 
    // ==================================================
    int* priorForkP = priorParI + 1;
    double* lambdaP =  priorParD + 2*(*kmaxP);
    double* deltaP =  priorParD + 2*(*kmaxP) + 1;
    double* xiP =  priorParD + 2*(*kmaxP) + 2;    
    double* kappaP =  priorParD + 2*(*kmaxP) + 3;    
    double* zetaP =  priorParD + 2*(*kmaxP) + 4;    
    double* gP =  priorParD + 2*(*kmaxP) + 5;    
    double* hP =  priorParD + 2*(*kmaxP) + 6;    
    double invkappaP[1] = {1/(*kappaP)};
    double xiInvkappaP[1] = {(*xiP)*(*invkappaP)};
    double sqrtkappaP[1] = {sqrt(*kappaP)};
    double halfl2pikappaP[1] = {LOG_SQRT_2PI + 0.5*log(*kappaP)};
    double* etaP = otherpM +0;                                   
    double lgammazetaP[1] = {lgammafn(*zetaP)};
    double llambdaP[1] = {log(*lambdaP)};

    if (*priorForkP == Fixed_k) for (i = 0; i < 4; i++) MHinfoM[i] = 1;

    // Functions and ways to generate a proposal vector for the split-combine and birth-death moves
    // ==============================================================================================
    // *uMtemp ........... pointer to the first component of uM which is used in a particular 
    //                     split/combine or birth/death move
    //                     with split and birth moves uMtemp[0], uMtemp[1], uMtemp[2] are used to propose new values
    //                     with combine and death moves uMtemp[-3], uMtemp[-2], uMtemp[-1] are filled by values of u
    //                       that would lead to a reversible move
    // *nuMtemp ....... length of uMtemp 
    //
    // logdu .......... log(density) of the canonical proposal vector u for split/combine move
    // ru ............. random generator of the canonical proposal vector u for split/combine move 
    //                  (used only during a development of this function )
    // *priorParmu .... parameters for a prior distribution of a canonical seed u
    //                  (used only during a development of this function)
    //
    // *mrcorrP ....... do we use a moody ring for correlated or uncorrelated variables
    // *mrepsilonP .... epsilon (time dependence) for a moody ring
    // *mrdeltaP ...... delta (between components dependence) for a moody ring
    // 
    // transuSC ....... transformation of the canonical proposal vector for split-combine move
    // invtransuSC .... inverse of the transformation of the canonical proposal vector for split-combine move
    // logJtransuSC ... logarithm of a Jacobian of a transformation of the canonical proposal vector for split-combine move
    // *transParmuSC .. parameters of a transformation from u -> v for split-combine move
    //
    // transuBD ....... transformation of the canonical proposal vector for birth-death move
    // invtransuBD .... inverse of the transformation of the canonical proposal vector for birth-death move
    // logJtransuBD ... logarithm of a Jacobian of a transformation of the canonical proposal vector for birth-death move
    // *transParmuBD .. parameters of a transformation from u -> v for birth-death move
    //
    double* uMtemp = uM + 3*(*kP);
    int nuMtemp[1] = {3*(*kmaxP -(*kP))};
    double (*logdu) (const double*, const double*);
    void (*ru) (double*, const double*);
    logdu = logdUnif;
    ru = rUnif;
    double priorParmu[1] = {0.0};                 // or whatever else, it effectively nowhere used

    int mrcorrP[1];
    double* mrepsilonP = revJumpParD + 0;         
    double* mrdeltaP = revJumpParD + 1;
    switch (revJumpParI[0]){
      case 0:     // basic algorithm
        *mrcorrP = 0;
        *mrepsilonP = 0.5;
        *mrdeltaP = 0.5;
	break;

      case 1:     // independent auxiliary variables
        *mrcorrP = 0;
        *mrdeltaP = 0.5;
        break;

      case 2:     // correlated auxiliary variables
        *mrcorrP = 1;
        break;

      default:
        throw returnR("C++ Error: Unknown algorithm to generate canonical proposal variables for reversible jumps.", 1);
    }

    void (*transuSC) (double*, const double*, const double*);
    void (*invtransuSC) (double*, const double*, const double*);
    double (*logJtransuSC) (const double*, const double*, const double*);
    double* transParmuSC = revJumpParD + 2;                                // next six components
    switch (revJumpParI[1]){
      case 0:      // Richardson-Green transformation
        transuSC = transBeBeBe;
        invtransuSC = invtransBeBeBe;
        logJtransuSC = logJtransBeBeBe;
        break;

      case 1:      // Brooks-Giudici-Roberts transformation
        transuSC = transBrooks;
        invtransuSC = invtransBrooks;
        logJtransuSC = logJtransBrooks;
        break;

      case 2:       // identity transformation
        transuSC = transId;
        invtransuSC = transId;        
        logJtransuSC = logJtransId;
        break;

      default:
        throw returnR("C++ Error: Unknown transformation to propose parameters during the split-combine move.", 1);
    }

    void (*transuBD) (double*, const double*, const double*);
    void (*invtransuBD) (double*, const double*, const double*);
    double (*logJtransuBD) (const double*, const double*, const double*);
    double* transParmuBD = revJumpParD + 8;              // next six components
    switch (revJumpParI[2]){
      case 0:       // Richardson-Green transformation
        transuBD = transBeNG;
        invtransuBD = invtransBeNG;
        logJtransuBD = logJtransBeNG;    
        transParmuBD[0] = 1.0;            transParmuBD[1] = *kP;       transParmuBD[2] = *xiP;  
        transParmuBD[3] = *sqrtkappaP;    transParmuBD[4] = *zetaP;    transParmuBD[5] = *etaP;  
        break;

      default:
        throw returnR("C++ Error: Unknown transformation to propose parameters during the birth-death move.", 1);
   
    }


    // ========== Main simulation ==========
      // nstored ........ number of values currently stored in *SimM matrices
      // write_flag ..... flag for writing to files 'a' = append
      // nullthIter ..... index for initial values (usually 0)
      // iter ........... counter of iterations
      // iterTotal ...... total number of iterations done previously (thinnig included)
      // iterTotalNow ... total number of iterations done with this function call
      //                  -> this number is used to compute acceptance rates
      // backs .......... how many times the carriage must be returned to print the next iteration number?
      // doadapt ........ it is set to 1 if witer = nthin
      //
    int nstored = 0;
    char write_flag = 'a';
    int nullthIter = *iterM;
    int iter;
    int iterTotal = nullthIter*nthin;
    int iterTotalNow = 0;
    int backs = 0;
    int doadapt = 1;

    Rprintf("Iteration ");
    for (iter=nullthIter + 1; iter <= nullthIter + niter; iter++){
      for (int witer = 1; witer <= nthin; witer++){  // thinning cycle
        iterTotal++;                                                // = (iter-1)*nthin + witer
        iterTotalNow++;                                             // = (iter-1)*nthin + witer - nullthIter*nthin
//        cout << " ==============================\n";
//        cout << "Iteration = " << iter << ",  within = " << witer << endl;

	updateData(YsM, regresResM, Y1A, Y2A, statusA, rM, wM, muM, invsigma2M, Eb0, kP, nP, errorTypeP, randomIntP);
        if ((*Eb0dependMix) && (*randomIntP))
          logLikelihood(loglikM + 0, loglikobs, nP, regresResM, YsM, kP, rM, wM, muM, invsigma2M, Eb0, errorTypeP, randomIntP, logdtrans);

  	updateWeights(&wM, &propwM, mixMomentM, 
                      loglikM + 0, &loglikobs, &proploglikobs,
                      loglikM + 1, &randomllcl, &proprandomllcl,
                      regresResM, YsM, betaM, bM, Dcm, kP, mixtureNM, muM, invsigma2M, rM, deltaP,
                      Eb0dependMix, randomIntP, nP, nclusterP, nrandomP, indbinXA, logdtrans); 

	updateMeans(muM, mixMomentM, regresResM, betaM, bM, Dcm, kP, mixtureNM, wM, invsigma2M, invrM, xiInvkappaP, invkappaP, 
                      Eb0dependMix, randomIntP, nP, nclusterP, nrandomP, indbinXA); 

	updateVars(invsigma2M, mixMomentM, Eb0, regresResM, kP, mixtureNM, wM, muM, rM, zetaP, etaP, randomIntP, nP);

	updateEta(etaP, kP, invsigma2M, zetaP, gP, hP);

	updateAlloc(rM, invrM, mixtureNM, wM, muM, invsigma2M, kP, regresResM, Eb0, randomIntP, nP);

        logLikelihood(loglikM + 0, loglikobs, nP, regresResM, YsM, kP, rM, wM, muM, invsigma2M, Eb0, errorTypeP, randomIntP, logdtrans);
        if ((*Eb0dependMix) && (*randomIntP))
	  randomLogLikelihood(loglikM + 1, randomllcl, &ZERO_INT, nclusterP, nclusterP, bM, betaM, Dcm, Eb0, indbinXA);

        if (*priorForkP != Fixed_k){
  	  if (*kP < *kmaxP)
            moodyRing(uMtemp, uM + 0, mrepsilonP, mrdeltaP, nuMtemp, mrcorrP, &ZERO_INT);

	  birthDeath(MHinfoM + 2, MHinfoM + 3, kP, 
                     loglikM + 0, &loglikobs, &proploglikobs,
                     loglikM + 1, &randomllcl, &proprandomllcl,
		     wM, muM, invsigma2M, mixMomentM, rM, invrM, mixtureNM,
	             propkP,
                     uMtemp, logdu, transuBD, invtransuBD, logJtransuBD,
                     regresResM, YsM, bM, betaM, Dcm, 
	             kmaxP, piBirthM, logpiBirthM, logpiDeathM,
	             deltaP, xiP, invkappaP, sqrtkappaP, halfl2pikappaP, zetaP, etaP, lgammazetaP, llambdaP,
                     priorParmu, transParmuBD, priorForkP, Eb0dependMix, randomIntP, 
                     indbinXA, nP, nclusterP, logdtrans);
          uMtemp = uM + 3*(*kP);
          *nuMtemp = 3*(*kmaxP -(*kP));

	  for (i = 0; i < *kP; i++){
            if (!R_finite(wM[i]) || !R_finite(muM[i]) || !R_finite(invsigma2M[i])){
              REprintf("\nk=%d", *kP);
              REprintf(";  weights="); for (j=0; j < *kP; j++) REprintf("%e, ", wM[j]);
              REprintf(";  means="); for (j=0; j < *kP; j++) REprintf("%e, ", muM[j]);
              REprintf(";  invvars="); for (j=0; j < *kP; j++) REprintf("%e, ", invsigma2M[j]);
              REprintf("\n");
              for(i = 0; i < backs; i++) Rprintf("\b");
              Rprintf("%d", iter - 1);
              writeToFiles(iterAwork, loglikAwork, mixtureAwork, mixMomentAwork, betaAwork, bAwork, DAwork, rAwork, YsAwork, 
                           otherpAwork, uAwork, MHinfoAwork, MHinfo2Awork, regresResAwork,
                           nstored, dir, write_flag, ncmSM, nMHinfo, nMHinfo2, nD, kmaxP, nXP, nrandomP, nclusterP, nP, 
                           storebP, storeyP, storerP, storeuP, storeMHbP, storeregresResP);  
              nstored = 0;
              throw returnR("Trap after birth-death move, NaN appears in the definition of the mixture", 1);
            }
          }

	  splitCombine(MHinfoM + 0, MHinfoM + 1,
	               loglikM + 0, &loglikobs, &proploglikobs,
	               kP, &wM, &muM, &invsigma2M, rM, &invrM, &mixtureNM,
	               propkP, &propwM, &propmuM, &propinvsigma2M, proprM, &propinvrM, &propmixtureNM,
	               uMtemp, logdu, transuSC, invtransuSC, logJtransuSC,
	               regresResM, YsM, Eb0,
	               kmaxP, randomIntP, piSplitM, logpiSplitM, logpiCombineM,
	               deltaP, xiP, invkappaP, halfl2pikappaP, zetaP, etaP, lgammazetaP, llambdaP,
                       priorParmu, transParmuSC, priorForkP, logdtrans, nP);
          uMtemp = uM + 3*(*kP);
          *nuMtemp = 3*(*kmaxP -(*kP));

	  for (i = 0; i < *kP; i++){
            if (!R_finite(wM[i]) || !R_finite(muM[i]) || !R_finite(invsigma2M[i])){
              REprintf("\nk=%d", *kP);
              REprintf(";  weights="); for (j=0; j < *kP; j++) REprintf("%e, ", wM[j]);
              REprintf(";  means="); for (j=0; j < *kP; j++) REprintf("%e, ", muM[j]);
              REprintf(";  invvars="); for (j=0; j < *kP; j++) REprintf("%e, ", invsigma2M[j]);
              REprintf("\n");
              for(i = 0; i < backs; i++) Rprintf("\b");
              Rprintf("%d", iter - 1);
              writeToFiles(iterAwork, loglikAwork, mixtureAwork, mixMomentAwork, betaAwork, bAwork, DAwork, rAwork, YsAwork, 
                           otherpAwork, uAwork, MHinfoAwork, MHinfo2Awork, regresResAwork,
                           nstored, dir, write_flag, ncmSM, nMHinfo, nMHinfo2, nD, kmaxP, nXP, nrandomP, nclusterP, nP, 
                           storebP, storeyP, storerP, storeuP, storeMHbP, storeregresResP);  
              nstored = 0;
              throw returnR("Trap after split-combine move, NaN appears in the definition of the mixture", 1);
            }
          }
        }

        if (*nrandomP){ 
          updateRandom(bb, regresResM, propregresResM, 
                       loglikM + 1, randomllcl, proprandomllcl,                               
                       loglikM + 0, loglikobs, proploglikobs,
                       betaM, Eb0, Dcm, YsM, XA, ZZt, invclusteriA, indbinXA, logdtrans, randomIntP, nXP, nP, 
                       errorTypeP, kP, rM, wM, muM, invsigma2M, tolerCholP);

	  updateCovMatRandom(Dcm, bM, betaM, Eb0, bb->priorD, bb->dfD, bb->scaleD, indbinXA, nclusterP, nP, tolerCholP, tolerQRP);
	  randomLogLikelihood(loglikM + 1, randomllcl, &ZERO_INT, nclusterP, nclusterP, bM, betaM, Dcm, Eb0, indbinXA);
        }

 	if (*nXP){ 
          doadapt = 1*(witer == nthin);
          updateRegres(betaMH, &regresResM, &propregresResM,
                       loglikM + 0, &loglikobs, &proploglikobs,
                       loglikM + 1, &randomllcl, &proprandomllcl,
                       YsM, Eb0, bM, Dcm, XA, XXt, indbA, indbinXA, logdtrans, &iter, &doadapt, randomIntP, nrandomP, nP, nclusterP,
                       errorTypeP, kP, rM, wM, muM, invsigma2M, tolerCholP);
        } 
      }    // end of the thinning cycle

      // Store values in arrays
      storeInArrays(iterAwork, loglikAwork, mixtureAwork, mixMomentAwork, betaAwork, bAwork, DAwork, rAwork, YsAwork, 
                    otherpAwork, uAwork, MHinfoAwork, MHinfo2Awork, regresResAwork,
                    iter, loglikM, kP, wM, muM, invsigma2M, mixMomentM, betaM, bM, DM, Dcm->det, rM, YsM, 
                    etaP, uM, MHinfoM, MHinfo2M, regresResM,
                    nstored, iterTotalNow, kmaxP, ncmSM, nMHinfo, nMHinfo2, nD, nXP, nrandomP, nclusterP, nP, betaMH->nBlock(),
                    storebP, storeyP, storerP, storeuP, storeMHbP, storeregresResP);  

      nstored++;


      // Write to files
      if (!(iter % nwrite)){
        for (i = 0; i < backs; i++) Rprintf("\b");
        Rprintf("%d", iter);
        backs = int(log10(double(iter))) + 1;
        writeToFiles(iterAwork, loglikAwork, mixtureAwork, mixMomentAwork, betaAwork, bAwork, DAwork, rAwork, YsAwork, 
                     otherpAwork, uAwork, MHinfoAwork, MHinfo2Awork, regresResAwork,
                     nstored, dir, write_flag, ncmSM, nMHinfo, nMHinfo2, nD, kmaxP, nXP, nrandomP, nclusterP, nP, 
                     storebP, storeyP, storerP, storeuP, storeMHbP, storeregresResP);  
        write_flag = 'a'; 
        nstored = 0;
      } 
    }    // end of the simulation

    // Write not yet written values from the last supcycle
    if (nstored){
      for(i = 0; i < backs; i++) Rprintf("\b");
      Rprintf("%d", iter - 1);
      writeToFiles(iterAwork, loglikAwork, mixtureAwork, mixMomentAwork, betaAwork, bAwork, DAwork, rAwork, YsAwork, 
                   otherpAwork, uAwork, MHinfoAwork, MHinfo2Awork, regresResAwork,
                   nstored, dir, write_flag, ncmSM, nMHinfo, nMHinfo2, nD, kmaxP, nXP, nrandomP, nclusterP, nP, 
                   storebP, storeyP, storerP, storeuP, storeMHbP, storeregresResP);  
      nstored = 0;
    }
    Rprintf("\n");

    // Do some adjustments on last sampled values before returning them back to R
    *iterM = iter - 1;
    for (i = 0; i < *nP; i++) rM[i]++;        // C++ -> R
    mixtureM[0] = double(*kP);
    for (j = 0; j < *kP; j++){
      mixtureM[1 + j] = wM[j];
      mixtureM[1 + (*kmaxP) + j] = muM[j];
      mixtureM[1 + 2*(*kmaxP) + j] = 1/invsigma2M[j];
    }
    for (j = *kP; j < *kmaxP; j++){
      mixtureM[1 + j] = 0.0;
      mixtureM[1 + (*kmaxP) + j] = 0.0;
      mixtureM[1 + 2*(*kmaxP) + j] = 0.0;
    }
    
    PutRNGstate();

    // Cleaning
    delete [] iterAwork;       delete [] loglikAwork;    delete [] mixtureAwork;
    delete [] betaAwork;       delete [] DAwork;         delete [] otherpAwork;
    delete [] MHinfoAwork;     delete [] MHinfo2Awork;
    delete [] regresResAwork;  delete [] mixMomentAwork;

    if (*storebP) delete [] bAwork;
    if (*storerP) delete [] rAwork;
    if (*storeyP) delete [] YsAwork;
    if (*storeuP) delete [] uAwork;    

    delete [] clusteriA;       delete [] invclusteriA;
    delete [] statusA;         delete [] indbinXA;

    for (j = 0; j < *nP; j++) delete [] ZZt[j];          
    delete [] ZZt;                                       delete [] diagIZZt;
    for (i = 0; i < betaMH->nBlock(); i++){
      if (betaMH->typeUpd[i] == Gibbs){
        for (j = 0; j < *nP; j++) delete [] XXt[i][j];
        delete [] XXt[i];
        delete [] diagIXXt[i];
      }
    }
    delete [] XXt;                                       delete [] diagIXXt;    

    delete betaMH;             delete Dcm;               delete bb;

    delete [] randomllcl;      delete [] proprandomllcl;
    delete [] loglikobs;       delete [] proploglikobs;
    delete [] regresResM;      delete [] propregresResM;
    delete [] wM;              delete [] propwM;
    delete [] muM;             delete [] propmuM;        
    delete [] invsigma2M;      delete [] propinvsigma2M;  
                               delete [] proprM;         
    delete [] mixtureNM;       delete [] propmixtureNM;   
    delete [] invrM;           delete [] propinvrM;

    delete [] piSplitM;        delete [] logpiSplitM;    delete [] logpiCombineM;
    delete [] piBirthM;        delete [] logpiBirthM;    delete [] logpiDeathM;
  
    delete [] MHinfoM;         delete [] MHinfo2M;

    return;

  }  // end of try
  catch(returnR rr){
    *errP = rr.errflag();
    PutRNGstate();
    return;
  }

}   // end of function bayesAFT1

}     // end of extern "C"



