// Header for functions used to run Bayesian survival regression models

// 03/11/2003: start working on it
// ========================================================================

#ifndef BAYESSURVREG_H
#define BAYESSURVREG_H

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <iomanip>
#include <climits>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>

#include "AK_Error.h"
#include "List.h"
#include "templatefun.h"
#include "classes.h"


// Types of updates
enum updateTypes {RandomWalk, AdaptiveM, Gibbs};

// Types of error distributions
enum errorTypes {Mixture, Spline, PolyaTree, WhoKnows};

// Types of priors for k
enum priorForkTypes {Poisson, Uniform, Fixed};

// Types of priors for matrix D
enum priorForDTypes {InvWishart, SDUniform};

const double SPI = 2.5066282746310007;                          // sqrt(2*pi)
const double ROOT_2 = 1.414213562373095048801688724210;         // sqrt(2)
const double LOG_SQRT_2PI = 0.918938533204672741780329736406;	// log(sqrt(2*pi)) 
const double LOG_2 = 0.69314718055994529;                       // log(2)
const double ZERO = 1e-50;
const double invFLT_MAX = 1/FLT_MAX;                                  

const double SCALE_ZERO = 1e-10;                                // zero scale (used in updateEta, updateVars)
const double NORM_ZERO = 1e-16;                                 // qnorm(1 - 1e-16) is still non infty in R
//const double QNORM_ONE = 10.0;                                  // qnorm(1), instead of infty (used in updateData)
//const double QNORM_ZERO = -10.0;                                // qnorm(1), instead of -infty (used in updateData)
                                                                  // approximately -40 = qnorm(1e-323)
                                                                  //               8.2 = qnorm(1 - 1e-16)
                                                                  //  these are last values that do not produce inf/-inf in R
const int ONE_INT = 1;
const int ZERO_INT = 0;
const double ZERO_DOUBLE = 0.0;
//const int AKINT_MAX = INT_MAX;
const int AKINT_MAX = 100;


// =========================================================
// Log of the Jacobian of the transformation
//   * for different transformations
// =========================================================
inline double
logdlogtrans(const double y){ return -y;}

inline double
logdidenttrans(const double y){ return 0;}


// ==================================================================================
// ***** arrays2Mat.cpp: Functions to create various matrices from input arrays *****
// ==================================================================================
void
create2pointer(double**& point,    const int lpoint,  const int llayer2,  
               const double init,  const bool doinit);

void
create3pointer(double***& point,   const int lpoint,  const int* llayer2,  const int llayer3,  
               const double init,  const bool doinit);

void
create3pointer(double***& point,   const int lpoint,  const int llayer2,  const int* llayer3,  
               const double init,  const bool doinit);

void
createDataShort(int* nwithinA,          int* clusteriA,            List<int>* invclusteriA,
                const double* XA,
                double** ZZt,           int* diagIZZt,             int* indbinXA,     
                const int* nP,          const int* nclusterP,
                const int* nXP,         const int* nfixedP,        const int* nrandomP,
                const int* randomIntP,  const int* indbA);

void 
createData(int* nwithinA,          int* clusteriA,            List<int>* invclusteriA, 
           int* statusA,           double*& Y1A,              double*& Y2A,           
           double** ZZt,           int* diagIZZt,             int* indbinXA,     
           double*** XXt,          int** diagIXXt,      
           const double* XA,       double* YA,
           const int* nP,          const int* nclusterP,      const int* nYP,
           const int* nXP,         const int* nfixedP,        const int* nrandomP,
           const int* randomIntP,  const int* indbA,      
           const int& nBlockbeta,  const int* nInBlockbeta,   int** indBlockbeta,        const int* typeUpdbeta);

void 
createParam(const int* nP,                 const int* kmaxP,             const double* mixtureA,        
            double* wM,                    double* muM,                  double* invsigma2M,
            int* rM,                       List<int>* invrM,             int* mixtureNM,  
            double* propwM,                double* propmuM,              double* propinvsigma2M,
            int* proprM,                   List<int>* propinvrM,         int* propmixtureNM);

void
createPriors(const int *kmaxP,     const double* priorParD,
             double* piSplitM,       
             double* logpiSplitM,  double* logpiCombineM,
             double* piBirthM,       
             double* logpiBirthM,  double* logpiDeathM);


// =======================================================
// ***** miscellaneous.cpp: Helping functions *****
// =======================================================
void
giveMixtureN(int* mixtureNM, 
             const int* kP,   const int* rM,  const int* nP);

void
giveMixtureN(int* mixtureNM, 
             const int* kP,  const List<int>* invrM);

void
Y2T(double** T,  const double* Y,  const int iter,  const int* nP,  double (*itrans)(const double));

void
giveSigmaAndInvsigma2(double* sigmaM,  double* invsigma2M,  const double* sigma2M,  const int* kP);

void
regresResidual(double* regresResA,
               const double* YsA,    const double* betaA,     const double* bA,
               const double* XA,     const int* clusteriA,    const int* randomIntP,   
               const int* indbA,     const int* nP,           const int* nXP,         const int* nrandomP);

void 
regresPredictor(double* regresPredA,
                const double* betaA,  const double* bA,
                const double* XA,     const int* clusteriA,    const int* randomIntP,   
                const int* indbA,     const int* nP,           const int* nXP,         const int* nrandomP);

void
regresResidual(double* regresResA,
               const double* YsA,   const double* newYsA,
               const int* nP);

void
regresResidual(double* regresResA,
               const double* betaA,   const double* newbetaA,   const int* indnewA,  
               const int* nnewP,      const double* XA,         const int* indbA,         
               const int* nP);

void
regresResidual(double* regresResA,
               const double* bA,      const double* newbA,  const int* indnewA,
               const int* nnewP,      const double* XA,     const int* clusteriA,    const int* randomIntP,   
               const int* indbinXA,   const int* nP,        const int* nXP,          const int* nrandomP);

void
regresResidual(double* regresResA,
               const double* bA,         const double* bclA,   const int* cl,
               const List<int>* indobs,  const double* XA,     const int* randomIntP,
               const int* indbinXA,      const int* nP,        const int* nXP,         const int* nrandomP);

void
regresResidual(double* regresResA,
               const double* bA,         const double* bclA,   
               const int* indnewA,       const int* nnewP,     const int* cl,
               const List<int>* indobs,  const double* XA,     const int* randomIntP,
               const int* indbinXA,      const int* nP,        const int* nXP,         const int* nrandomP);

void 
storeInArrays(int* iterA,              double* loglikA,            double* mixtureA,           double* mixMomentA,
              double* betaA,           double* bA,                 double* DA,             
              int* rA,                 double* YSA,                double* otherpA,
              double* uA,              double* MHinfoA,            double* MHinfo2A,           double* regresResA,
              const int& iterindex,    const double* logLikP,      const int* kP,      
              const double* wM,        const double* muM,          const double* invsigma2M,   const double* mixMomentM,
              const double* betaM,     const double* bM,           const double* DM,           const double& Ddet,
              const int* rM,           const double* YsM,          const double* etaP,
              const double* uM,        const int* MHinfoM,         const int* MHinfo2M,        const double* regresResM,
              const int& nstored,      const int& iterNow,
              const int* kmaxP,        const int& ncmSM,           const int& nMHinfo,         const int& nMHinfo2,       const int& nD,
              const int* nXP,          const int* nrandomP,        const int* nclusterP,       const int* nP,
              const int& nBetaBlocks,
              const int* storebP,      const int* storeyP,         const int* storerP,         const int* storeuP,    
              const int* storeMHbP,    const int* storeregresResP);

void
writeToFiles(const int* iterA,      const double* loglikA,       const double* mixtureA,  const double* mixMomentA,
             const double* betaA,   const double* bA,            const double* DA,
             const int* rA,         const double* YSA,           const double* otherpA,
             const double* uA,      const double* MHinfoA,       const double* MHinfo2A,  const double* regresResA,
             const int& nstored,    const std::string& dir,      const char& write_flag, 
             const int& ncmSM,      const int& nMHinfo,          const int& nMHinfo2,
             const int& nD,         const int *kmaxP,
             const int* nXP,        const int* nrandomP,         const int* nclusterP,    const int* nP,
             const int* storebP,    const int* storeyP,          const int* storerP,      const int* storeuP,   
             const int* storeMHbP,  const int* storeregresResP);

void
writeToFiles2(double** ETsA,          double** TsA,           
              double*** SurvA,        double*** hazardA,        double*** cumhazardA,
              const int& nstored,     const std::string& dir,   const char& write_flag, 
              const int *nP,          const int* ngridM,  
              const int* storeETP,    const int* storeTP,     
              const int* storeSurvP,  const int* storehazardP,   const int* storecumhazardP);

void
writeToFiles3(double** ETquant,         double** Tquant,       
              double*** Survquant,      double*** hazardquant,      double*** cumhazardquant,
              const int& nquant,        const std::string& dir,     const char& write_flag, 
              const int *nP,            const int* ngridM,  
              const int* predictETP,    const int* predictTP,  
              const int* predictSurvP,  const int* predicthazardP,  const int* predictcumhazardP);

void
readMixtureFromFiles(double* array,               const int nR,   const int kmax,
                     const std::string& dir,  
                     const std::string& kname,    
                     const std::string& wname,    
                     const std::string& muname,  
                     const std::string& sigma2name,  
                     const int skip);

void
writeMixtureToFiles(const double* mixtureA,   const int nR,                   const int kmax,
            const std::string& dir,           const std::string& wname,  
            const std::string& muname,        const std::string& sigma2name,
            const char &flag,
	    const int& prec = 6,              const int& width = 0);

std::string
giveString(const int num);


// =======================================================
// ***** logLikelihood.cpp:   *****
// =======================================================
void
logLikelihood(double* loglikelhood,       double* loglikobs,
              const int* nP,
              const double* regresResM,   const double* YsM,
              const int* kP,              const int* rM,                
              const double* wM,           const double* muM,      const double* invsigma2M,
              const double* Eb0,          const int* errorTypeP,  const int* randomIntP,
              double (*logdtrans)(const double));

void
logLikelihood(double* loglikelhood,       double* loglikobs,
              const List<int>* obsUpd,    const int* nP,
              const double* regresResM,   const double* YsM,
              const int* kP,              const int* rM,                
              const double* wM,           const double* muM,      const double* invsigma2M,
              const double* Eb0,          const int* errorTypeP,  const int* randomIntP,
              double (*logdtrans)(const double));

void
clusterlogLikelihood(double* clusterloglik,  const double* loglikobs,  const int* cl,  const List<int>* obsInCluster);


// ========================================================
// ***** randomLogLikelihood.cpp:    *****
// ========================================================
void
randomLogLikelihood(double* randomloglik,    double* randomllcl,
                    const int* clusterUpd,   const int* nclusterUpd,  const int* nclusterP,
                    const double* bM,        const double* betaM,     const covMatrix* Dcm,
                    const double* Eb0,       const int* indbinXA);

void
randomLogLikelihood(double* randomloglik,   double* randomllcl,
                    const int* clusterUpd,  const int* nclusterP, 
                    const double* bclM,     const double* betaM,    const covMatrix* Dcm,
                    const double* Eb0,      const int* indbinXA);


// =======================================================
// ***** mixMoments.cpp: *****
// =======================================================
void
mixMoments(double* mixMomentM,  const int* kP,             const double* wM,
           const double* muM,   const double* invsigma2M,  const bool onlySD = false);

void
mixMean(double* mixMeanP,  const int* kP,  const double* wM,  const double* muM);


// =======================================================
// ***** updateData.cpp: *****
// =======================================================
void
updateData(double* YsM,             double* regresResM,
           const double* Y1M,       const double* Y2M,     const int* statusM,
           int* rM,                 double* cumwM,         const double* muM,     const double* invsigma2M,
           const double* Eb0,       const int* kP,         const int* nP,        
           const int* errorTypeP,   const int* randomIntP);


// =======================================================
// ***** updateWeights.cpp: *****
// =======================================================
void
updateWeights(double** wM,               double** propwM,            double* mixMomentM,
              double* loglikelhood,      double** loglikobs,         double** proploglikobs,
              double* randomloglik,      double** randomllcl,        double** proprandomllcl,
              const double* regresResM,  const double* YsM,          const double* betaM,        
              const double* bM,          const covMatrix* Dcm,
              const int* kP,             const int* mixtureNM,  
              const double* muM,         const double* invsigma2M,   const int* rM,
              const double* deltaP,
              const int* Eb0dependMix,   const int* randomIntP,      
              const int* nP,             const int* nclusterP,       const int* nrandomP,
              const int* indbinXA,       double (*logdtrans)(double));


// =======================================================
// ***** updateMeans.cpp: *****
// =======================================================
void
updateMeans(double* muM,                double* mixMomentM,
            const double* regresResM,   const double* betaM,        const double* bM,
            const covMatrix* Dcm,
            const int* kP,              const int* mixtureNM,
            const double* wM,           const double* invsigma2M,   const List<int>* invrM,  
            const double* xiInvkappaP,  const double* invkappaP,  
            const int* Eb0dependMix,    const int* randomIntP,      
            const int* nP,              const int* nclusterP,       const int* nrandomP,
            const int* indbinXA);


// =======================================================
// ***** updateVars.cpp: *****
// =======================================================
void
updateVars(double* invsigma2M,        double* mixMomentM,      double* Eb0,
           const double* regresResM,  
           const int* kP,             const int* mixtureNM,
           const double* wM,          const double* muM,       const int* rM, 
           const double* zetaP,       const double* etaP,  
           const int* randomIntP,     const int* nP);


// =======================================================
// ***** updateEta.cpp: *****
// =======================================================
void
updateEta(double* etaP,
          const int* kP,        const double* invsigma2M,
          const double* zetaP,  const double* gP,           const double* hP);


// =======================================================
// ***** updateAlloc.cpp: *****
// =======================================================
void
updateAlloc(int* rM,                List<int>* invrM,          int* mixtureNM, 
            const double* wM,       const double* muM,         const double* invsigma2M,
            const int* kP,          const double* regresResM,  const double* Eb0,
            const int* randomIntP,  const int* nP);

void
updateAlloc(int* rM,                
            const double* wM,       const double* muM,         const double* invsigma2M,
            const int* kP,          const double* regresResM,  const double* Eb0,
            const int* randomIntP,  const int* nP);


// =======================================================
// ***** birthDeath.cpp: *****
// =======================================================
void
birthDeath(int* acceptedP,             int* birthP,                   int* kP, 
           double* loglikelhood,       double** loglikobs,            double** proploglikobs,
           double* randomloglik,       double** randomllcl,           double** proprandomllcl,
           double* wM,                 double* muM,                   double* invsigma2M,
           double* mixMomentM,
           int* rM,                    List<int>* invrM,              int* mixtureNM, 
           int* propkP,                          
           double *uM,
           double (*logdu) (const double*, const double*),
           void (*transu) (double*, const double*, const double*),
           void (*invtransu) (double*, const double*, const double*),
           double (*logJtransu) (const double*, const double*, const double*),
           const double* regresResM,   const double* YsM,
           const double* bM,           const double* betaM,           const covMatrix* Dcm,
           const int* kmaxP,              
           const double* piBirthM,     const double* logpiBirthM,     const double* logpiDeathM,
           const double* deltaP,       const double* xiP,             const double* invkappaP, 
           const double* sqrtkappaP,   const double* halfl2pikappaP,  const double* zetaP,
           const double* etaP,         const double* lgammazetaP,     const double* llambdaP,
           const double* priorParmu,   double* transParmu,          
           const int* priorForkP,      const int* Eb0dependMix,       const int* randomIntP,
           const int* indbinXA,        const int* nP,                 const int* nclusterP,
           double (*logdtrans) (double));

void
proposeDeath(int& jdeath,        double* vM,
             const int& numE,    const int* emptyCompM,
             const double* wM,   const double* muM,        const double* invsigma2M);

double
logdtransBirthDeath(const double* vM,          const double* uM,  
                    const double* priorParmu,  const double* transParmu,          
                    double (*logdu) (const double*, const double*),
                    double (*logJtransu) (const double*, const double*, const double*),
                    const bool* RichardsonGreenP);

double
logPostRatioJacobianBirthDeath(const int* shortkP,       const double* vM,
                               const int* nP,
                               const double* deltaP,
                               const double* xiP,        const double* invkappaP,  const double* halfl2pikappaP,
                               const double* zetaP,      const double* etaP,       const double* lgammazetaP,
                               const double* llambdaP,   const int* priorForkP,    const bool* RichardsonGreenP);

void
moveParamsBirthDeath(int& jdeath,                   
                     double* wM,          double* muM,        double* invsigma2M,
                     int* rM,             List<int>* invrM,   int* mixtureNM,
                     const int* propkP,   const double* vM,   const int* birthP);

int
numEmpty(int* emptyCompM, const int *kP,  const int* mixtureNM);


// =======================================================
// ***** splitCombine.cpp: *****
// =======================================================
void
splitCombine(int* acceptedP,                int* splitP,                 
             double* loglikelhood,          double** loglikobs,          double** proploglikobs,
             int* kP,                             
             double** wM,                   double** muM,                double** invsigma2M,
             int* rM,                       List<int>** invrM,           int** mixtureNM,
             int* propkP,                         
             double** propwM,               double** propmuM,            double** propinvsigma2M,
             int* proprM,                   List<int>** propinvrM,       int** propmixtureNM,
             double* uM,              
             double (*logdu) (const double*, const double*),
             void (*transu) (double*, const double*, const double*),
             void (*invtransu) (double*, const double*, const double*),
             double (*logJtransu) (const double*, const double*, const double*),
             const double* regresResM,      const double* YsM,            const double* Eb0,
             const int* kmaxP,              const int* randomIntP,
             const double* piSplitM,        const double* logpiSplitM,    const double* logpiCombineM,
             const double* deltaP,          const double* xiP,            const double* invkappaP,                     
             const double* halfl2pikappaP,  const double* zetaP,          const double* etaP,
             const double* lgammazetaP,     const double* llambdaP,       
             const double* priorParmu,      const double* transParmu,
             const int* priorForkP,         double (*logdtrans) (double), const int* nP);

void
proposeSplit(int* acceptedP,
             double* propwM,      double* propmuM,     double* propinvsigma2M,
             const double* wM,    const double* muM,   const double* invsigma2M,
             const double* vM,    const int jsplit,    const int* kP);

void
proposeCombine(int* acceptedP,
               double* vM,
               double* propwM,    double* propmuM,     double* propinvsigma2M,
               const double* wM,  const double* muM,   const double* invsigma2M,
               const int jsplit,  const int* propkP);

double
allocSplit(int* proprM,               List<int>* propinvrM,    int* propmixtureNM,
           const int* rM,             const List<int>* invrM,  const  int* mixtureNM,
           const double* propwM,      const double* propmuM,   const double* propinvsigma2M,
           const int jsplit,          const int* kP,           
           const double* regresResM,  const double* Eb0,       const int* randomIntP);

double
allocCombine(int* proprM,               List<int>* propinvrM,    int* propmixtureNM,
             const int* rM,             const List<int>* invrM,  const int* mixtureNM,
             const double* wM,          const double* muM,       const double* invsigma2M,
             const int jsplit,          const int* propkP,
             const double* regresResM,  const double* Eb0,       const int* randomIntP);

double 
logPostRatioSplitCombine(const int jsplit,              const int* shortkP,
                         const double* longwM,          const double* shortwM,     
                         const double* longmuM,         const double* shortmuM,
                         const double* longinvsigma2M,  const double* shortinvsigma2M,
                         const int* longmixtureNM,      const int* shortmixtureNM,
                         const double* deltaP,          const double* xiP,              const double* invkappaP,                   
                         const double* halfl2pikappaP,  const double* zetaP,            const double* etaP,                        
                         const double* lgammazetaP,     const double* llambdaP,         const int* priorForkP);

double
logJacobianSplitCombine(const double w,
                        const double mu0,           const double mu1,
                        const double invsigma20,    const double invsigma21,
                        const double invsigma2,     const double* vM);


// =======================================================
// ***** updateRegres.cpp: *****
// =======================================================
void
updateRegres(MHblocks* betaMH,             double** regresResM,   double** propregresResM,   
             double* loglikelhood,         double** loglikobs,    double** proploglikobs,
             double* randomloglik,         double** randomllcl,   double** proprandomllcl,
             const double* YsM,            const double* Eb0,     const double* bM,      
             const covMatrix* Dcm,         const double* XA,      double*** XXt, 
             const int* indbA,             const int* indbinXA,
             double (*logdtrans)(double),  const int* iterTotal,  const int* doadapt,
             const int* randomIntP,        const int* nrandomP,
             const int* nP,                const int* nclusterP,
             const int* errorTypeP,        const int* kP,         const int* rM,
             const double* wM,             const double* muM,     const double* invsigma2M,
             const double* tolerChol);

void
GIBBSproposalFixed(double* betaM,            double* regresResM,
                   double* covpar,           double* ichicovpar,
                   const double* priormean,  const double* priorinvvar,
                   const double* Eb0,        const int* randomIntP,
                   const double* XA,         double** XXtb,              const int* diagIXXtb,
                   const double* wM,         const double* muM,          const double* invsigma2M,  const int* rM,
                   const int* indUpd,        const int* npar,            const int* nupdate,        const int* nP,
                   const int* errorTypeP,    const double* tolerChol);

void
GIBBSproposalMeanRandom(double* betaM,
                        double* covpar,           double* ichicovpar,
                        const double* priormean,  const double* priorinvvar,
                        const double* Eb0,        const double* bM,
                        const covMatrix* Dcm,     const int* diagIXXtb,
                        const int* indbA,         const int* indbinXA,
                        const int* indUpd,        const int* npar,            
                        const int* nupdate,       const int* nrandomP,         const int* nclusterP,
                        const int* errorTypeP,    const double* tolerChol);


// =======================================================
// ***** updateRandom.cpp: *****
// =======================================================
void
updateRandom(bblocks* bb,                   double* regresResM,    double* propregresResM,
            double* randomloglik,           double* randomllcl,    double* proprandomllcl,
            double* loglikelhood,           double* loglikobs,     double* proploglikobs,
            const double* betaM,            const double* Eb0,     const covMatrix* Dcm,
            const double* YsM,              const double* XA,      double** ZZt,
            const List<int>* invclusteriA,  const int* indbinXA,   double (*logdtrans)(double),
            const int* randomIntP,          const int* nXP,        const int* nP,
            const int* errorTypeP,          const int* kP,         const int* rM,
            const double* wM,               const double* muM,     const double* invsigma2M,
            const double* tolerChol);

void
GIBBSproposalRandom(double* bM,                     double* regresResM, 
                    double* covpar,                 double* ichicovpar,
                    const double* betaM,            const double* Eb0,     const covMatrix* Dcm,
                    const double* XA,               double** ZZt,        
                    const double* wM,               const double* muM,     const double* invsigma2M,  const int* rM,
                    const int* randomIntP,          const int* nrandomP,       
                    const int* nXP,                 const int* nclusterP,  const int* nP,             const int* errorTypeP,
                    const List<int>* invclusteriA,  const int* indbinXA,   const int* indUpd,
                    const double* tolerChol);


// =======================================================
// ***** MHproposal.cpp: *****
// =======================================================
void
AMproposal(double* proppar,              double* chcovpar,     
           const double* covpar,         const double* par,          const int* indUpd,
           const int* npar,              const int* nupdate,         const int* diagI,             
           const double* halfRangeUnif,  const double* weightUnif,   
           const double* eps,            const double* sdNum,        const double* tolerChol);

void
AMadapt(double* covpar,     double* meanpar,     const double* par,
        const int* indUpd,  const int* nupdate,  const int* diagI,   const int* iter,
        const double* eps,  const double* sdNum);

void
MHproposal(double* proppar,              
           const double* chcovpar,       const double* par,          const int* indUpd,
           const int* npar,              const int* nupdate,         const int* diagI,             
           const double* halfRangeUnif,  const double* weightUnif);


// =======================================================
// ***** updateCovMatRandom.cpp: *****
// =======================================================
void
updateCovMatRandom(covMatrix* Dcm,           const double* bM,             
                   const double* betaM,      const double* Eb0,
                   const int* priorD,        const double* priordf,   const double* priorscaleMat, 
                   const int* indbinXA,      const int* nclusterP,    const int* nP,      
                   const double* tolerChol,  const double* tolerQR);

void
sumSquares(double* sumSq,          const double* bM,       
           const double* betaM,    const double* Eb0,
           const int* indbinXA,    const int* diagI,       
           const int* nclusterP,   const int* nrandomP,    const int* lSS);


// =======================================================
// ***** predictSurv.cpp: *****
// =======================================================
void
predictSurv(double*** SM,              double*** lambdaM,         double*** LambdaM,
            const int iter,            double** gridM,            double** loggridM,  const double* time0P,
            const double* regresPredM,
            const int* rM,             const double* wM,          const double* muM,  const double* sigmaM,
            const double* Eb0,         const int* kP,             const int* nP,      const int* ngridM,
            const int* errorTypeP,     const int* randomIntP,
            const int* hazardP,        const int* cumhazardP);

void
predictData(double* YsM,             const double* regresPredM,
            int* rM,                 double* cumwM,             const double* muM,     const double* sigmaM,
            const double* Eb0,       const int* kP,             const int* nP,        
            const int* errorTypeP,   const int* randomIntP);

void
predictRandom(double* bM,  
              const double* betaM,     const double* Eb0,    const covMatrix* Dcm,
              const int* nrandomP,     const int* nclusterP,  
              const int* indbinXA,     const int* indUpd);

void
predictET(double** ET,            const double* time0P, const int iter,
          const double* betaM,    const double* wM,     const double*muM,       const double* sigma2M,
          const covMatrix* Dcm,   const double* XA,
          const int* kP,          const int* nP,        const int* nXP,         const int* indbinXA,
          const int* randomIntP,  const int* nrandomP,  const int* errorTypeP);

void
Y2T(double** T,  const double* Y,  const double* time0P, const int iter,  const int* nP,  double (*itrans)(const double));


// =======================================================
// ***** quantile.cpp: *****
// =======================================================
void
cumsumQuantile1(double** quant,  double** newval,  const int nquant,  const int nobs,  
                const int iter);

void
cumsumQuantile2(double*** quant,  double*** newval,     const int nquant,  const int nobs,  const int* ngridM,  
                const int iter);

void
meanQuantile1(double** quant,           double** sample,   
              const double* quantileA,  const int* indquant1,  const int* indquant2,
              const int nobs,           const int nquant,      const int sampleSize);

void
meanQuantile2(double*** quant,          double*** sample,  
              const double* quantileA,  const int* indquant1,  const int* indquant2,
              const int nobs,           const int* ngridM,     const int nquant,      const int sampleSize);


extern "C"{

// =================================================================================
// ***** random.cpp: set of routines for random sampling and related functions *****
// =================================================================================
void
rltruncGamma(double* x,  
             const double* shape,  const double* rate,   const double* minx,
             const int* n,         const int* callFromR);

void
discreteUniformSampler(int* sampledj,
                       const int* kP,  const int* nP,  const int* callFromR);

void 
discreteSampler(int* sampledj, 
                double* propA,  const int* kP,  const int* nP,  const int* cumul,  const int* callFromR);

int 
findUniformIndex(const double u,  const int startInd,  const int endInd,  const int k);

int 
findIndex(const double u,  const int startInd,  const int endInd,  const double* ValuesA);


// =======================================================
// ***** moodyRing.cpp: Routines for a moody ring *****
// =======================================================
void
moodyRing(double* uA,              double* moodP,             
          const double* timeDepP,  const double* componentDepP,     
          const int* nuP,          const int* corrP,
          const int* callFromR);

void
corr_moodyRing(double* uA,              double* moodA,                double* initmoodP,
               const double* timeDepP,  const double* componentDepP,  const int* nuP,
               const int* nP,           const int* callFromR);

void
indep_moodyRing(double* uA,               double* inituA,
                const double* timeDepP,   const int* nuP,
                const int* nP,            const int* callFromR);


// ======================================================================================
// ***** cholesky.cpp: Routines for a Cholesky decomposition and related functions *****
// ======================================================================================
void
cholesky(double* C,  int* rankC,  const int* nC,  const int* diagI, const double* toler);

void
cholesky2(double* C,  int* rankC,  const int* nC,  const double* toler);

void 
chinv(double* C , const int* nC,  const int* diagI, const int* onlyCholInv);

void 
chinv2(double* C , double* ichol, const int* nC,  const int* diagI);

void
chposDef(const double* C,      double* cholC,      int* rankC,           
         int* attempt,         const int* nC,      const int* diagI,  
         const double* toler,  const double* eps,  const int* nattempt = &AKINT_MAX);


// ======================================================================================
// ***** qrdecomp.cpp: Routines for a QR decomposition of a matrix
//                    - grabbed from dqrdc2.f and blas.f
// ======================================================================================
void 
dqrdc2CPP(double* x, const int* n, const int* p, const double* tol, int* k, double* qraux, int* jpvt);

double
ddotCPP(const int n, double* dx, const int incx, double* dy, const int incy);

void
daxpyCPP(const int n, const double da, double* dx, const int incx, double* dy, const int incy);

void
dscalCPP(const int n, const double da, double* dx, const int incx);

double 
dnrm2CPP(const int n, const double* x, const int incx);


// ======================================================================================
// ***** mvtdist.cpp: Routines for dealing with the multivariate distributions
// ======================================================================================
void
LxMxtL(double* LAL,     const double* L,   const double* A,
       const int* nA,   const int* diagI);

void
tLxMxL(double* LAL,     const double* L,   const double* A,
       const int* nA,   const int* diagI);

void
Mxa(double* Ma,     const double* a,  const double* A,  const int* inda,  
    const int* na,  const int* nA,    const int* diagI);

void
Mxa2(double* Ma,     const double* a,  const double* A,  const int* indA,  
     const int* na,  const int* nA,    const int* diagI);

void
Wxa(double* Wa,   const double* a,   const double* A,  const int* indr,  const int* indc,
    const int* na,  const int* nA,  const int* nrow,   const int* diagI);

void
axMxa(double* aMa,    const double* a,  const double* A,  const int* inda,  
      const int* na,  const int* nA,    const int* diagI);

void
dmvtnorm(double* dens,      const double* x,      const double* mean,  const double* vari,  
         const int* indx,   const int* indxrepl, 
         const int* nx,     const int* nmean,     const int* nxrepl,    
         const int* nP,     const int* diagI,     const int* logP);

void
rmvtnorm(double* x,        const double* mean,   const double* L,  
         const int* indx,  const int* indxrepl,
         const int* nx,    const int* nmean,     const int* nxrepl,      
         const int* nP,    const int* diagI,     const int* callFromR);

void
rmvtnorm2(double* x,        const double* mean,   const double* iLi,  
          const int* indx,  const int* indxrepl,
          const int* nx,    const int* nmean,     const int* nxrepl,      
          const int* nP,    const int* diagI,     const int* callFromR);

void
rmvtiunif(double* x,        const double* mean,   const double* halfRange,  
          const int* indx,  const int* indxrepl,  
          const int* nx,    const int* nmean,     const int* nxrepl,
          const int* nP,    const int* callFromR);

void
rwishart(double* w,
         const int* p,      const double* nu,   const double* L,
         const int* diagI,  const int* nP,      const int* callFromR);

void
rinvwishart(double* w,
            const int* p,      const double* nu,   const double* L,
            const int* diagI,  const int* nP,      const int* callFromR);

void
rwishart2(double* w,
          const int* p,      const double* nu,   const double* iLi,
          const int* diagI,  const int* nP,      const int* callFromR);


// =======================================================
// ***** propVector.cpp: *****
// =======================================================
double
logdUnif(const double* u,  const double* priorParmu);

void
rUnif(double* u,  const double* priorParmu);


void
transId(double* v, const double* u, const double* transParmu);

double
logJtransId(const double* u, const double* v, const double* transParmu);


void
transBeBeBe(double* v,  const double* u,  const double* transParmu);

void
invtransBeBeBe(double* u,  const double* v,  const double* transParmu);

double
logJtransBeBeBe(const double* u, const double* v, const double* transParmu);


void
transBrooks(double* v, const double* u, const double* transParmu);

void
invtransBrooks(double* u, const double* v, const double* transParmu);

double
logJtransBrooks(const double* u, const double* v, const double* transParmu);

void
transBeNG(double* v, const double* u, const double* transParmu);

void
invtransBeNG(double* u, const double* v, const double* transParmu);

double
logJtransBeNG(const double* u, const double* v, const double* transParmu);


// =======================================================
// ***** bayesDensity.cpp: *****
// =======================================================
void
bayesDensity(double* aver,      double* staver,      double* intercept,     double* scale,      int* Mk,
             char** dirP,       const double* grid,  const double* stgrid,
             const int* kmax,   const int* M,        const int* skip,
             const int* ngrid,  const int* nstgrid,
             const int* type,   int* errP);

}   // end of extern "C"

#endif    
