// =============================================================
// Header for in_output.cpp: Input/output functions       *****
// =============================================================
#ifndef _IN_OUTPUT_H_
#define _IN_OUTPUT_H_

#include <R.h>

#include <iostream>
#include <iomanip>
#include <string>

#include "AK_Error.h"
#include "templatefun.h"

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
readMixtureFromFiles(double* array,                  int* nread,
                     const int& nR,                  const int& skip,     const int& by,     const int& kmax,
                     const std::string& dir,  
                     const std::string& kname,    
                     const std::string& wname,    
                     const std::string& muname,  
                     const std::string& sigma2name);

void
writeMixtureToFiles(const double* mixtureA,   const int nR,                   const int kmax,
            const std::string& dir,           const std::string& wname,  
            const std::string& muname,        const std::string& sigma2name,
            const char &flag,
	    const int& prec = 6,              const int& width = 0);

std::string
giveString(const int num);

#endif
