// Smaller helping function for the Bayesian AFT model

// 24/11/2003: start working on it
// 06/08/2004: function 'readMixtureFromFiles' added
// 08/08/2004: function 'writeMixtureToFiles' added

#include "bayessurvreg.h"

using namespace std;

// ================================================================================
// Function to compute numbers of observations belonging to each mixture component
// ** overloaded function with two prototypes
// ================================================================================
void
giveMixtureN(int* mixtureNM,                   // result (at least *kP x 1)
             const int* kP,                    // number of mixture components
             const int* rM,                    // vector of pertainence indicators  (*nP x 1)
             const int* nP                     // number of observations             
             )
{
  for (int j = 0; j < *kP; j++) mixtureNM[j] = 0;
  for (int obs = 0; obs < *nP; obs++){
      mixtureNM[rM[obs]]++;
  }

  return;
}   // end of function giveMixtureN


// ================================================================================
void
giveMixtureN(int* mixtureNM,           // result (at least *kP x 1)
             const int* kP,            // number of mixture components
             const List<int>* invrM    // pointers to lists indicating observations belonging to different mixture components
             )
{
  for (int j = 0; j < *kP; j++) mixtureNM[j] = invrM[j].length();

  return;
}   // end of function giveMixtureN


// ================================================================================
// invsigma2 -> sigma
// ================================================================================
void
giveSigmaAndInvsigma2(double* sigmaM,  double* invsigma2M,  const double* sigma2M,  const int* kP)
{
  int j;

  for (j = 0; j < *kP; j++){
    if (sigma2M[j] <= 0.0){
      sigmaM[j] = 0.0;
      invsigma2M[j] = FLT_MAX;
    }
    else{
      sigmaM[j] = sqrt(sigma2M[j]);
      invsigma2M[j] = 1 / sigma2M[j];
    }
  }
  
  return;
}


// ================================================================================
// Function to compute regression residuals
//  r = y - x'beta - z'b
//
// ** overloaded function with six prototypes
//
// ================================================================================
void
regresResidual(double* regresResA,
               const double* YsA,    const double* betaA,     const double* bA,
               const double* XA,     const int* clusteriA,    const int* randomIntP,   
               const int* indbA,     const int* nP,           const int* nXP,         const int* nrandomP)
{
  int i, j;

  for (i = 0; i < *nP; i++){
    regresResA[i] = YsA[i];
    if (*randomIntP) regresResA[i] -= bA[(*nrandomP)*clusteriA[i]];     // subtract a random intercept
    for (j = 0; j < *nXP; j++){
      if (indbA[j] == -1)
        regresResA[i] -= XA[(*nP)*j + i] * betaA[j];
      else
        regresResA[i] -= XA[(*nP)*j + i] * bA[(*nrandomP)*clusteriA[i] + indbA[j]];
    }
  }

  return;
}


// ================================================================================
// ***** regresPredictor *****
//
// Function to compute (x'beta + z'b)
//
// ================================================================================
void 
regresPredictor(double* regresPredA,
                const double* betaA,  const double* bA,
                const double* XA,     const int* clusteriA,    const int* randomIntP,   
                const int* indbA,     const int* nP,           const int* nXP,         const int* nrandomP)
{
  int i, j;

  for (i = 0; i < *nP; i++){
    regresPredA[i] = 0.0;
    if (*randomIntP) regresPredA[i] += bA[(*nrandomP)*clusteriA[i]];     // subtract a random intercept
    for (j = 0; j < *nXP; j++){
      if (indbA[j] == -1)
        regresPredA[i] += XA[(*nP)*j + i] * betaA[j];
      else
        regresPredA[i] += XA[(*nP)*j + i] * bA[(*nrandomP)*clusteriA[i] + indbA[j]];
    }
  }

  return;
}



// ================================================================================
// Function to update regression residuals when y changes
//  r = y - x'beta - z'b
// ================================================================================
//
// YsA ...... old value of Y
// newYsA ... new value of Y
//
void
regresResidual(double* regresResA,
               const double* YsA,   const double* newYsA,
               const int* nP)
{
  int i;

  for (i = 0; i < *nP; i++){
    regresResA[i] -= (YsA[i] - newYsA[i]);
  }

  return;
}


// ================================================================================
// Function to update regression residuals when beta changes
// NOTE: changes in beta's that correspond to random effects do not have any impact
//       to regression residuals
//  r = y - x'beta - z'b
// ================================================================================
//
// betaA ...... array of ALL old betas
// newbetaA ... array of ALL new betas ONLY                     
// indnewA .... indeces of new betas in the whole beta vector   length = nnewP
//
void
regresResidual(double* regresResA,
               const double* betaA,   const double* newbetaA,   const int* indnewA,  
               const int* nnewP,      const double* XA,         const int* indbA,         
               const int* nP)
{
  int i, j;

  for (i = 0; i < *nP; i++){
    for (j = 0; j < *nnewP; j++){
      if (indbA[indnewA[j]] == -1) regresResA[i] += XA[(*nP)*indnewA[j] + i] * (betaA[indnewA[j]] - newbetaA[indnewA[j]]);
    }
  }

  return;
}


// ================================================================================
// Function to update regression residuals when b (random effects) changes
//  * change of same random effects for all clusters
//  r = y - x'beta - z'b
// ================================================================================
//
// bA ........ array of ALL old values of b                 length = (nclusterP * nrandomP)
// newbA ..... array of NEW values of b ONLY                length = (nclusterP * nnewP)
// indnewA ... indeces of new b's in the whole sequence
//
void
regresResidual(double* regresResA,
               const double* bA,      const double* newbA,  const int* indnewA,
               const int* nnewP,      const double* XA,     const int* clusteriA,    const int* randomIntP,   
               const int* indbinXA,   const int* nP,        const int* nXP,          const int* nrandomP)
{
  int i, j;

  int start = 0;
  for (i = 0; i < *nP; i++){
    if (*randomIntP && (indnewA[0] == 0)){ 
      start = 1;
      regresResA[i] += (bA[(*nrandomP)*clusteriA[i]] - newbA[(*nnewP)*clusteriA[i]]);
    }
    for (j = start; j < *nnewP; j++){
      regresResA[i] += XA[(*nP)*indbinXA[indnewA[j]] + i] * (bA[(*nrandomP)*clusteriA[i] + indnewA[j]] - newbA[(*nnewP)*clusteriA[i] + j]);
    }
  }

  return;
}


// ===================================================================================================
// Function to update regression residuals when all b (random effects) changes for a specific cluster
//  r = y - x'beta - z'b
// ===================================================================================================
//
// bA ....... array of ALL old values of b                          (nrandomP x nclusterP)
// bclA ..... array of new values of b for a specific cluster only  (nrandomP)
// cl ....... index of a cluster for which b changed                (1)
// indobs ... indeces of observations belonging to a cluster cl     
//
void
regresResidual(double* regresResA,
               const double* bA,         const double* bclA,   const int* cl,
               const List<int>* indobs,  const double* XA,     const int* randomIntP,
               const int* indbinXA,      const int* nP,        const int* nXP,         const int* nrandomP)
{
  int nobs = indobs->length();
  int i, j, obs;

  for (i = 0; i < nobs; i++){
    obs = (*indobs)[i];
    if (*randomIntP) 
      regresResA[obs] += (bA[(*nrandomP)*(*cl)] - bclA[0]);
    for (j = *randomIntP; j < *nrandomP; j++){
      regresResA[obs] += XA[(*nP)*indbinXA[j] + obs] * (bA[(*nrandomP)*(*cl) + j] - bclA[j]);
    }
  }

  return;
}

// ===================================================================================================
// Function to update regression residuals when some b (random effects) changes for a specific cluster
//  r = y - x'beta - z'b
// ===================================================================================================
//
// bA ....... array of ALL old values of b                             (nrandomP x nclusterP)
// bclA ..... array of BOTH new values of b and unchanged values of b  
//            for a specifiv cluster                                   (nrandomP)
// indnewA .. indeces of changed b's                                   (nnewP)
// nnewP .... number of changed b's                                    (1)
// cl ....... index of a cluster for which b changed                   (1)
// indobs ... indeces of observations belonging to a cluster cl     
//
void
regresResidual(double* regresResA,
               const double* bA,         const double* bclA,   
               const int* indnewA,       const int* nnewP,     const int* cl,
               const List<int>* indobs,  const double* XA,     const int* randomIntP,
               const int* indbinXA,      const int* nP,        const int* nXP,         const int* nrandomP)
{
  int nobs = indobs->length();
  int i, j, obs;

  int start = 0;
  for (i = 0; i < nobs; i++){
    obs = (*indobs)[i];
    if (*randomIntP && (indnewA[0] == 0)){ 
      start = 1;
      regresResA[obs] += (bA[(*nrandomP)*(*cl)] - bclA[0]);
    }
    for (j = start; j < *nnewP; j++){
      regresResA[obs] += XA[(*nP)*indbinXA[indnewA[j]] + obs] * (bA[(*nrandomP)*(*cl) + indnewA[j]] - bclA[indnewA[j]]);
    }
  }

  return;
}



// ==================================================================================================
// ****** storeInArrays: store simulated values in arrays,
//                       version for bayessurvreg1
// ==================================================================================================
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
              const int* storeMHbP,    const int* storeregresResP)
{
  int i, j, l, obs;

  iterA[nstored] = iterindex;

  loglikA[nstored*2] = logLikP[0];
  loglikA[nstored*2 + 1] = logLikP[1];
  
  mixtureA[nstored*(ncmSM) + 0] = *kP;
  for (j = 0; j < *kP; j++){
    mixtureA[nstored*(ncmSM) + 1 + j] = wM[j];
    mixtureA[nstored*(ncmSM) + 1 + (*kmaxP) + j] = muM[j];           
    mixtureA[nstored*(ncmSM) + 1 + 2*(*kmaxP) + j] = 1/invsigma2M[j];
  }
  //  for (j = 0; j < *kmaxP - (*kP); j++){
  //    mixtureA[nstored*(ncmSM) + 1 + (*kP) + j] = FLT_MAX;
  //    mixtureA[nstored*(ncmSM) + 1 + (*kmaxP) + (*kP) + j] = FLT_MAX;           
  //    mixtureA[nstored*(ncmSM) + 1 + 2*(*kmaxP) + (*kP) + j] = FLT_MAX;
  //  }

  mixMomentA[nstored*2] = mixMomentM[0];
  mixMomentA[nstored*2 + 1] = mixMomentM[1];

  if (*nXP) for (i = 0; i < *nXP; i++) betaA[nstored*(*nXP) + i] = betaM[i];

  if (*storebP) for (i = 0; i < *nclusterP; i++)
    for (j = 0; j < *nrandomP; j++) bA[nstored*(*nrandomP)*(*nclusterP) + (*nrandomP)*i + j] = bM[(*nrandomP)*i + j];

  if (*nrandomP){
    DA[nstored*nD] = Ddet;
    for (i = 0; i < nD - 1; i++) DA[nstored*nD + i + 1] = DM[i];
  }

  if (*storeyP) for (obs = 0; obs < *nP; obs++) YSA[nstored*(*nP) + obs] = YsM[obs];

  if (*storerP) for (obs = 0; obs < *nP; obs++) rA[nstored*(*nP) + obs] = rM[obs] + 1;         // give R index
  otherpA[nstored] = *etaP;

  if (*storeuP) for (i = 0; i < 3*(*kmaxP); i++) uA[nstored*3*(*kmaxP) + i] = uM[i];

  for (i = 0; i < nMHinfo; i++){
    MHinfoA[nstored * nMHinfo + i] = double(MHinfoM[i])/double(iterNow);
  }

  if (*storeMHbP){
    for (i = 0; i < nMHinfo2; i++){
      MHinfo2A[nstored * nMHinfo2 + i] = double(MHinfo2M[i])/double(iterNow);
    }
  }

  if (*storeregresResP){
    for (obs = 0; obs < *nP; obs++) regresResA[nstored*(*nP) + obs] = regresResM[obs];
  }

  return;
}


// ==================================================================================================
// ****** writeToFiles: write simulated values to files
//                      version for bayessurvreg1
// ==================================================================================================
// temporary???: give here fixed file names
void
writeToFiles(const int* iterA,      const double* loglikA,       const double* mixtureA,  const double* mixMomentA,
             const double* betaA,   const double* bA,            const double* DA,
             const int* rA,         const double* YSA,           const double* otherpA,
             const double* uA,      const double* MHinfoA,       const double* MHinfo2A,  const double* regresResA,
             const int& nstored,    const std::string& dir,      const char& write_flag, 
             const int& ncmSM,      const int& nMHinfo,          const int& nMHinfo2,
             const int& nD,         const int* kmaxP,
             const int* nXP,        const int* nrandomP,         const int* nclusterP,    const int* nP,
             const int* storebP,    const int* storeyP,          const int* storerP,      const int* storeuP,   
             const int* storeMHbP,  const int* storeregresResP)
{
  int kmax = (ncmSM - 1)/3;
  writeToFile(iterA, nstored, 1, dir, "/iteration.sim", write_flag);
  writeToFile(loglikA, nstored, 2, dir, "/loglik.sim", write_flag);
  writeMixtureToFiles(mixtureA, nstored, kmax, dir, "/mweight.sim", "/mmean.sim", "/mvariance.sim", write_flag);
  writeTwoToFile(mixtureA, nstored, ncmSM, 0, mixMomentA, nstored, 2, dir, "/mixmoment.sim", write_flag);
  if (*nXP) writeToFile(betaA, nstored, *nXP, dir, "/beta.sim", write_flag);

  if (*nrandomP){
    writeToFile(DA, nstored, nD, dir, "/D.sim", write_flag);     

    if (*storebP) writeToFile(bA, nstored, *nrandomP*(*nclusterP), dir, "/b.sim", write_flag);
    else          writeToFile(bA, 1, *nrandomP*(*nclusterP), dir, "/b.sim", 'o');
  }

  if (*storeyP) writeToFile(YSA, nstored, *nP, dir, "/Y.sim", write_flag);             
  else          writeToFile(YSA, 1, *nP, dir, "/Y.sim", 'o');             

  if (*storerP) writeToFile(rA, nstored, *nP, dir, "/r.sim", write_flag);   // here: R index was computed already by 'storeInArrays' function         
  else{
    int* rAtemp = new int[*nP];      if (rAtemp == NULL) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    for (int obs = 0; obs < *nP; obs++) rAtemp[obs] = rA[obs] + 1;         // give R index
    writeToFile(rAtemp, 1, *nP, dir, "/r.sim", 'o');
    delete [] rAtemp;
  }

  writeToFile(otherpA, nstored, 1, dir, "/otherp.sim", write_flag);
  writeToFile(MHinfoA, nstored, nMHinfo, dir, "/MHinfo.sim", write_flag);
  if (*storeMHbP) writeToFile(MHinfo2A, nstored, nMHinfo2, dir, "/MHbinfo.sim", write_flag);

  if (*storeuP) writeToFile(uA, nstored, 3*(*kmaxP), dir, "/u.sim", write_flag);
  else          writeToFile(uA, 1, 3*(*kmaxP), dir, "/u.sim", 'o');

  if (*storeregresResP) writeToFile(regresResA, nstored, *nP, dir, "/regresres.sim", write_flag);

  return;
}  // end of function writeToFiles


// ==================================================================================================
// ****** writeToFiles2: write simulated values to files
//                       version for 'predictive'
// ==================================================================================================
void
writeToFiles2(double** ETsA,          double** TsA,           
              double*** SurvA,        double*** hazardA,        double*** cumhazardA,
              const int& nstored,     const std::string& dir,   const char& write_flag, 
              const int *nP,          const int* ngridM,  
              const int* storeETP,    const int* storeTP,     
              const int* storeSurvP,  const int* storehazardP,   const int* storecumhazardP)
{
  string file;
  int obs;

  if (*storeETP) writeToFile2(ETsA, *nP, nstored, dir, "/predET.sim", write_flag);
  if (*storeTP) writeToFile2(TsA, *nP, nstored, dir, "/predT.sim", write_flag);
  if (*storeSurvP){
    for (obs = 0; obs < *nP; obs++){
      file = string("/predS") + giveString(obs + 1);
      writeToFile2(SurvA[obs], ngridM[obs], nstored, dir, file, write_flag);
    }
  }
  if (*storehazardP){
    for (obs = 0; obs < *nP; obs++){
      file = string("/predhazard") + giveString(obs + 1);
      writeToFile2(hazardA[obs], ngridM[obs], nstored, dir, file, write_flag);
    }
  }
  if (*storecumhazardP){
    for (obs = 0; obs < *nP; obs++){
      file = string("/predcumhazard") + giveString(obs + 1);
      writeToFile2(cumhazardA[obs], ngridM[obs], nstored, dir, file, write_flag);
    }
  }
  return;
}  


// ==================================================================================================
// ****** writeToFiles3: write computed quantiles and predictive means to files
//                       version for 'predictive'
// ==================================================================================================
void
writeToFiles3(double** ETquant,         double** Tquant,       
              double*** Survquant,      double*** hazardquant,      double*** cumhazardquant,
              const int& nquant,        const std::string& dir,     const char& write_flag, 
              const int *nP,            const int* ngridM,  
              const int* predictETP,    const int* predictTP,  
              const int* predictSurvP,  const int* predicthazardP,  const int* predictcumhazardP)
{
  string file;
  int obs, i;

  if (*predictETP) writeToFile2(ETquant, *nP, nquant + 1, dir, "/quantET.sim", write_flag);
  if (*predictTP) writeToFile2(Tquant, *nP, nquant + 1, dir, "/quantT.sim", write_flag);
  if (*predictSurvP){
    for (obs = 0; obs < *nP; obs++){
      file = string("/quantS") + giveString(obs + 1);
      writeToFile2(Survquant[obs], ngridM[obs], nquant + 1, dir, file, write_flag);
    }
  }
  if (*predicthazardP){
    for (obs = 0; obs < *nP; obs++){
      file = string("/quanthazard") + giveString(obs + 1);
      writeToFile2(hazardquant[obs], ngridM[obs], nquant + 1, dir, file, write_flag);
    }
  }
  if (*predictcumhazardP){
    for (obs = 0; obs < *nP; obs++){
      file = string("/quantcumhazard") + giveString(obs + 1);
      writeToFile2(cumhazardquant[obs], ngridM[obs], nquant + 1, dir, file, write_flag);
    }
  }

  return;
}


// =======================================================================
// ***** readMixtureFromFiles: Function to read mixture from four files 
// =======================================================================
// The file is read by ROWS and written to an array
//
// array........ array where to read it        (length = nR * (1+3*kmax))
// nR .......... number of rows in files THAT ARE TO BE READ
// dir ......... directory where all files are stored
// kname ....... name of the file where 'k' is stored in the first column
// wname ....... name of the file where mixture weights are stored (1 mixture/row)
// muname ...... name of the file where mixture means are stored (1 mixture/row)
// sigma2name .. name of the file where mixture variances are stored (1 mixture/row)
// skip ........ number of rows that are to be skipped at the beginnig of each file
//
void
readMixtureFromFiles(double* array,               const int nR,   const int kmax,
                     const std::string& dir,  
                     const std::string& kname,    
                     const std::string& wname,    
                     const std::string& muname,  
                     const std::string& sigma2name,  
                     const int skip)
{
  try{
    int size = nR * kmax;
    if (size <= 0){
      throw returnR("C++ Error: Files of null size are to be read.", 99);
    }

    std::string kpath = dir + kname;
    std::string wpath = dir + wname;
    std::string mupath = dir + muname;
    std::string sigma2path = dir + sigma2name;

    std::string errmes, mess;
    std::ifstream kfile(kpath.c_str(), std::ios::in);    
    std::ifstream wfile(wpath.c_str(), std::ios::in);    
    std::ifstream mufile(mupath.c_str(), std::ios::in);    
    std::ifstream sigma2file(sigma2path.c_str(), std::ios::in);    

    if (!kfile){
      errmes = std::string("C++ Error: Could not open ") + kpath;
      throw returnR(errmes, 99);
    } 
    if (!wfile){
      errmes = std::string("C++ Error: Could not open ") + wpath;
      throw returnR(errmes, 99);
    } 
    if (!mufile){
      errmes = std::string("C++ Error: Could not open ") + mupath;
      throw returnR(errmes, 99);
    } 
    if (!sigma2file){
      errmes = std::string("C++ Error: Could not open ") + sigma2path;
      throw returnR(errmes, 99);
    } 

    int kcurrent;
    char ch;
    Rprintf("Reading mixture files. \n");

    // Skip rows at the beginning that are to be skipped
    for (int i = 0; i < skip; i++){
      kfile.get(ch);        while (ch != '\n') kfile.get(ch);
      wfile.get(ch);        while (ch != '\n') wfile.get(ch);
      mufile.get(ch);       while (ch != '\n') mufile.get(ch);
      sigma2file.get(ch);   while (ch != '\n') sigma2file.get(ch);
    }

    // Read values that are to be read
    for (int i = 0; i < nR; i++){
      if (kfile.eof()){
        errmes = std::string("C++ Error: Reached end of file ") + kpath + " before "
                 + char(nR) + std::string(" values were read.");
        throw returnR(errmes, 99);
      }
      kfile >> kcurrent;
      if (kcurrent > kmax) throw returnR("C++ Error: k value higher than kmax was read.", 99);
      array[i*(1+3*kmax)+0] = kcurrent;
      kfile.get(ch);                     while (ch != '\n') kfile.get(ch);

      if (wfile.eof()){
        errmes = std::string("C++ Error: Reached end of file ") + wpath + " before "
                 + char(nR) + std::string(" sets of mixture weights were read.");
        throw returnR(errmes, 99);
      }
      for (int j = 0; j < kcurrent; j++) wfile >> array[i*(1+3*kmax) + 1 + j];
      wfile.get(ch);                 while (ch != '\n') wfile.get(ch);

      if (mufile.eof()){
        errmes = std::string("C++ Error: Reached end of file ") + mupath + " before "
                 + char(nR) + std::string(" sets of mixture means were read.");
        throw returnR(errmes, 99);
      }
      for (int j = 0; j < kcurrent; j++) mufile >> array[i*(1+3*kmax) + 1 + kmax + j];
      mufile.get(ch);                while (ch != '\n') mufile.get(ch);

      if (mufile.eof()){
        errmes = std::string("C++ Error: Reached end of file ") + sigma2path + " before "
                 + char(nR) + std::string(" sets of mixture variances were read.");
        throw returnR(errmes, 99);
      }
      for (int j = 0; j < kcurrent; j++) sigma2file >> array[i*(1+3*kmax) + 1 + 2*kmax + j];
      sigma2file.get(ch);            while (ch != '\n') sigma2file.get(ch);
    }

    kfile.close();
    wfile.close();
    mufile.close();
    sigma2file.close();
    return;
  }  // end of try
  catch(returnR){
    throw;
  }  
}    // end of the function readMixtureFromFiles


// ========================================================================================
// ***** writeMixtureToFiles: Write mixture to files (except number of mixture components)
// ========================================================================================
//
// mixtureA ........ array where the mixture is stored
//
void
writeMixtureToFiles(const double* mixtureA,   const int nR,                   const int kmax,
            const std::string& dir,           const std::string& wname,  
            const std::string& muname,        const std::string& sigma2name,
            const char &flag,
            const int& prec,                  const int& width)
{
  try{
    std::string wpath = dir + wname;
    std::string mupath = dir + muname;
    std::string sigma2path = dir + sigma2name;

    std::ofstream wout;
    std::ofstream muout;
    std::ofstream sigma2out;
    bool err = false;

    std::string errmess;
    if (flag == 'n') {
      std::fstream wtemp(wpath.c_str(), std::ios::in);
      std::fstream mutemp(mupath.c_str(), std::ios::in);
      std::fstream sigma2temp(sigma2path.c_str(), std::ios::in);
      if (!wtemp) {
	wout.open(wpath.c_str(), std::ios::out);
      } 
      else {
	wtemp.close();
	err = true;
      }
      if (!mutemp) {
	muout.open(mupath.c_str(), std::ios::out);
      } 
      else {
	mutemp.close();
	err = true;
      }
      if (!sigma2temp) {
	sigma2out.open(sigma2path.c_str(), std::ios::out);
      } 
      else {
	sigma2temp.close();
	err = true;
      }
    }
    else if (flag == 'o'){
           wout.open(wpath.c_str(), std::ios::out | std::ios::trunc);
           muout.open(mupath.c_str(), std::ios::out | std::ios::trunc);
           sigma2out.open(sigma2path.c_str(), std::ios::out | std::ios::trunc);
         }
         else if (flag == 'a'){
                wout.open(wpath.c_str(), std::ios::out | std::ios::app);
                muout.open(mupath.c_str(), std::ios::out | std::ios::app);
                sigma2out.open(sigma2path.c_str(), std::ios::out | std::ios::app);
	      }
              else {
                errmess = std::string("C++ Error: Incorrect flag for writing to ") + wpath + ". ";
		returnR error(errmess, 99);
                throw error;
              }

    if (!wout || err) {
      errmess = std::string("C++ Error: Could not open ") + wpath + " for writing. ";
      returnR error(errmess, 99);
      throw error; 
    }
    if (!muout || err) {
      errmess = std::string("C++ Error: Could not open ") + mupath + " for writing. ";
      returnR error(errmess, 99);
      throw error; 
    }
    if (!sigma2out || err) {
      errmess = std::string("C++ Error: Could not open ") + sigma2path + " for writing. ";
      returnR error(errmess, 99);
      throw error; 
    }

    std::ostringstream s;
    unsigned int mlen = width;

    // Write to files
    // ==============
    int k_now;
    int lmix = 1 + 3*kmax;              // length of one mixture in mixtureA

    // * Write mixture weights 
    // =======================
    /* Passes up to 5 rows to get things to line up nicely */
    for (int i = 0; i < nR && i < 5; i++) {
      k_now = int(mixtureA[i*lmix]);
      for (int j = 0; j < k_now; j++) {
	s.str("");        
        if (mixtureA[i*lmix + 1 + j] >= FLT_MAX){
	  s << std::setw(width)
	    << std::setiosflags(std::ios::fixed)
	    << "1e50" << "   ";
        }
        else{
          if (mixtureA[i*lmix + 1 + j] < 1 && mixtureA[i*lmix + 1 + j] > -1) 
            s << std::scientific << std::setw(width) << std::setprecision(prec) << mixtureA[i*lmix + 1 + j] << "   ";
          else
            s << std::fixed << std::setw(width) << std::setprecision(prec) << mixtureA[i*lmix + 1 + j] << "   ";
        }
	if (s.str().length() > mlen) mlen = s.str().length();
      }
    }

    /* Write */
//    s.str("");
    for (int i = 0; i < nR; i++) {
      k_now = int(mixtureA[i*lmix]);
      for (int j = 0; j < k_now; j++){
        if (mixtureA[i*lmix + 1 + j] >= FLT_MAX){                                                 
          wout << std::setw(mlen) << "1e50";
          wout << "   ";
	}
        else{
          if (mixtureA[i*lmix + 1 + j] < 1 && mixtureA[i*lmix + 1 + j] > -1){
            wout << std::scientific << std::setw(mlen) << std::setprecision(prec) << mixtureA[i*lmix + 1 + j];
            wout << "   ";
          }
          else{
            wout << std::fixed << std::setw(mlen) << std::setprecision(prec) << mixtureA[i*lmix + 1 + j];
            wout << "   ";
          }
        }
      }
      wout << endl;
    }
//    wout << s.str();
    wout.close();

    // * Write mixture means
    // =======================
    /* Passes up to 5 rows to get things to line up nicely */
    for (int i = 0; i < nR && i < 5; i++) {
      k_now = int(mixtureA[i*lmix]);
      for (int j = 0; j < k_now; j++) {
	s.str("");        
        if (mixtureA[i*lmix + 1 + kmax + j] >= FLT_MAX){
	  s << std::setw(width)
	    << std::setiosflags(std::ios::fixed)
	    << "1e50" << "   ";
        }
        else{
          if (mixtureA[i*lmix + 1 + kmax + j] < 1 && mixtureA[i*lmix + 1 + kmax + j] > -1) 
            s << std::scientific << std::setw(width) << std::setprecision(prec) << mixtureA[i*lmix + 1 + kmax + j] << "   ";
          else
            s << std::fixed << std::setw(width) << std::setprecision(prec) << mixtureA[i*lmix + 1 + kmax + j] << "   ";
        }
	if (s.str().length() > mlen) mlen = s.str().length();
      }
    }

    /* Write */
//    s.str("");
    for (int i = 0; i < nR; i++) {
      k_now = int(mixtureA[i*lmix]);
      for (int j = 0; j < k_now; j++){
        if (mixtureA[i*lmix + 1 + kmax + j] >= FLT_MAX){                                                 
          muout << std::setw(mlen) << "1e50";
          muout << "   ";
	}
        else{
          if (mixtureA[i*lmix + 1 + kmax + j] < 1 && mixtureA[i*lmix + 1 + kmax + j] > -1){
            muout << std::scientific << std::setw(mlen) << std::setprecision(prec) << mixtureA[i*lmix + 1 + kmax + j];
            muout << "   ";
          }
          else{
            muout << std::fixed << std::setw(mlen) << std::setprecision(prec) << mixtureA[i*lmix + 1 + kmax + j];
            muout << "   ";
          }
        }
      }
      muout << endl;
    }
//    muout << s.str();
    muout.close();

    // * Write mixture variances
    // ==========================
    /* Passes up to 5 rows to get things to line up nicely */
    for (int i = 0; i < nR && i < 5; i++) {
      k_now = int(mixtureA[i*lmix]);
      for (int j = 0; j < k_now; j++) {
	s.str("");        
        if (mixtureA[i*lmix + 1 + 2*kmax + j] >= FLT_MAX){
	  s << std::setw(width)
	    << std::setiosflags(std::ios::fixed)
	    << "1e50" << "   ";
        }
        else{
          if (mixtureA[i*lmix + 1 + 2*kmax + j] < 1 && mixtureA[i*lmix + 1 + 2*kmax + j] > -1) 
            s << std::scientific << std::setw(width) << std::setprecision(prec) << mixtureA[i*lmix + 1 + 2*kmax + j] << "   ";
          else
            s << std::fixed << std::setw(width) << std::setprecision(prec) << mixtureA[i*lmix + 1 + 2*kmax + j] << "   ";
        }
	if (s.str().length() > mlen) mlen = s.str().length();
      }
    }

    /* Write */
//    s.str("");
    for (int i = 0; i < nR; i++) {
      k_now = int(mixtureA[i*lmix]);
      for (int j = 0; j < k_now; j++){
        if (mixtureA[i*lmix + 1 + 2*kmax + j] >= FLT_MAX){                                                 
          sigma2out << std::setw(mlen) << "1e50";
          sigma2out << "   ";
	}
        else{
          if (mixtureA[i*lmix + 1 + 2*kmax + j] < 1 && mixtureA[i*lmix + 1 + 2*kmax + j] > -1){
            sigma2out << std::scientific << std::setw(mlen) << std::setprecision(prec) << mixtureA[i*lmix + 1 + 2*kmax + j];
            sigma2out << "   ";
          }
          else{
            sigma2out << std::fixed << std::setw(mlen) << std::setprecision(prec) << mixtureA[i*lmix + 1 + 2*kmax + j];
            sigma2out << "   ";
          }
        }
      }
      sigma2out << endl;
    }
//    sigma2out << s.str();
    sigma2out.close();

  }  // end of try
  catch(returnR){
    throw;
  }  
}    // end of function writeMixtureToFiles


// =====================================================================================
// ****** giveString: Convert integer to std::string ended by ".sim"
// =====================================================================================
std::string
giveString(const int num){

  int sign = (num < 0 ? (-1) : 1);
  if (num == 0) return std::string("0");

  int anum = sign*num;
  int decimal = int(floor(log10(double(anum))) + 1);

  std::string result(".sim");
  int digit;
  int prevdivisor = 1;
  int divisor = 10;
  for (int i = 0; i < decimal; i++){
    digit = anum % divisor;
    digit = digit / prevdivisor;
    prevdivisor = divisor;
    divisor *= 10;
    switch (digit){
    case 0: result = std::string("0") + result; break;
    case 1: result = std::string("1") + result; break;
    case 2: result = std::string("2") + result; break;
    case 3: result = std::string("3") + result; break;
    case 4: result = std::string("4") + result; break;
    case 5: result = std::string("5") + result; break;
    case 6: result = std::string("6") + result; break;
    case 7: result = std::string("7") + result; break;
    case 8: result = std::string("8") + result; break;
    case 9: result = std::string("9") + result; break;
    }     
  }
  if (sign < 0) result = std::string("-") + result;

  return(result);
}

