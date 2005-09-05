// Input/output functions for the Bayesian AFT model, version 1

// 24/11/2003: 'storeInArrays'
//             'writeToFiles'
//             'writeToFiles2'
//             'writeToFiles3'
//             'giveString'
// 06/08/2004: 'readMixtureFromFiles' 
// 08/08/2004: 'writeMixtureToFiles' 

#include "in_output.h"

using namespace std;


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
  int i, j, obs;

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

  if (*storerP) writeToFile(rA, nstored, *nP, dir, "/r.sim", write_flag);      // here: R index was computed already by 'storeInArrays' function         
  else          writeAddToFile(rA, 1, *nP, 1, dir, "/r.sim", 'o');             // here: give R index

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
  int obs;

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
// array........ array where to read it        (length = nread * (1+3*kmax))
// nread ....... on OUTPUT: number of read mixtures
// nR .......... number of rows in files THAT ARE AVAILABLE after skipping possible header
//               (this will usually be MCMC sample size)
// skip ........ how many rows at the beginning of the file are to be skipped (not counting header)
// by .......... only every by-th mixture will be read
// dir ......... directory where all files are stored
// kname ....... name of the file where 'k' is stored in the first column
// wname ....... name of the file where mixture weights are stored (1 mixture/row)
// muname ...... name of the file where mixture means are stored (1 mixture/row)
// sigma2name .. name of the file where mixture variances are stored (1 mixture/row)
//
void
readMixtureFromFiles(double* array,                  int* nread,
                     const int& nR,                  const int& skip,     const int& by,     const int& kmax,
                     const std::string& dir,  
                     const std::string& kname,    
                     const std::string& wname,    
                     const std::string& muname,  
                     const std::string& sigma2name)
{
  try{
    int size = nR * kmax;
    if (size <= 0) throw returnR("C++ Error: Files of null size are to be read by 'readMixtureFromFiles'", 99);
    if (skip < 0) throw returnR("C++ Error: 'skip' parameter must be >= 0 in 'readMixtureFromFiles'", 1);
    if (by <= 0) throw returnR("C++ Error: 'by' parameter must be > 0 in 'readMixtureFromFiles'", 1);
    if (skip >= nR) throw returnR("C++ Error: too many rows are to be skipped by 'readMixtureFromFiles'", 1);

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

    // Skip rows at the beginning that are to be skipped (also header)
    int kcurrent;
    char ch;
    Rprintf("Reading mixture files. \n");
    for (int i = 0; i < skip + 1; i++){
      kfile.get(ch);        while (ch != '\n') kfile.get(ch);
      wfile.get(ch);        while (ch != '\n') wfile.get(ch);
      mufile.get(ch);       while (ch != '\n') mufile.get(ch);
      sigma2file.get(ch);   while (ch != '\n') sigma2file.get(ch);
    }

    // Read the first mixture
    *nread = 1;
    double* veld = array;
    if (kfile.eof()){
      errmes = std::string("C++ Error: Reached end of file ") + kpath + " before "
               + char(*nread) + std::string(" values were read.");
      throw returnR(errmes, 99);
    }
    kfile >> kcurrent;
    if (kcurrent > kmax) throw returnR("C++ Error: k value higher than kmax was read.", 99);
    *veld = kcurrent;
    kfile.get(ch);                     while (ch != '\n') kfile.get(ch);

    veld++;
    if (wfile.eof()){
      errmes = std::string("C++ Error: Reached end of file ") + wpath + " before "
               + char(*nread) + std::string(" sets of mixture weights were read.");
      throw returnR(errmes, 99);
    }
    for (int j = 0; j < kcurrent; j++) wfile >> *(veld + j);
    wfile.get(ch);                 while (ch != '\n') wfile.get(ch);

    veld += kmax;
    if (mufile.eof()){
      errmes = std::string("C++ Error: Reached end of file ") + mupath + " before "
               + char(*nread) + std::string(" sets of mixture means were read.");
      throw returnR(errmes, 99);
    }
    for (int j = 0; j < kcurrent; j++) mufile >> *(veld + j);
    mufile.get(ch);                while (ch != '\n') mufile.get(ch);

    veld += kmax;
    if (sigma2file.eof()){
      errmes = std::string("C++ Error: Reached end of file ") + sigma2path + " before "
               + char(*nread) + std::string(" sets of mixture variances were read.");
      throw returnR(errmes, 99);
    }
    for (int j = 0; j < kcurrent; j++) sigma2file >> *(veld + j);
    sigma2file.get(ch);            while (ch != '\n') sigma2file.get(ch);

    // Read remaining mixtures
    for (int i = skip + 1 + by; i <= nR; i += by){

      // Skip by-1 mixtures
      for (int ii = 0; ii < by - 1; ii++){
        kfile.get(ch);        while (ch != '\n') kfile.get(ch);
        wfile.get(ch);        while (ch != '\n') wfile.get(ch);
        mufile.get(ch);       while (ch != '\n') mufile.get(ch);
        sigma2file.get(ch);   while (ch != '\n') sigma2file.get(ch);
      }

      // Read the values
      (*nread)++;
      veld += kmax;
      if (kfile.eof()){
        errmes = std::string("C++ Error: Reached end of file ") + kpath + " before "
                 + char(*nread) + std::string(" values were read.");
        throw returnR(errmes, 99);
      }
      kfile >> kcurrent;
      if (kcurrent > kmax) throw returnR("C++ Error: k value higher than kmax was read.", 99);
      *veld = kcurrent;
      kfile.get(ch);                     while (ch != '\n') kfile.get(ch);

      veld++;
      if (wfile.eof()){
        errmes = std::string("C++ Error: Reached end of file ") + wpath + " before "
                 + char(*nread) + std::string(" sets of mixture weights were read.");
        throw returnR(errmes, 99);
      }
      for (int j = 0; j < kcurrent; j++) wfile >> *(veld + j);
      wfile.get(ch);                 while (ch != '\n') wfile.get(ch);

      veld += kmax;
      if (mufile.eof()){
        errmes = std::string("C++ Error: Reached end of file ") + mupath + " before "
                 + char(*nread) + std::string(" sets of mixture means were read.");
        throw returnR(errmes, 99);
      }
      for (int j = 0; j < kcurrent; j++) mufile >> *(veld + j);
      mufile.get(ch);                while (ch != '\n') mufile.get(ch);

      veld += kmax;
      if (sigma2file.eof()){
        errmes = std::string("C++ Error: Reached end of file ") + sigma2path + " before "
                 + char(*nread) + std::string(" sets of mixture variances were read.");
        throw returnR(errmes, 99);
      }
      for (int j = 0; j < kcurrent; j++) sigma2file >> *(veld + j);
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

