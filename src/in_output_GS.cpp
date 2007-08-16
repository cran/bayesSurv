// Some helping functions for a Bayesian estimation of the density
//  using the G-splines
//
// 23/10/2004: 'storeInArrays_bayesHistogram'
//             'writeToFiles_bayesHistogram'
//             'findClosestKnot'
// 01/11/2004: 'checkGsplineFiles'
//             'readGsplineFromFiles'
//             'closeGsplineFiles'
// 03/11/2004: 'openFilesForOutput'
//             'openFiles_bayesHistogram'
// 23/01/2005: 'readGsplineFromFiles2'
//             'openRegresFiles'
//             'readRegresFromFiles'
//             'closeRegresFiles'
// 30/01/2005: 'writeToFiles_random' 
// 02/02/2005: 'writeToFiles_Gspl_intcpt'
//             'open_File_toRead'
//             'adjust_intercept'
// 06/02/2005: 'readGsplineFromFiles3'
// 16/05/2005: 'readMean_and_Scale'
//             'adjust_intercept' slightly changed 
// 26/11/2005: 'openGsplineFiles_forTau'
// 27/11/2005: 'readGsplineFiles_forTau'
//             'readGsplineFiles_forMarginal'
// 12/12/2006: 'writeToFiles_RandomEff32'
// 13/12/2006: 'openD32File'
//
#include "in_output_GS.h"

using namespace std;

// ====================================================================================
// ***** openFiles_bayesHistogram: open files to write simulated values          *****/
// *****                   version for bayesHistogram                            *****/
// ====================================================================================
//
void
openFiles_bayesHistogram(
    std::ofstream& sigmafile,          std::ofstream& lambdafile,
    std::ofstream& mixmomentfile,      std::ofstream& mweightfile,      std::ofstream& mlogweightfile,
    std::ofstream& mmeanfile,          std::ofstream& Yfile,            std::ofstream& rfile,              
    std::ofstream& logposterfile,
    const std::string& sigmapath,      const std::string& lambdapath,
    const std::string& mixmomentpath,  const std::string& mweightpath,  const std::string& mlogweightpath,
    const std::string& mmeanpath,      const std::string& Ypath,        const std::string& rpath,          
    const std::string& logposterpath,
    const int& n_censored,             const char& write_flag)
{
  openFile(sigmafile, sigmapath, write_flag);
  openFile(lambdafile, lambdapath, write_flag);
  openFile(mixmomentfile, mixmomentpath, write_flag);
  openFile(mweightfile, mweightpath, write_flag);
  openFile(mlogweightfile, mlogweightpath, write_flag);
  openFile(mmeanfile, mmeanpath, write_flag);
  if (n_censored) openFile(Yfile, Ypath, write_flag);
  openFile(rfile, rpath, write_flag);
  openFile(logposterfile, logposterpath, write_flag);

  return;
}

// ====================================================================================
// ***** closeFiles_bayesHistogram: close files after writing simulated values   *****/
// *****                   version for bayesHistogram                            *****/
// ====================================================================================
//
void
closeFiles_bayesHistogram(
    std::ofstream& sigmafile,          std::ofstream& lambdafile,
    std::ofstream& mixmomentfile,      std::ofstream& mweightfile,      std::ofstream& mlogweightfile,
    std::ofstream& mmeanfile,          std::ofstream& Yfile,            std::ofstream& rfile,              
    std::ofstream& logposterfile,
    const int& n_censored)
{
  sigmafile.close();            
  lambdafile.close();           
  mixmomentfile.close();        
  mweightfile.close();          
  mlogweightfile.close();       
  mmeanfile.close();            
  if (n_censored) Yfile.close();
  rfile.close();                
  logposterfile.close();        

  return;
}


// ====================================================================================
// ***** writeToFiles_bayesHistogram: write simulated values to files            *****/
// *****                   version for bayesHistogram                            *****/
// ====================================================================================
//
// writeAll ..... 0/1, write Y and r even if they are not to be stored
//                * this serves for writing Y and r after some given number of iterations to be able
//                  to start McMC again in case of some problems
//                * possible 1 will always be reset to zero on return
//
// check_k_effect ... if positive (DEFAULT) then only components with high enough weights are written to files
//                    otherwise all components are written
//
void
writeToFiles_bayesHistogram(
    const Gspline* gg,     
    const int* rM,                     const double* YsM,               double* log_poster,
    const int& l_momentsA,             const int& l_lambdaA,            const int& l_log_poster,
    const int* nP,                     const int* storeaP,              const int* storeyP,              
    const int* storerP,                const int& n_censored,
    const int* writeAll,               int* workI,                      double* workD,
    std::ofstream& sigmafile,          std::ofstream& lambdafile,
    std::ofstream& mixmomentfile,      std::ofstream& mweightfile,      std::ofstream& mlogweightfile,
    std::ofstream& mmeanfile,          std::ofstream& Yfile,            std::ofstream& rfile,              
    std::ofstream& logposterfile,
    const double& null_weight,         const int& prec,                 const int& width,
    const int& check_k_effect)
{
  int i, K0;
  double sumexpa;

  /*** Middle knots, basis standard deviations and distance between two knots ***/  
  writeFiveToFile_1(gg->gammaP(), gg->sigmaP(), gg->deltaP(), gg->intcptP(), gg->scaleP(), 
                    gg->dim(), gg->dim(), gg->dim(), gg->dim(), gg->dim(),
                    sigmafile, prec, width);

  /*** Lambda ***/
  writeToFile_1(gg->lambdaP(), l_lambdaA, lambdafile, prec, width);

  /*** Mixture weights, indeces of means and effective number of components ***/
  static int k_effect_write;
  static double gewicht;
  static double *pworkD, *expaP;
  static int *pworkI;
  if (check_k_effect){
    k_effect_write = 0;
    pworkD = workD;
    pworkI = workI;
    for (i = 0; i < gg->k_effect(); i++){
      gewicht = gg->w(gg->ind_w_effect(i));
      if (gewicht >= null_weight){
        *pworkD = gewicht;
  /*      for (j = 0; j < gg->dim(); j++){                                                                   */
  /* 	    muA[nstored*l_muA + k_effect_write*gg->dim() + j] = gg->mu_component(j, gg->ind_w_effect(i));    */
  /*      }                                                                                                  */
        switch (gg->dim()){
        case 1: *pworkI = gg->ind_w_effect(i) - gg->K(0);
	        break;
        case 2: *pworkI = gg->ind_w_effect(i) % gg->length(0) - gg->K(0);
         	pworkI++;
	        *pworkI = gg->ind_w_effect(i) / gg->length(0) - gg->K(1);
                break;
        default: throw returnR("C++ Error: Unimplemented part (dim > 2) of the function writeToFiles_bayesHistogram", 1);
        }
        pworkD++;
        pworkI++;
        k_effect_write++;
      } 
    }
  }
  else{
    k_effect_write = gg->total_length();
    K0 = gg->K(0);
    sumexpa = gg->sumexpa();

    if (gg->dim() > 1) throw returnR("C++ Error: check_k_effect must be > 0 if dim > 1 in writeToFiles_bayesHistogram", 1);
    pworkI = workI;
    pworkD = workD;
    expaP  = gg->expaP();
    for (i = 0; i < k_effect_write; i++){
      *pworkI = i - K0; 
      *pworkD = *expaP/sumexpa;
      pworkI++;
      pworkD++;
      expaP++;
    }
  }
  writeToFile_1(workD, k_effect_write, mweightfile, prec, width);
  writeToFile_1(workI, gg->dim()*k_effect_write, mmeanfile, prec, width);


  /*** Mixture moments ***/
  gg->moments(workD, workD + gg->dim());
  writeTwoToFile_1(&k_effect_write, workD, 0, l_momentsA, mixmomentfile, prec, width);
  
  /*** Mixture a coefficients ***/
  if (*storeaP || *writeAll){
    writeToFile_1(gg->aP(), gg->total_length(), mlogweightfile, prec, width);
  }

  /*** Sampled (augmented) observations ***/
  if ((*storeyP || *writeAll) && n_censored){
    writeToFile_1(YsM, (*nP) * gg->dim() , Yfile, prec, width);
  }

  /*** Labels of components ***/
  if (*storerP || *writeAll){      
    writeAddToFile_1(rM, *nP, 1, rfile, prec, width);     /** give R indeces **/
  } 

  /*** Log-posterior ***/
  for (i = 0; i < (gg->equal_lambda() ? 1 : gg->dim()); i++) log_poster[1 + i] = gg->penalty(i);
  writeToFile_1(log_poster, l_log_poster, logposterfile, prec, width);

  return;
}


// ====================================================================================
// ***** storeInArrays_bayesHistogram: store simulated values in working arrays  *****/
// *****                   version for bayesHistogram                            *****/
//   currently not used
// ====================================================================================
//
// muA ............. indeces of means corresponding to non-zero weights
//                   * indeces on scale -K,...,0,...,K will be stored here
//                   * if (dim == 2) then pairs of indeces are stored close to each other
//
// null_weight ..... only mixture components with w >= null_weight will be recorded
//
void
storeInArrays_bayesHistogram(
    int* iterA,             int* k_effectA,    double* momentsA,  
    double* weightsA,       int* muA,          double* gamma_sigma_deltaA,
    double* lambdaA,        int* rA,           double* YsA,
    const int& l_momentsA,
    const int& l_weightsA,  const int& l_muA,  const int& l_gamma_sigma_deltaA,
    const int& l_lambdaA,   const int* nP,
    const int& iterindex,   const Gspline* gg, 
    const int* rM,          const double* YsM,
    const int& nstored,     const double& null_weight,
    const int* storeyP,     const int* storerP)
{
  int i, j;
  double gewicht;

  iterA[nstored] = iterindex;

  /*** Mixture weights, indeces of means and effective number of components ***/
  int k_effect_write = 0;
  for (i = 0; i < gg->k_effect(); i++){
    gewicht = gg->w(gg->ind_w_effect(i));
    if (gewicht >= null_weight){
      weightsA[nstored*l_weightsA + k_effect_write] = gewicht;
/*      for (j = 0; j < gg->dim(); j++){                                                                   */
/*	  muA[nstored*l_muA + k_effect_write*gg->dim() + j] = gg->mu_component(j, gg->ind_w_effect(i));    */
/*      }                                                                                                  */
      switch (gg->dim()){
      case 1: muA[nstored*l_muA + k_effect_write*gg->dim()] = gg->ind_w_effect(i) - gg->K(0);
	      break;
      case 2: muA[nstored*l_muA + k_effect_write*gg->dim()]     = gg->ind_w_effect(i) % gg->length(0) - gg->K(0);
              muA[nstored*l_muA + k_effect_write*gg->dim() + 1] = gg->ind_w_effect(i) / gg->length(0) - gg->K(1);
              break;
      default: throw returnR("C++ Error: Unimplemented part (dim > 2) of the function storeInArrays_bayesHistogram", 1);
      }
      k_effect_write++;
    }
  }
  k_effectA[nstored] = k_effect_write;

  /*** Middle knots, basis standard deviations and distance between two knots ***/
  for (j = 0; j < gg->dim(); j++){
    gamma_sigma_deltaA[nstored*l_gamma_sigma_deltaA + j]               = gg->gamma(j);
    gamma_sigma_deltaA[nstored*l_gamma_sigma_deltaA + gg->dim() + j]   = gg->sigma(j);
    gamma_sigma_deltaA[nstored*l_gamma_sigma_deltaA + 2*gg->dim() + j] = gg->delta(j);
  }

  /*** Lambda precision parameters ***/
  for (j = 0; j < l_lambdaA; j++){
    lambdaA[nstored*l_lambdaA + j] = gg->lambda(j);
  }

  /*** Mixture moments ***/
  int beg = nstored*l_momentsA;
  gg->moments(momentsA + beg, momentsA + beg + gg->dim());

  /*** Sampled (augmented) observations ***/
  if (*storeyP){
    for (i = 0; i < *nP; i++) YsA[nstored*(*nP) + i] = YsM[i];
  }

  /*** Labels of components ***/
  if (*storerP){      /** give R indeces **/
    for (i = 0; i < *nP; i++) rA[nstored*(*nP) + i] = rM[i] + 1;
  } 

  return;
}   /** end of function storeInArrays_bayesHistogram **/


// ====================================================================================
// ***** writeToFiles2_bayesHistogram: write simulated values to files           *****/
// *****                   version for bayesHistogram                            *****/
// *****                     (currently not used)                                *****/
// ====================================================================================
void
writeToFiles2_bayesHistogram(
    const int* iterA,        const int* k_effectA,    const double* momentsA,
    const double* weightsA,  const int* muA,          const double* gamma_sigma_deltaA,
    const double* lambdaA,   const int* rA,           const double* YsA,
    const int& nstored,      const std::string& dir,  const char& write_flag,
    const int& l_momentsA,
    const int& l_weightsA,   const int& l_muA,        const int& l_gamma_sigma_deltaA,
    const int& l_lambdaA,    const int* dimP,         const int* nP,
    const int* storeyP,      const int* storerP,      const int* n_censored)
{
  writeToFile(iterA, nstored, 1, dir, "/iteration.sim", write_flag);
  writeToFile(gamma_sigma_deltaA, nstored, l_gamma_sigma_deltaA, dir, "/gamma_sigma_delta.sim", write_flag);
  writeToFile(lambdaA, nstored, l_lambdaA, dir, "/lambda.sim", write_flag);
  writeTwoToFile(k_effectA, nstored, 1, 0, momentsA, nstored, l_momentsA, dir, "/mixmoment.sim", write_flag);
  writeRaggedToFile(weightsA, nstored, l_weightsA, k_effectA, 1, dir, "/mweight.sim", write_flag);
  writeRaggedToFile(muA, nstored, l_muA, k_effectA, *dimP, dir, "/mmean.sim", write_flag); 

  if (*n_censored){
    if (*storeyP) writeToFile(YsA, nstored, (*dimP)*(*nP), dir, "/Y.sim", write_flag);             
    else          writeToFile(YsA, 1, (*dimP)*(*nP), dir, "/Y.sim", 'o');             
  }

  if (*storerP) writeToFile(rA, nstored, *nP, dir, "/r.sim", write_flag);      // here: R index was computed already by 'storeInArrays' function         
  else          writeAddToFile(rA, 1, *nP, 1, dir, "/r.sim", 'o');             // here: give R index

  return;
}   /** end of function writeToFiles_bayesHistogram **/


// ==============================================================================================================
// ***** openGsplineFiles: Open files where sampled G-spline are stored for reading
//                          and skip first 'skip' rows that the user wishes to skip
// ==============================================================================================================
void
openGsplineFiles(std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,      std::ifstream& sigmafile,
                 const std::string& kpath,  const std::string& wpath,  const std::string& mupath,  const std::string& sigmapath,
                 const int& skip)
{
  open_File_toRead(kfile, kpath, skip);
  open_File_toRead(wfile, wpath, skip);
  open_File_toRead(mufile, mupath, skip);
  open_File_toRead(sigmafile, sigmapath, skip);

  return;
}


// ===========================================================================================================================
// ***** openGsplineFiles_forTau: Open files where sampled G-spline are stored for reading (to be used by sampledKendallTau)
//                                and skip first 'skip' rows that the user wishes to skip
// ===========================================================================================================================
void
openGsplineFiles_forTau(std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,
                        const std::string& kpath,  const std::string& wpath,  const std::string& mupath,
                        const int& skip)
{
  open_File_toRead(kfile, kpath, skip);
  open_File_toRead(wfile, wpath, skip);
  open_File_toRead(mufile, mupath, skip);

  return;
}


// ==================================================================================
// ***** open_File_toRead: Open file for reading and skip first 'skip' rows
//   * written first to be used by 'bayesGspline' to open a file
//     with adjustment intercept
//   * then incorporated into 'openGsplineFiles'
// ==================================================================================
void
open_File_toRead(std::ifstream& file,  const std::string& path,  const int& skip)
{
  std::string errmes;

  file.open(path.c_str(), std::ios::in);
  if (!file){
    errmes = std::string("Error: Could not open ") + path;
    throw returnR(errmes, 99);
  }

  char ch;
  for (int i = 0; i < skip; i++){
    if (file.eof()){
      int ihelp = i + 1;
      errmes = std::string("Error: Reached end of file ") + path + " before "
               + char(ihelp) + std::string(" rows were skipped.");
      throw returnR(errmes, 99);
    }
    file.get(ch);        
    while (ch != '\n') file.get(ch);
  }

  return;
}


// ==================================================================================
// ***** readGsplineFromFiles: Function to read one sampled G-spline from four files 
//   version used by 'bayesGspline'
// ==================================================================================
//
// k_effect .................... effective number of mixture components read
// w[total_length] ............. read weights, only first k_effect components of this array are filled
// ind_mu[dim][total_length] ... indeces of means (on the scale -K,...,K) corresponding to w
// mu[dim][total_length] ....... means corresponding to w
// gamma[dim] .................. middle knots in each dimension
// sigma[dim] .................. basis standard deviations in each dimension
// delta[dim] .................. delta parameter in each dimension
// intcpt[dim] ................. intercept in each dimension
// scale[dim] .................. scale parameter in each dimension
// skip ........................ number of rows that should be skipped before reading real data
// row ......................... how many rows of the data will be read at the end of the function call
// dim ......................... dimension of the G-spline
// total_length ................ maximal number of components in the G-spline
// kfile ....................... usually stream associated with 'mixmoment*.sim'
// wfile ....................... usually stream associated with 'mweight*.sim'
// mufile ...................... usually stream associated with 'mmean*.sim'
// sigmafile ................... usually stream associated with 'gspline*.sim'
//
void
readGsplineFromFiles(int* k_effect,             double* w,                 int** ind_mu,               double** mu,
                     double* gamma,             double* sigma,             double* delta,
                     double* intcpt,            double* scale,
                     const int& skip,           const int& row,            const int& dim,             const int& total_length,
                     std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,      std::ifstream& sigmafile,
                     const std::string& kpath,  const std::string& wpath,  const std::string& mupath,  const std::string& sigmapath)
{
  try{
    static int j, dd, k_current, ihelp;
    static char ch;
    static std::string errmes;

    /**  Skip rows that are to be skipped  **/
    for (j = 0; j < skip; j++){
      kfile.get(ch);        while (ch != '\n') kfile.get(ch);
      wfile.get(ch);        while (ch != '\n') wfile.get(ch);
      mufile.get(ch);       while (ch != '\n') mufile.get(ch);
      sigmafile.get(ch);    while (ch != '\n') sigmafile.get(ch);
    }

    /**  Read effective k **/
    if (kfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + kpath + " before "
               + char(ihelp) + std::string(" values were read.");
      throw returnR(errmes, 99);
    }
    kfile >> k_current;
    if (k_current > total_length) throw returnR("C++ Error: k value higher than indicated total_length of the G-spline was read.", 99);
    *k_effect = k_current;
    kfile.get(ch);                     while (ch != '\n') kfile.get(ch);

    /** Read G-spline weights **/
    if (wfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + wpath + " before "
               + char(ihelp) + std::string(" sets of G-spline weights were read.");
      throw returnR(errmes, 99);
    }
    for (j = 0; j < k_current; j++) wfile >> w[j];
    wfile.get(ch);                 while (ch != '\n') wfile.get(ch);

    /** Read G-spline intercept, standard deviations and distance between two knots **/
    if (sigmafile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + sigmapath + " before "
               + char(ihelp) + std::string(" sets of G-spline intercepts/std. deviations were read.");
      throw returnR(errmes, 99);
    }
    for (dd = 0; dd < dim; dd++) sigmafile >> gamma[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> sigma[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> delta[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> intcpt[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> scale[dd];
    sigmafile.get(ch);            while (ch != '\n') sigmafile.get(ch);

    /** Read indeces of G-spline means and calculate means **/
    if (mufile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + mupath + " before "
               + char(ihelp) + std::string(" sets of G-spline means were read.");
      throw returnR(errmes, 99);
    }
    for (j = 0; j < k_current; j++){
      for (dd = 0; dd < dim ; dd++){
        mufile >> ind_mu[dd][j];
        mu[dd][j] = gamma[dd] + ind_mu[dd][j]*delta[dd];
      }
    }
    mufile.get(ch);                while (ch != '\n') mufile.get(ch);

    return;    
  }  // end of try
  catch(returnR){
    throw;
  }  
}


// ===============================================================================================
// ***** readGsplineFromFiles_forMarginal: Function to read one sampled G-spline from four files 
//   version used by 'marginal_bayesGspline'
//
// --> transform also ind_mu into the scale 0,..., 2*K
// --> compute marginal weights
//  (this is different from 'readGsplineFromFiles' function!!!)
//
// ===============================================================================================
//
// w_temp[total_length] ........ read weights, only first k_effect components of this array are filled
// w[dim][2*KK[j]+1] ........... weights in each margin
// mu[dim][2*KK[j]+1] .......... means in each margin
// gamma[dim] .................. middle knots in each dimension
// sigma[dim] .................. basis standard deviations in each dimension
// delta[dim] .................. delta parameter in each dimension
// intcpt[dim] ................. intercept in each dimension
// scale[dim] .................. scale parameter in each dimension
// KK[dim] ..................... numbers of knots at each side of the reference knot
// skip ........................ number of rows that should be skipped before reading real data
// row ......................... how many rows of the data will be read at the end of the function call
// dim ......................... dimension of the G-spline
// total_length ................ maximal number of components in the G-spline
// kfile ....................... usually stream associated with 'mixmoment*.sim'
// wfile ....................... usually stream associated with 'mweight*.sim'
// mufile ...................... usually stream associated with 'mmean*.sim'
// sigmafile ................... usually stream associated with 'gspline*.sim'
//
void
readGsplineFromFiles_forMarginal
   (double* w_temp,            double** w,                double** mu,
    double* gamma,             double* sigma,             double* delta,
    double* intcpt,            double* scale,
    const int* KK,             
    const int& skip,           const int& row,            const int& total_length,
    std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,      std::ifstream& sigmafile,
    const std::string& kpath,  const std::string& wpath,  const std::string& mupath,  const std::string& sigmapath)
{
  try{
    const int dim = 2;
    static int ind0, ind1;
    static int j, dd, k_current, ihelp;
    static char ch;
    static std::string errmes;

    /** Reset marginal weights  **/
    for (dd = 0; dd < dim; dd++){
      for (j = 0; j < 2*KK[dd] + 1; j++){
        w[dd][j] = 0.0;
      }
    }

    /**  Skip rows that are to be skipped  **/
    for (j = 0; j < skip; j++){
      kfile.get(ch);        while (ch != '\n') kfile.get(ch);
      wfile.get(ch);        while (ch != '\n') wfile.get(ch);
      mufile.get(ch);       while (ch != '\n') mufile.get(ch);
      sigmafile.get(ch);    while (ch != '\n') sigmafile.get(ch);
    }

    /**  Read effective k **/
    if (kfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + kpath + " before "
               + char(ihelp) + std::string(" values were read.");
      throw returnR(errmes, 99);
    }
    kfile >> k_current;
    if (k_current > total_length) throw returnR("C++ Error: k value higher than indicated total_length of the G-spline was read.", 99);
    kfile.get(ch);                     while (ch != '\n') kfile.get(ch);

    /** Read G-spline weights **/
    if (wfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + wpath + " before "
               + char(ihelp) + std::string(" sets of G-spline weights were read.");
      throw returnR(errmes, 99);
    }
    for (j = 0; j < k_current; j++) wfile >> w_temp[j];
    wfile.get(ch);                 while (ch != '\n') wfile.get(ch);

    /** Read G-spline intercept, standard deviations and distance between two knots **/
    if (sigmafile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + sigmapath + " before "
               + char(ihelp) + std::string(" sets of G-spline intercepts/std. deviations were read.");
      throw returnR(errmes, 99);
    }
    for (dd = 0; dd < dim; dd++) sigmafile >> gamma[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> sigma[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> delta[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> intcpt[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> scale[dd];
    sigmafile.get(ch);            while (ch != '\n') sigmafile.get(ch);

    /** Calculate means (knots) **/
    for (dd = 0; dd < dim; dd++){
      mu[dd][0] = gamma[dd] - KK[dd]*delta[dd];
      for (j = 1; j < 2*KK[dd]+1; j++) mu[dd][j] = mu[dd][j-1] + delta[dd];
    }

    /** Read indeces of G-spline means and calculate marginal weights **/
    if (mufile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + mupath + " before "
               + char(ihelp) + std::string(" sets of G-spline means were read.");
      throw returnR(errmes, 99);
    }
    for (j = 0; j < k_current; j++){
      mufile >> ind0;
      mufile >> ind1;
      ind0 += KK[0];
      ind1 += KK[1];
      w[0][ind0] += w_temp[j];
      w[1][ind1] += w_temp[j];
    }
    mufile.get(ch);                while (ch != '\n') mufile.get(ch);

    return;    
  }  // end of try
  catch(returnR){
    throw;
  }  
}


// ==================================================================================
// ***** readGsplineFromFiles_forTau: Function to read one sampled G-spline from four files 
//   version used by 'sampledKendallTau'
//
// --> transform also ind_mu into the scale 0,..., 2*K
//     (this is different from 'readGsplineFromFiles' function!!!)
//
// ==================================================================================
//
// k_effect .................... effective number of mixture components read
// w[total_length] ............. read weights, only first k_effect components of this array are filled
// ind_mu[dim][total_length] ... indeces of means (on the scale -K,...,K) corresponding to w
// skip ........................ number of rows that should be skipped before reading real data
// row ......................... how many rows of the data will be read at the end of the function call
// dim ......................... dimension of the G-spline
// KK[dim] ..................... numbers of knots at each side of the reference knot
// total_length ................ maximal number of components in the G-spline
// kfile ....................... usually stream associated with 'mixmoment*.sim'
// wfile ....................... usually stream associated with 'mweight*.sim'
// mufile ...................... usually stream associated with 'mmean*.sim'
//
void
readGsplineFromFiles_forTau(int* k_effect,             double* w,                 int** ind_mu,
                            const int& skip,           const int& row,            const int& dim,             
                            const int* KK,             const int& total_length,
                            std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile, 
                            const std::string& kpath,  const std::string& wpath,  const std::string& mupath)
{
  try{
    static int j, dd, k_current, ihelp;
    static char ch;
    static std::string errmes;

    /**  Skip rows that are to be skipped  **/
    for (j = 0; j < skip; j++){
      kfile.get(ch);        while (ch != '\n') kfile.get(ch);
      wfile.get(ch);        while (ch != '\n') wfile.get(ch);
      mufile.get(ch);       while (ch != '\n') mufile.get(ch);
    }

    /**  Read effective k **/
    if (kfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + kpath + " before "
               + char(ihelp) + std::string(" values were read.");
      throw returnR(errmes, 99);
    }
    kfile >> k_current;
    if (k_current > total_length) throw returnR("C++ Error: k value higher than indicated total_length of the G-spline was read.", 99);
    *k_effect = k_current;
    kfile.get(ch);                     while (ch != '\n') kfile.get(ch);

    /** Read G-spline weights **/
    if (wfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + wpath + " before "
               + char(ihelp) + std::string(" sets of G-spline weights were read.");
      throw returnR(errmes, 99);
    }
    for (j = 0; j < k_current; j++) wfile >> w[j];
    wfile.get(ch);                 while (ch != '\n') wfile.get(ch);

    /** Read indeces of G-spline means  **/
    if (mufile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + mupath + " before "
               + char(ihelp) + std::string(" sets of G-spline means were read.");
      throw returnR(errmes, 99);
    }
    for (j = 0; j < k_current; j++){
      for (dd = 0; dd < dim ; dd++){
        mufile >> ind_mu[dd][j];
        ind_mu[dd][j] += KK[dd];
      }
    }
    mufile.get(ch);                while (ch != '\n') mufile.get(ch);

    return;    
  }  // end of try
  catch(returnR){
    throw;
  }  
}


// ============================================================================================
// ***** readMean_and_Scale: read mixture mean and standard deviation from the file mixmoment
// ============================================================================================
void
readMean_and_Scale(double* E_gx,                  double* sd_gx,
                   const int& skip,               const int& row,  
                   const int& dim,
                   std::ifstream& mixmomentfile,  const std::string& mixmomentpath)
{
  try{
    static int j, ihelp;
    static char ch;
    static std::string errmes;

    if (dim > 1) throw returnR("C++ Error: Function readMean_and_Scale not implemented for dim > 1.", 99);

    /**  Skip rows that are to be skipped  **/
    for (j = 0; j < skip; j++){
      mixmomentfile.get(ch);        
      while (ch != '\n') mixmomentfile.get(ch);
    }

    /** Read G-spline average and variance **/
    if (mixmomentfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + mixmomentpath + " before "
               + char(ihelp) + std::string(" sets of G-spline intercepts/std. deviations were read.");
      throw returnR(errmes, 99);
    }
    mixmomentfile >> j;    /** read number of components **/
    mixmomentfile >> *E_gx;
    mixmomentfile >> *sd_gx;
    if (*sd_gx <= 0) throw returnR("Error: non-positive variance read.", 99);
    *sd_gx = sqrt(*sd_gx);

    mixmomentfile.get(ch);            
    while (ch != '\n') mixmomentfile.get(ch);
  }
  catch(returnR){
    throw;
  }  
}


// ============================================================================================
// ***** adjust_intercept: Function to read adjustment constant which is added or subtracted
//                         from the G-spline intercept
//    * it performs also addition/subtraction
// ============================================================================================
//
// E_gx .......... mean of the original G-spline, it is also adjusted
// skip .......... see 'readGsplineFromFiles' above
// row ........... see 'readGsplineFromFiles' above
// version ....... 30 or 31 (see bayesGspline.cpp for explanation)
// file .......... file where to take adjustment constant
// path .......... path to the file with the adjustment constant
//                 * this will usually ends with "mixmoment_b.sim" or "mixmoment_b2.sim"
//
void
adjust_intercept(double* intcpt,        const int* version,
                 double* E_gx,
                 const int& skip,       const int& row,
                 std::ifstream& file,   const std::string& path)
{
  static int j, ihelp;
  static char ch;
  static double intcpt_adj;
  static std::string errmes;

  /**  Skip rows that are to be skipped  **/
  for (j = 0; j < skip; j++){
    file.get(ch);
    while (ch != '\n') file.get(ch);
  }

  /**  Read adjustment constant **/
  if (file.eof()){
    ihelp = row + 1;
    errmes = std::string("C++ Error: Reached end of file ") + path + " before "
             + char(ihelp) + std::string(" values were read.");
    throw returnR(errmes, 99);
  }
  
  /** Skip the first column with k **/
  file >> intcpt_adj;
  
  /** Read the adjustment constant  and skip the rest of the row **/
  file >> intcpt_adj;
  file.get(ch);
  while (ch != '\n') file.get(ch);

  /** Adjust G-spline intercept **/
  switch (*version){
  case 30:        /** We are computing density of the error term **/
    *intcpt += intcpt_adj;
    *E_gx += intcpt_adj;
    break;
  case 31:        /** We are computing density of the random intercept -> center it  **/
    *intcpt -= intcpt_adj;
    *E_gx -= intcpt_adj;
    break;
  default:
    throw returnR("Error: Strange version appeared in 'adjust_intercept' function", 1);
  }

  return;  
}


// ==================================================================================
// ***** readGsplineFromFiles2: Function to read one sampled G-spline from four files 
//   version used by 'predictive_GS' for the error term
// ==================================================================================
//
// k_effect ...................... effective number of mixture components read
// w_marg[dim][2*K[dd]+1] ....... weights for marginal distributions in each dimension
// mu_sig_marg[dim][2*K[dd]+1] .. means (knots) divided by the basis std. deviation  for marginal distribution in each dimension
// gamma[dim] ................... middle knots in each dimension
// sigma[dim] ................... basis standard deviations in each dimension
// delta[dim] ................... delta parameter in each dimension
// intcpt[dim] .................. intercept in each dimension
// scale[dim] ................... scale parameter in each dimension
// delta_sig[dim] ............... delta/sigma in each dimension
// skip ......................... number of rows that should be skipped before reading real data
// row .......................... how many rows of the data will be read at the end of the function call
// dim .......................... dimension of the G-spline
// total_length ................. maximal number of components in the G-spline
// GsplK[dim] ................... "K" parameters of the G-spline
// kfile ........................ usually stream associated with 'mixmoment*.sim'
// wfile ........................ usually stream associated with 'mweight*.sim'
// mufile ....................... usually stream associated with 'mmean*.sim'
// sigmafile .................... usually stream associated with 'gspline*.sim'
//
void
readGsplineFromFiles2(int* k_effect,             double** w_marg,           double** mu_sig_marg,
                      double* gamma,             double* sigma,             double* delta,
                      double* intcpt,            double* scale,             double* delta_sig,
                      const int& skip,           const int& row,            const int& dim,             
                      const int& total_length,   const int* GsplK,
                      std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,      std::ifstream& sigmafile,
                      const std::string& kpath,  const std::string& wpath,  const std::string& mupath,  const std::string& sigmapath)
{
  try{
    static int j, dd, k_current, ihelp;
    static double tmp;
    static char ch;
    static std::string errmes;

    /**  Skip rows that are to be skipped  **/
    for (j = 0; j < skip; j++){
      kfile.get(ch);        while (ch != '\n') kfile.get(ch);
      wfile.get(ch);        while (ch != '\n') wfile.get(ch);
      mufile.get(ch);       while (ch != '\n') mufile.get(ch);
      sigmafile.get(ch);    while (ch != '\n') sigmafile.get(ch);
    }

    /**  Read effective k **/
    if (kfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + kpath + " before "
               + char(ihelp) + std::string(" values were read.");
      throw returnR(errmes, 99);
    }
    kfile >> k_current;
    if (k_current > total_length) throw returnR("C++ Error: k value higher than indicated total_length of the G-spline was read.", 99);
    *k_effect = k_current;
    kfile.get(ch);                     while (ch != '\n') kfile.get(ch);

    /** Read G-spline intercept, standard deviations and distance between two knots **/
    if (sigmafile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + sigmapath + " before "
               + char(ihelp) + std::string(" sets of G-spline intercepts/std. deviations were read.");
      throw returnR(errmes, 99);
    }
    for (dd = 0; dd < dim; dd++) sigmafile >> gamma[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> sigma[dd];
    for (dd = 0; dd < dim; dd++){
      sigmafile >> delta[dd];
      delta_sig[dd] = delta[dd]/sigma[dd];
    }
    for (dd = 0; dd < dim; dd++) sigmafile >> intcpt[dd];
    for (dd = 0; dd < dim; dd++) sigmafile >> scale[dd];
    sigmafile.get(ch);            while (ch != '\n') sigmafile.get(ch);


    /** Read G-spline weights                              **/
    /** Compute directly weights of marginal distributions **/
    /** Compute also means (knots) in each margin          **/
    if (wfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + wpath + " before "
               + char(ihelp) + std::string(" sets of G-spline weights were read.");
      throw returnR(errmes, 99);
    }
    if (mufile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + mupath + " before "
               + char(ihelp) + std::string(" sets of G-spline means were read.");
      throw returnR(errmes, 99);
    }      
    for (dd = 0; dd < dim; dd++){
      mu_sig_marg[dd][0] = (gamma[dd] - GsplK[dd]*delta[dd]) / sigma[dd];
      w_marg[dd][0]  = 0.0;
      for (j = 1; j < 2*GsplK[dd] + 1; j++){
        w_marg[dd][j]      = 0.0;                                       /* reset marginal     weights             */
        mu_sig_marg[dd][j] = mu_sig_marg[dd][j-1] + delta_sig[dd];      /* compute mean/sigma for each margin      */
      }
    }
    for (j = 0; j < k_current; j++){
      wfile >> tmp;                                        /* read the weight                    */
      for (dd = 0; dd < dim; dd++){
        mufile >> ihelp;                                   /* read index of the component        */
        w_marg[dd][ihelp + GsplK[dd]] += tmp;              /* update appropriate marginal weight */
      }
    }
    wfile.get(ch);                 while (ch != '\n') wfile.get(ch);
    mufile.get(ch);                while (ch != '\n') mufile.get(ch);

    return;    
  }  // end of try
  catch(returnR){
    throw;
  }  
}

// ==================================================================================
// ***** readGsplineFromFiles3: Function to read one sampled G-spline from four files 
//   version used by 'predictive_GS' for the random intercept
// ==================================================================================
//
// !!!!! useful only for univariate G-spline !!!!!!
//
// k_effect ...................... effective number of mixture components read
// cum_w[k_effect] ............... cumulative weights
// prop_mu[k_effect] ............. means corresponding to cumulative weights (intcpt + scale*knot)
//                                 !!! this is suitable only for univariate G-spline !!!
// sig_scale[dim] ................ sigma*scale for each dimension
// skip .......................... number of rows that should be skipped before reading real data
// row ........................... how many rows of the data will be read at the end of the function call
// dim ........................... dimension of the G-spline
// total_length .................. maximal number of components in the G-spline
// GsplK[dim] .................... "K" parameters of the G-spline
// kfile ......................... usually stream associated with 'mixmoment_b*.sim'
// wfile ......................... usually stream associated with 'mweight_b*.sim'
// mufile ........................ usually stream associated with 'mmean_b*.sim'
// sigmafile ..................... usually stream associated with 'gspline_b*.sim'
//
void
readGsplineFromFiles3(int* k_effect,             double* cum_w,             double* prop_mu,            double* sig_scale,
                      const int& skip,           const int& row,            const int& dim,             const int& total_length,   
                      std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,      std::ifstream& sigmafile,
                      const std::string& kpath,  const std::string& wpath,  const std::string& mupath,  const std::string& sigmapath)
{
  try{
    static int j, k_current, ihelp;
    static double tmp, gamma, sigma, intcpt, scale, delta, intcpt_scale_gamma, scale_delta;
    static char ch;
    static std::string errmes;

    if (dim > 1) throw returnR("Error: 'readGsplineFromFiles3' is not implemented for dimension higher than 1", 1);

    /**  Skip rows that are to be skipped  **/
    for (j = 0; j < skip; j++){
      kfile.get(ch);        while (ch != '\n') kfile.get(ch);
      wfile.get(ch);        while (ch != '\n') wfile.get(ch);
      mufile.get(ch);       while (ch != '\n') mufile.get(ch);
      sigmafile.get(ch);    while (ch != '\n') sigmafile.get(ch);
    }

    /**  Read effective k **/
    if (kfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + kpath + " before "
               + char(ihelp) + std::string(" values were read.");
      throw returnR(errmes, 99);
    }
    kfile >> k_current;
    if (k_current > total_length) throw returnR("C++ Error: k value higher than indicated total_length of the G-spline was read.", 99);
    *k_effect = k_current;
    kfile.get(ch);                     
    while (ch != '\n') kfile.get(ch);

    /** Read G-spline intercept, standard deviations and distance between two knots **/
    if (sigmafile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + sigmapath + " before "
               + char(ihelp) + std::string(" sets of G-spline intercepts/std. deviations were read.");
      throw returnR(errmes, 99);
    }
    sigmafile >> gamma;
    sigmafile >> sigma;
    sigmafile >> delta;
    sigmafile >> intcpt;
    sigmafile >> scale;
    sigmafile.get(ch);            
    while (ch != '\n') sigmafile.get(ch);
    sig_scale[0]        = sigma * scale;
    intcpt_scale_gamma  = intcpt + scale * gamma;
    scale_delta         = scale * delta;

    /** Read G-spline weights                              **/
    /** Compute directly cumulative weights                **/
    /** Compute also corresponding propoal means           **/
    if (wfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + wpath + " before " 
               + char(ihelp) + std::string(" sets of G-spline weights were read.");
      throw returnR(errmes, 99);
    }
    if (mufile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + mupath + " before "
               + char(ihelp) + std::string(" sets of G-spline means were read.");
      throw returnR(errmes, 99);
    }      
    wfile >> cum_w[0];
    mufile >> ihelp;
    prop_mu[0] = intcpt_scale_gamma + ihelp*scale_delta;
    for (j = 1; j < k_current; j++){
      wfile >> tmp;
      cum_w[j] = cum_w[j-1] + tmp;
      mufile >> ihelp;
      prop_mu[j] = intcpt_scale_gamma + ihelp*scale_delta;
    }    
    wfile.get(ch);                 
    while (ch != '\n') wfile.get(ch);
    mufile.get(ch);                
    while (ch != '\n') mufile.get(ch);

    return;    
  }  // end of try
  catch(returnR){
    throw;
  }  
}


// ==================================================================================
// ***** closeGsplineFiles
// ==================================================================================
void
closeGsplineFiles(std::ifstream& kfile,  std::ifstream& wfile,  std::ifstream& mufile,  std::ifstream& sigmafile)
{
  kfile.close();
  wfile.close();
  mufile.close();
  sigmafile.close();
  return;
}


// ==============================================================================================================
// ***** openRegresFiles: Open files where sampled betas and D's are stored for reading
//                          and skip first 'skip' rows that the user wishes to skip
//     * used by predictive_GS in the case of either no or normal random effects
// ==============================================================================================================
void
openRegresFiles(std::ifstream& betafile,      std::ifstream& Dfile,
                const std::string& betapath,  const std::string& Dpath,
		const int& skip,              const int& nbeta,          const int& nRandom,   const bool& reff_NORMAL)
{
  try{
    std::string errmes;
    char ch;
    int i;

    if (nbeta){
      betafile.open(betapath.c_str(), std::ios::in);
      if (!betafile){
        errmes = std::string("C++ Error: Could not open ") + betapath;
        throw returnR(errmes, 99);
      } 
      for (i = 0; i < skip; i++){
        betafile.get(ch);        
        while (ch != '\n') betafile.get(ch);
      }
    }

    if (nRandom && reff_NORMAL){
      Dfile.open(Dpath.c_str(), std::ios::in);
      if (!Dfile){
        errmes = std::string("C++ Error: Could not open ") + Dpath;
        throw returnR(errmes, 99);
      } 
      for (i = 0; i < skip; i++){
        Dfile.get(ch);        
        while (ch != '\n') Dfile.get(ch);
      }
    }
    return;
  }
  catch(returnR){
    throw;
  }   
}


// ==============================================================================================================
// ***** openD32File: Open file where sampled D's are stored for reading
//                          and skip first 'skip' rows that the user wishes to skip
//     * used by "version=32" predictive_GS
// ==============================================================================================================
void
openD32File(std::ifstream& D32file,  const std::string& D32path,  const int& skip)
{
  try{
    std::string errmes;
    char ch;
    int i;

    D32file.open(D32path.c_str(), std::ios::in);
    if (!D32file){
      errmes = std::string("C++ Error: Could not open ") + D32path;
      throw returnR(errmes, 99);
    } 
    for (i = 0; i < skip; i++){
      D32file.get(ch);        
      while (ch != '\n') D32file.get(ch);
    }

    return;
  }
  catch(returnR){
    throw;
  }   
}


/*** =========================================================================== ***/
/*** readDfromFile:  Read D matrix from the file                                 ***/
/***                 and update accordingly corresponding elements               ***/
/***                                                                             ***/
/*** =========================================================================== ***/
//
// skip:     Number of rows to be skipped before reading
// row:      How many rows of the data will be read at the end of the function call
// Dfile:    Opened ifstream to the file with D values
//           We are reading from the current row and its assumed that the row contains
//           det(D), D(0,0), D(1,0), D(1,1)
// Dpath:    String with the path to the D file (just to generate error messages in the case they are needed)
//
void
readDfromFile(RandomEff32::RE *data,  const int &skip,  const int &row,  std::ifstream &Dfile,   const std::string &Dpath)
{
  static std::string errmes;
  static double tmp, *DP;
  static int j, ihelp;
  static char ch;
  //const int lD = 3;

  for (j = 0; j < skip; j++){
    Dfile.get(ch);        
    while (ch != '\n') Dfile.get(ch);
  }
  if (Dfile.eof()){
    ihelp = row + 1;
    errmes = std::string("Error: Reached end of file ") + Dpath + " before "
             + char(ihelp) + std::string(" sets of random effects covariance matrices were read.");
    throw returnR(errmes, 99);
  }

  Dfile >> tmp;                /* skip the first column with determinant */
  DP = data->_D;
  for (j = 0; j < data->_lD; j++){
    Dfile >> (*DP);
    DP++;
  }
  RandomEff32::updateAfterChangeD(data);
  Dfile.get(ch);            
  while (ch != '\n') Dfile.get(ch);

  return;
}


// ==================================================================================
// ***** readRegresFromFiles
//  * used in predictive_GS in the case that random effects are NORMAL
// ==================================================================================
void
readRegresFromFiles(BetaGamma* bg,                CovMatrix* DD,
                    const int& skip,              const int& row,
                    std::ifstream& betafile,      std::ifstream& Dfile,
                    const std::string& betapath,  const std::string& Dpath,  const bool& reff_NORMAL)
{
  static int j, ihelp;
  static double tmp;
  static std::string errmes;
  static char ch;

  if (bg->nbeta()){
    for (j = 0; j < skip; j++){
      betafile.get(ch);        
      while (ch != '\n') betafile.get(ch);
    }
    if (betafile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + betapath + " before "
               + char(ihelp) + std::string(" sets of regression parameters were read.");
      throw returnR(errmes, 99);
    }
    for (j = 0; j < bg->nbeta(); j++){
      betafile >> tmp;
      bg->new_beta(j, tmp);
    }
    betafile.get(ch);            
    while (ch != '\n') betafile.get(ch);
  }

  if (bg->nRandom() && reff_NORMAL){
    for (j = 0; j < skip; j++){
      Dfile.get(ch);        
      while (ch != '\n') Dfile.get(ch);
    }
    if (Dfile.eof()){
      ihelp = row + 1;
      errmes = std::string("C++ Error: Reached end of file ") + Dpath + " before "
               + char(ihelp) + std::string(" sets of random effects covariance matrices were read.");
      throw returnR(errmes, 99);
    }
    Dfile >> tmp;                /* skip the first column with determinant */
    for (j = 0; j < DD->larray(); j++){
      Dfile >> tmp;
      DD->new_covm(j, tmp);
    }
    DD->update_after_change_covm();
    Dfile.get(ch);            
    while (ch != '\n') Dfile.get(ch);
  }
  return;
}


// ==================================================================================
// ***** closeRegresFiles
//  * used in predictive_GS in the case that random effects are NORMAL
// ==================================================================================
void
closeRegresFiles(std::ifstream& betafile,  std::ifstream& Dfile,
                 const int& nbeta,         const int& nRandom,    const bool& reff_NORMAL)
{
  if (nbeta) betafile.close();
  if (reff_NORMAL && nRandom) Dfile.close();
  return;
}


// ==================================================================================
// ***** writeToFiles_random
// ==================================================================================
void
writeToFiles_random(
    const CovMatrix* Dm,     const RandomEff* bb, 
    const int* storebP,      const int* writeAll,
    std::ofstream& Dfile,    std::ofstream& bbfile,
    const int& prec,         const int& width)
{

  /*** Sampled covariance matrix of random effects (together with its determinant) ***/
  static double detD;
  detD = Dm->det();
  writeTwoToFile_1(&detD, Dm->covmP(), 0, Dm->larray(), Dfile, prec, width);

  /*** Sampled random effects ***/
  if (*storebP || *writeAll){
    writeToFile_1(bb->bMP(), bb->lbMarray(), bbfile, prec, width);
  }

  return;
}

// ==================================================================================
// ***** writeToFiles_Gspl_intcpt
// ==================================================================================
void
writeToFiles_Gspl_intcpt(
    const RandomEff* bb, 
    const int* storebP,      const int* writeAll,
    std::ofstream& bbfile,
    const int& prec,         const int& width)
{

  /*** Sampled random intercept ***/

  if (*storebP || *writeAll){
    writeToFile_1(bb->bMP(), bb->lbMarray(), bbfile, prec, width);
  }

  return;
}


// ==================================================================================
// ***** writeToFiles_RandomEff32
// ==================================================================================
void
writeToFiles_RandomEff32(
    const RandomEff32::RE *db,
    const int *storedP,         const int *storebP,      const int *writeAll,
    std::ofstream &Dfile,       std::ofstream &ddfile,   std::ofstream &bbfile,
    const int &prec,            const int &width)
{

  /*** Sampled covariance matrix of random effects (together with its determinant) ***/
  writeTwoToFile_1(&db->_detD, db->_D, 0, db->_lD, Dfile, prec, width);

  /*** Sampled random effects ***/
  if (*writeAll){
    writeToFile_1(db->_d, db->_nCluster, ddfile, prec, width);
    writeToFile_1(db->_b, db->_nCluster, bbfile, prec, width);
  }
  else{
    if (*storedP) writeToFile_1(db->_d, db->_nCluster, ddfile, prec, width);
    if (*storebP) writeToFile_1(db->_b, db->_nCluster, bbfile, prec, width);
  }
  return;
}
