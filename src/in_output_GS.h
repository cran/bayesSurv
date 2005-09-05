// =======================================================================================================
// ***** Header for in_output_GS.cpp: input/output applications for G-spline applications          *****//
// =======================================================================================================
#ifndef _IN_OUTPUT_GS_H_
#define _IN_OUTPUT_GS_H_

#include <R.h>
//#include <Rmath.h>

#include <iostream>
#include <iomanip>
#include <string>

#include "AK_Error.h"
#include "openFile.h"
#include "Gspline.h"
#include "in_output.h"
#include "templatefun.h"
#include "templatefun_GS.h"
#include "classBetaGamma.h"
#include "classCovMatrix.h"
#include "classRandomEff.h"

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
    const int& n_censored,             const char& write_flag);

void
closeFiles_bayesHistogram(
    std::ofstream& sigmafile,          std::ofstream& lambdafile,
    std::ofstream& mixmomentfile,      std::ofstream& mweightfile,      std::ofstream& mlogweightfile,
    std::ofstream& mmeanfile,          std::ofstream& Yfile,            std::ofstream& rfile,              
    std::ofstream& logposterfile,
    const int& n_censored);

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
    const double& null_weight,         const int& prec,                 const int& width);

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
    const int* storeyP,     const int* storerP);

void
writeToFiles2_bayesHistogram(
    const int* iterA,        const int* k_effectA,    const double* momentsA,
    const double* weightsA,  const int* muA,          const double* gamma_sigma_deltaA,
    const double* lambdaA,   const int* rA,           const double* YsA,
    const int& nstored,      const std::string& dir,  const char& write_flag,
    const int& l_momentsA,
    const int& l_weightsA,   const int& l_muA,        const int& l_gamma_sigma_deltaA,
    const int& l_lambdaA,    const int* dimP,         const int* nP,
    const int* storeyP,      const int* storerP,      const int* n_censored);

void
openGsplineFiles(std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,      std::ifstream& sigmafile,
                 const std::string& kpath,  const std::string& wpath,  const std::string& mupath,  const std::string& sigmapath,
                 const int& skip);

void
open_File_toRead(std::ifstream& file,  const std::string& path,  const int& skip);

void
readGsplineFromFiles(int* k_effect,             double* w,                 int** ind_mu,               double** mu,
                     double* gamma,             double* sigma,             double* delta,
                     double* intcpt,            double* scale,
                     const int& skip,           const int& row,            const int& dim,             const int& total_length,
                     std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,      std::ifstream& sigmafile,
                     const std::string& kpath,  const std::string& wpath,  const std::string& mupath,  const std::string& sigmapath);

void
readMean_and_Scale(double* E_gx,                  double* sd_gx,
                   const int& skip,               const int& row,  
                   const int& dim,
                   std::ifstream& mixmomentfile,  const std::string& mixmomentpath);

void
adjust_intercept(double* intcpt,        const int* version,
                 double* E_gx,
                 const int& skip,       const int& row,
                 std::ifstream& file,   const std::string& path);

void
readGsplineFromFiles2(int* k_effect,             double** w_marg,           double** mu_sig_marg,
                      double* gamma,             double* sigma,             double* delta,
                      double* intcpt,            double* scale,             double* delta_sig,
                      const int& skip,           const int& row,            const int& dim,             
                      const int& total_length,   const int* GsplK,
                      std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,      std::ifstream& sigmafile,
                      const std::string& kpath,  const std::string& wpath,  const std::string& mupath,  const std::string& sigmapath);

void
readGsplineFromFiles3(int* k_effect,             double* cum_w,             double* prop_mu,            double* sig_scale,
                      const int& skip,           const int& row,            const int& dim,             const int& total_length,
                      std::ifstream& kfile,      std::ifstream& wfile,      std::ifstream& mufile,      std::ifstream& sigmafile,
                      const std::string& kpath,  const std::string& wpath,  const std::string& mupath,  const std::string& sigmapath);

void
closeGsplineFiles(std::ifstream& kfile,  std::ifstream& wfile,  std::ifstream& mufile,  std::ifstream& sigmafile);

void
openRegresFiles(std::ifstream& betafile,      std::ifstream& Dfile,
                const std::string& betapath,  const std::string& Dpath,
		const int& skip,              const int& nbeta,          const int& nRandom,   const bool& reff_NORMAL);

void
readRegresFromFiles(BetaGamma* bg,                CovMatrix* DD,
                    const int& skip,              const int& row,
                    std::ifstream& betafile,      std::ifstream& Dfile,
                    const std::string& betapath,  const std::string& Dpath,  const bool& reff_NORMAL);

void
closeRegresFiles(std::ifstream& betafile,  std::ifstream& Dfile,
                 const int& nbeta,         const int& nRandom,    const bool& reff_NORMAL);

void
writeToFiles_random(
    const CovMatrix* Dm,     const RandomEff* bb, 
    const int* storebP,      const int* writeAll,
    std::ofstream& Dfile,    std::ofstream& bbfile,
    const int& prec,         const int& width);

void
writeToFiles_Gspl_intcpt(
    const RandomEff* bb, 
    const int* storebP,      const int* writeAll,
    std::ofstream& bbfile,
    const int& prec,         const int& width);

#endif
