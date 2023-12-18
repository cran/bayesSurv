// Headers for MHblocks.cpp
//
// 22/03/2004: everything put separately into this file
//
#ifndef MH_BLOCKS_H_
#define MH_BLOCKS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "constants.h"
#include "cholesky.h"
#include "templatefun.h"
#include "printArray.h"

// ======================================================================================
// ***** MHblocks.cpp: a class to store information concerning the parameters that
//                     will be possibly updated in blocks using a Metropolis-Hastings step
// =======================================================================================
class MHblocks
{
  private:
    int _nBlocks;            // number of blocks
    int _nParams;            // number of parameters in all blocks
    int _maxnInBlock;        // size of the biggest block
    int _isdprior;

  public:
    double* par;            // EXT array of length _nParams with current values of sampled params
    double* proppar;        // %%array of length _nParams to store a proposal during the MH step
                            //   outside the MH steps, this should be an exact copy of 'par'!!!!!

    double* meanpar;        // EXT array of length _nParams with means of up to now sampled values
    double* halfRangeUnif;  // EXT array of length _nParams with specification of a uniform component of the proposal

    double* priorMean;      // EXT array of length _nParams with prior means for each parameter
    double* priorSD;        // EXT array of length _nParams with prior standard deviations for each parameter 
    double* priorInvVar;    // %%array of length _nParams with prior inverse-variances for each parameter

    int* typeUpd;           // EXT array of length _nBlocks giving the type of update (Gibbs = 2, AM = 1, MH = 0)
    int* nInBlock;          // EXT array of length _nBlocks giving a number of parameters within each block
    int* nRandomB;          // %%array of length _nBlocks giving a number of means of random effects in given block
    int* nFixedB;           // %%array of length _nBlocks giving a number of fixed regr. params in a given block

    int** indBlock;         // %%array of length _nBlocks of EXT arrays, each of length nInBlock giving the indeces
                            //   of parameters from the block in the total parameter vector
    int** diagI;            // %%array of length _nBlocks of %%arrays with diagonal indeces of a diagonal 
    double** covpar;        // %%array of length _nBlocks of EXT arrays with a proposal covariance matrix itself
    double** chcovpar;      // %%array of length _nBlocks of %%arrays with a Cholesky decomposition of a proposal cov. matrix 
    double* logdprior;      // %%array of length _nBlocks with log prior density evaluated at each block
    double* weightUnif;     // EXT array of length _nBlocks with weights of the uniform component in the proposal for each block
    double* eps;            // EXT array of length _nBlocks with epsilon numbers from the (AM) algorithm for each block 
    double* sdNum;          // array of length _nBlocks with s(d) numbers from the (AM) algorithm for each block
    int* sumAccept;         // EXT array of length _nBlocks with a number of accepted proposals for each block

    MHblocks();

    MHblocks(const MHblocks& mh);

    MHblocks(double* parP,             int* parmI,        double* parmD,                   
             int* sumAcceptP,          const int* indbA,
             const double* tolerChol,  const int* iter,   const int* wantdprior);

    ~MHblocks();
 
    MHblocks& 
    operator=(const MHblocks& mh);

    inline int
    nBlock() const { return _nBlocks;}

    inline int
    nParam() const { return _nParams;}

    void
    print() const;

};    // end of the class MHblocks


#endif

