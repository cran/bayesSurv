// ==========================================================================================
// ***** Header for bblocks.cpp: a class to store information concerning the random effects
// ==========================================================================================
#ifndef _B_BLOCKS_H_
#define _B_BLOCKS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "constants.h"
#include "cholesky.h"
#include "templatefun.h"
#include "printArray.h"

class bblocks
{
  private:
    int _nRandom;    // number of random effects
    int _nBlocks;    // number of blocks
    int _nCluster;   // number of clusters

  public:
    double* bM;      // pointer to an EXTernal array, current values of random effects  (_nRandom * _nCluster)
    int* priorD;     // pointer to an EXTernal array, prior for matrix D                           (1)
    double* dfD;     // pointer to an EXTernal array, prior df for matrix D if Wishart prior used  (1)
                     //   or 1/range^2 if SDUniform prior used
    double* scaleD;  // pointer to an EXTernal array, prior scale matrix for D if Wishart prior used  (1)
                     //   or upper limit if SDUniform prior used
    int* typeUpd;    // pointer to an EXTernal array, type of update                    (1)
    int* nInBlock;   // pointer to an EXTernal array,                                   (_nBlocks)
 
    int** indBlock;        // array of length _nBlocks with pointers to EXTernal arrays with indeces
                           //   of random effects from a current block
    int** diagI;           // array of length _nBlocks with pointers to arrays with diagonal indeces
                           //   of covariance matrices stored in arrays
    double** covpar;       // array of length _nBlocks with pointers to EXTernal arrays
                           //   where a covariance matrix of a proposal is stored
    double** chcovpar;     // array of length _nBlocks with pointers to arrays
                           //   where a Cholesky decompositions of the proposal covariance matrices are stored
                           //   or only working space

    double* weightUnif;    // pointer to an EXTernal array with weights of a uniform component
                           //   for MH proposal for each block                          (_nBlocks)
    double* halfRangeUnif; // pointer to an EXTernal array with 0.5range of a uniform component
                           //   for MH proposal for each paramater                      (_nRandom)
    int* sumAccept;        // pointer to an EXTernal array with a number of accepted moves
                           //                                                           (_nBlocks * _nCluster)


    bblocks();

    bblocks(const bblocks& bb);

    bblocks(double* bbM,             int* parmI,   double* parmD,
            int* sumAcceptP, 
            const double* tolerChol);

    ~bblocks();

    bblocks&
    operator=(const bblocks& bb);   

    inline int
    nRandom() const { return _nRandom;}

    inline int
    nBlock() const { return _nBlocks;}

    inline int
    nCluster() const { return _nCluster;}

    void
    print() const;

};   // end of the class bblocks

#endif
