// Headers for non-template classes
//
// 22/03/2004: everything put separately into this file
//
#ifndef CLASSES_H
#define CLASSES_H


// ======================================================================================
// ***** covMatrix.cpp: a class to store covariance matrices
// ======================================================================================
class covMatrix
{
 private:
  int _nrow;          // number of rows
  int _larray;        // length of the lower triangle
  int _rank;          // rank of the matrix  

 public:
  double* covm;       // lower triangle of the matrix stored in an array
                      // This points always to an external array!!!

  double* ichicovm;   // lower triangle of the inverse of the Cholesky decomposition of the inverse of the covariance matrix stored in an array
  double* icovm;      // lower triangle of the inverse of this matrix stored in an array
  int* diagI;         // indeces of diagonal elements in above mentioned arrays

  double* qr;         // QR decomposition as returned by the function dqrdc2CPP
  double* qraux;      // further information required to recover the orthogonal part of the QR decomposition
  int* jpvt;          // pivot information from QR decomposition
  double det;         // determinant of the matrix

  covMatrix();

  covMatrix(const covMatrix& cm);  

  covMatrix(double* covmA, const int* nr, const double* tolChol, const double* tolQR);

  ~covMatrix();

  covMatrix& 
  operator=(const covMatrix& mh);

  inline int
  nrow() const { return _nrow;}

  inline int
  larray() const { return _larray;}

  inline int
  rank() const { return _rank;}

  inline bool
  isNull() const { if (_nrow == 0) return true;
                   else            return false;}

  void
  print() const;

  void
  updateFromInv(const double* tolChol, const double* tolQR);

  void
  update(const double* tolChol, const double* tolQR);

};   // end of the class covMatrix


// ======================================================================================
// ***** bblocks.cpp: a class to store information concerning the random effects
// ======================================================================================
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

