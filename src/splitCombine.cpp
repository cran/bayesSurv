// Function to perform split-combine move

// 28/11/2003: start woking on it
// 16/02/2004: rewritten to be used with auxiliary variables
// 13/03/2004: possibility of having a mean of random intercept dependent on a mixture overall mean was considered
//             (no considerable changes in splitCombine move were needed since due to moments-matching,
//             neither overall mixture mean nor overall mixture variance change with split-combine move)
//
#include "bayessurvreg.h"

using namespace std;

// ********** splitCombine **********
//
// PARAMETERS:
//
// acceptedP ............... +1 if the move will be accepted, +0 if the move will be rejected
// splitP .................. +1 if the split move was proposed, +0 if the combine move was proposed
// loglikelhood ............ on INPUT: current value of the log-likelihood, 
//                           on OUTPUT: (new) value of the log-likelihood
// loglikobs ...............
// proploglikobs ...........
// kP ...................... on INPUT: current number of mixture components
//                           on OUTPUT: number of mixture components after the split-combine move
// wM ...................... on INPUT: current mixture weights
//                           on OUTPUT: mixture weights after the split-combine move
// muM ..................... on INPUT: current mixture means
//                           on OUTPUT: mixture means after the split-combine move
// invsigma2M .............. on INPUT: current mixture inverse-variances
//                           on OUTPUT: mixture inverse-variances after the split-combine move
// rM ...................... on INPUT: current allocations of the observations
//                           on OUTPUT: allocations of the observations after the split-combine move
// invrM ................... working space to store proposed allocations
// mixtureNM ............... on INPUT: current numbers of observations at each mixture component
//                           on OUTPUT: numbers of observations at each mixture component after the split-combine move
// propkP .................. working space to store proposed values of k
// propwM .................. working space to store proposed values of w
// propmuM ................. working space to store proposed values of mu
// propinvsigma2M .......... working space to store proposed values of invsigma2
// proprM .................. working space to store proposed values of r
// propinvrM ............... working space to store proposed values of invr
// propmixtureNM ........... working space to store proposed values of mixtureN
// uM ...................... canonical proposal vector of length 3
//                           (pointer to an appropriate part of full uM)
//                           INPUT value from uM[0], uM[1], uM[2] used if the split move is proposed
//                           OUTPUT value returns in uM[-3], uM[-2], uM[-1] computed canonical proposal vector 
//                                        if the combine move is proposed and accepted
// (*logdu) ................ function to compute a log-density of a canonical proposal vector
// (*transu) ............... function used to transform a canonical proposal vector u -> v
// (*invtransu) ............ inverse to that transformation (v -> u
// (*logJtransu) ........... log(Jacobian) of transformation v(u)
// regresResM .............. regression residuals (y - x'beta - z'b)   
// YsM .....................
// Eb0 ..................... mean of the random intercept
//                           (it never changes with split-combine move due to moments matching)
// randomIntP .............. 0/1 indicating a presence of a random intercept
// kmaxP ................... maximal number of mixture components
// piSplitM ................ probabilities of the split move (given k) ((*kmaxP + 1) x 1)
// logpiSplitM ............. log-probabilities of the split move (given k) ((*kmaxP + 1) x 1)
// logpiCombineM ........... log-probabilities of the combine move (given k) ((*kmaxP + 1) x 1)
// deltaP .................. hyperparameter
// xiP ..................... hyperparameter
// invkappaP ............... (transformed) hyperparameter
// halfl2pikappaP .......... (transformed) hyperparameter
// zetaP ................... hyperparameter
// etaP .................... hyperparameter
// lgammazetaP ............. (transformed) hyperparameter
// llambdaP ................ (transformed) hyperparameter
// priorParmu .............. parameters of the prior distribution of the canonical proposal vector u
// transParmu .............. parameters for transformation of a canonical proposal vector u
// priorForkP .............. indicator of the prior for k (0 = Poisson, 1 = uniform)
// (*logdtrans). ........... log of the Jacobian of the transformation of the response
// nP ...................... number of observations
//
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
             const int* priorForkP,         double (*logdtrans) (double), const int* nP)
{
  if (*priorForkP == Fixed) return;

  int errorType = Mixture;

  int i, j;
  double unif;
  int accepted = 1;
  int split = 0;
  double vM[3];           // transformed proposal vector
  double propuM[3];       // computed proposed canonical proposal vector (used in a combine move)

  // Determine whether split or combine move will be proposed
  // ========================================================
  unif = runif(0, 1);
  if (unif < piSplitM[*kP]){
    split = 1;
    *splitP += 1;
  }

  //  cout << (split ? "SPLIT" : "COMBINE") << endl;

  int jsplit;           // index of the component which is to be splitted
                        // OR the index of the first component which is to be combined with the component jsplit + 1
  double logPalloc;               // logarithm of the allocation probability
  double logJacobian;             // log(Jacobian) of the split move (for combine move, it must be multiplied by -1)
  double logPostRatio;            // log(posterior ratio) of the move
  double proploglikelhood;        // log-likelihood at the proposed parameters
  double logpiu;                  // log-density of the canonical proposal vector
  double logdtransu;              // log-Jacobian of the transformation u -> v
  double Paccept;                 // acceptance probability
  double logpiComSpl;             // log(pi_{k+1}(combine) / pi_k(split)) OR log(pi_{k-1}(split) / pi_k(combine))
  
  // Split move
  // ============
  if (split){
    *propkP = *kP + 1;    
    logpiComSpl = logpiCombineM[*propkP] - logpiSplitM[*kP];

    (*transu)(vM, uM, transParmu);    

    discreteUniformSampler(&jsplit, kP, &ONE_INT, &ZERO_INT);

    proposeSplit(&accepted, *propwM, *propmuM, *propinvsigma2M, *wM, *muM, *invsigma2M, vM, jsplit, kP);
    if (!accepted) return;

    logPalloc = allocSplit(proprM, *propinvrM, *propmixtureNM, rM, *invrM, *mixtureNM, 
                           *propwM, *propmuM, *propinvsigma2M, jsplit, kP, regresResM, Eb0, randomIntP);

    logJacobian = logJacobianSplitCombine((*wM)[jsplit], (*propmuM)[jsplit], (*propmuM)[jsplit + 1],
                                          (*propinvsigma2M)[jsplit], (*propinvsigma2M)[jsplit + 1], (*invsigma2M)[jsplit], vM);

    logLikelihood(&proploglikelhood, *proploglikobs, nP, regresResM, YsM, 
                  propkP, proprM, *propwM, *propmuM, *propinvsigma2M, Eb0, &errorType, randomIntP, logdtrans);

    logPostRatio = proploglikelhood - (*loglikelhood);
    logPostRatio += logPostRatioSplitCombine(jsplit, kP, 
                                *propwM, *wM, *propmuM, *muM, *propinvsigma2M, *invsigma2M, *propmixtureNM, *mixtureNM,
                                deltaP, xiP, invkappaP, halfl2pikappaP, zetaP, etaP, lgammazetaP, llambdaP, priorForkP);

    
    logpiu = (*logdu)(uM, priorParmu);
    logdtransu = (*logJtransu)(uM, vM, transParmu);

    Paccept = exp(logPostRatio + logpiComSpl - logPalloc - logpiu + logdtransu + logJacobian);

  }   // end of the split move
   
  // Combine move
  // =============
  else{
    *propkP = *kP - 1;    
    logpiComSpl = logpiSplitM[*propkP] - logpiCombineM[*kP];

    discreteUniformSampler(&jsplit, propkP, &ONE_INT, &ZERO_INT);

    proposeCombine(&accepted, vM, *propwM, *propmuM, *propinvsigma2M, *wM, *muM, *invsigma2M, jsplit, propkP);

    (*invtransu)(propuM, vM, transParmu);

    logPalloc = allocCombine(proprM, *propinvrM, *propmixtureNM, rM, *invrM, *mixtureNM, 
                             *wM, *muM, *invsigma2M, jsplit, propkP, regresResM, Eb0, randomIntP);

    logJacobian = logJacobianSplitCombine((*propwM)[jsplit], (*muM)[jsplit], (*muM)[jsplit + 1],
                                          (*invsigma2M)[jsplit], (*invsigma2M)[jsplit + 1], (*propinvsigma2M)[jsplit], vM);

    logLikelihood(&proploglikelhood, *proploglikobs, nP, regresResM, YsM, 
                  propkP, proprM, *propwM, *propmuM, *propinvsigma2M, Eb0, &errorType, randomIntP, logdtrans);

    logPostRatio = (*loglikelhood) - proploglikelhood;
    logPostRatio += logPostRatioSplitCombine(jsplit, propkP, 
                                *wM, *propwM, *muM, *propmuM, *invsigma2M, *propinvsigma2M, *mixtureNM, *propmixtureNM,
                                deltaP, xiP, invkappaP, halfl2pikappaP, zetaP, etaP, lgammazetaP, llambdaP, priorForkP);

    logpiu = (*logdu)(propuM, priorParmu);
    logdtransu = (*logJtransu)(propuM, vM, transParmu);

    Paccept = exp(-logPostRatio + logpiComSpl + logPalloc + logpiu - logdtransu - logJacobian);

  }   // end of combine move

  // Accept or reject the proposal
  // =============================
  if (Paccept < 1.0){
    unif = runif(0, 1);
    if (unif > Paccept){
      return;
    }  
  }

    // Proposal was accepted
  *acceptedP += 1;
  *kP = *propkP;
  *loglikelhood = proploglikelhood;
  changePointers(loglikobs, proploglikobs);
  changePointers(wM, propwM);
  changePointers(muM, propmuM);
  changePointers(invsigma2M, propinvsigma2M);
  for (i = 0; i < *nP; i++){ rM[i] = proprM[i];}     // this cannot be done using changePointers!!! Since rM was not created using new.
  changePointers(invrM, propinvrM);
  changePointers(mixtureNM, propmixtureNM);
  if (!split){
    for (i = 0; i < 3; i++) uM[i - 3] = propuM[i];
  }

  return;    
}   // end of function splitCombine


// ********** proposeSplit **********
//
// Compute the proposed values of w, mu, invsigma2 for the split move
//
// INPUT PARAMETERS:
// 
// wM
// muM ............ w, mu, invsigma2 before the split move is proposed
// invsigma2M
// vM ............. proposal vector of length 3
// jsplit ......... index of the splitted component 
// kP ............. number of components before the split move was proposed
//
// OUTPUT PARAMETERS:
//
// acceptedP ...... it is set to 0 if the adjacency condition is not satisfied
// propwM 
// propmuM ........ proposed values of w, mu, invsigma2
// propinvsigma2M
//
void
proposeSplit(int* acceptedP,
             double* propwM,      double* propmuM,     double* propinvsigma2M,
             const double* wM,    const double* muM,   const double* invsigma2M,
             const double* vM,    const int jsplit,    const int* kP)
{
    int j;

        // Reject directly if we are trying to split a component with zero weight
        //   (such component should never appear...)
    if (wM[jsplit] <= 0.0){
      *acceptedP = 0;
      return;
    }

        // Shift the end of the vectors w, mu, invsigma2
    for (j = *kP; j > jsplit + 1; j--){
        propwM[j] = wM[j - 1];
        propmuM[j] = muM[j - 1];
        propinvsigma2M[j] = invsigma2M[j - 1];
    }    

        // Compute proposed means
    propmuM[jsplit] = muM[jsplit] - vM[1] 
                         * sqrt((1-vM[0])/(vM[0]*(invsigma2M[jsplit])));
    propmuM[jsplit + 1] = muM[jsplit] + vM[1] 
                             * sqrt((vM[0])/((1-vM[0])*(invsigma2M[jsplit])));

        // Check the adjacency condition
    if ((jsplit > 0) && (propmuM[jsplit] <= muM[jsplit - 1])) *acceptedP = 0;
    if ((jsplit < *kP - 1) && (propmuM[jsplit + 1] >= muM[jsplit + 1])) *acceptedP = 0;
    if (!(*acceptedP)) return;

        // Compute proposed weights and variances
    propwM[jsplit] = wM[jsplit] * vM[0];
    propwM[jsplit + 1] = wM[jsplit] * (1 - vM[0]);
    double hulp = 1 / (1 - vM[1]*vM[1]);
    propinvsigma2M[jsplit] = invsigma2M[jsplit] * (vM[0] / vM[2]) * hulp;
    propinvsigma2M[jsplit + 1] = invsigma2M[jsplit] * ((1 - vM[0]) / (1 - vM[2])) * hulp;

        // Keep the beginning of the vectors w, mu, invsigma2
    for (j = jsplit - 1; j >= 0; j--){
        propwM[j] = wM[j];
        propmuM[j] = muM[j];
        propinvsigma2M[j] = invsigma2M[j];
    }

        // Check whether proposed weights are still computationally bigger than 0
        // otherwise, reject this move
    if (propwM[jsplit] <= 0.0 || propwM[jsplit + 1] <= 0.0) *acceptedP = 0;

    return;
}   // end of function proposeSplit


// ********** proposeCombine **********
//
// Compute the proposed values of w, mu, invsigma2 for the combine move
//
// INPUT PARAMETERS:
// 
// wM
// muM ............ w, mu, invsigma2 before the split move is proposed
// invsigma2M
// jsplit ......... index of the first combined component 
// propkP ......... number of components after the combine move
//
// OUTPUT PARAMETERS:
//
// vM ............. proposal vector of length 3
// propwM 
// propmuM ........ proposed values of w, mu, invsigma2
// propinvsigma2M
//
void
proposeCombine(int* acceptedP,
               double* vM,
               double* propwM,    double* propmuM,     double* propinvsigma2M,
               const double* wM,  const double* muM,   const double* invsigma2M,
               const int jsplit,  const int* propkP)
{
    int j;
 
        // Reject directly if at least one combined component has computationally zero weight
        //   (this should never happens)
    if (wM[jsplit] <= 0.0 || wM[jsplit + 1] <= 0.0){
      *acceptedP = 0;
      return;
    }

        // Keep the beginning of the vectors w, mu, invsigma2
    for (j = 0; j < jsplit; j++){
        propwM[j] = wM[j];
        propmuM[j] = muM[j];
        propinvsigma2M[j] = invsigma2M[j];
    }

        // Compute proposed weights, means, variances
    propwM[jsplit] = wM[jsplit] + wM[jsplit + 1];
    propmuM[jsplit] = (wM[jsplit] * muM[jsplit] + wM[jsplit + 1] * muM[jsplit + 1])/(propwM[jsplit]);
    propinvsigma2M[jsplit] = (propwM[jsplit]) / (wM[jsplit]*(muM[jsplit] * muM[jsplit] + (1/invsigma2M[jsplit]))
                                         + wM[jsplit+1]*(muM[jsplit+1] * muM[jsplit+1] + (1/invsigma2M[jsplit+1]))
                                         - propwM[jsplit]*propmuM[jsplit]*propmuM[jsplit]);

        // Shift the end of the vectors w, mu, invsigma2
    for (j = jsplit + 1; j < *propkP; j++){
        propwM[j] = wM[j + 1];
        propmuM[j] = muM[j + 1];
        propinvsigma2M[j] = invsigma2M[j + 1];
    }    
    propwM[*propkP] = 0.0;
    propmuM[*propkP] = 0.0;
    propinvsigma2M[*propkP] = 0.0;

       // Compute appropriate values of the corresponding proposal vector
    vM[0] = wM[jsplit]/propwM[jsplit];
    vM[1] = sqrt(propinvsigma2M[jsplit] * (wM[jsplit]/wM[jsplit+1])) * (propmuM[jsplit] - muM[jsplit]);
    vM[2] = (propinvsigma2M[jsplit]/invsigma2M[jsplit]) * (vM[0]/(1 - vM[1] * vM[1]));

    return; 
}   // end of function proposeCombine


// ********** allocSplit **********
//
// Propose reallocations to the splitted components when the split move was proposed
//   and compute a logarithm of the allocation probability
//
// RETURN:
//   logarithm of the allocation probability
//   * it returns -FLT_MAX if that probability is zero
//
// INPUT PARAMETERS:
//
// proprM
// propinvrM
// propmixtureNM
// jsplit ......... index of the splitted component 
// kP ............. number of components before the split move was proposed
//
// OUTPUT PARAMETERS:
//
// rM
// invrM
// mixtureNM
//
double
allocSplit(int* proprM,               List<int>* propinvrM,    int* propmixtureNM,
           const int* rM,             const List<int>* invrM,  const  int* mixtureNM,
           const double* propwM,      const double* propmuM,   const double* propinvsigma2M,
           const int jsplit,          const int* kP,           
           const double* regresResM,  const double* Eb0,       const int* randomIntP)
{ 
    int i, j, obs;
    double logPalloc = 0.0;
    double intcptadd = (*randomIntP ? (*Eb0) : 0.0);    

      // Shift rM to proprM (only observations that do not belong to a splitted component)
    for (j = *kP; j > jsplit + 1; j--){
      propinvrM[j] = invrM[j - 1]; 
      for (i = 0; i < invrM[j - 1].length(); i++){
        proprM[invrM[j - 1][i]] = j;
      }
      propmixtureNM[j] = mixtureNM[j - 1];
    }    
    for (j = jsplit - 1; j >= 0; j--){
      propinvrM[j] = invrM[j]; 
      for (i = 0; i < invrM[j].length(); i++){
        proprM[invrM[j][i]] = j;
      }
      propmixtureNM[j] = mixtureNM[j];
    }

    propinvrM[jsplit] = List<int>();
    propinvrM[jsplit + 1] = List<int>();
    propmixtureNM[jsplit] = 0;
    propmixtureNM[jsplit + 1] = 0;

      // Propose reallocations
    if (!invrM[jsplit].isEmpty()){      // Are we splitting a  non-empty component?

      double unif;
      double probs[2];
      double sumprobs;
      double winvsigmaM[2];
      winvsigmaM[0] = propwM[jsplit] * sqrt(propinvsigma2M[jsplit]);
      winvsigmaM[1] = propwM[jsplit + 1] * sqrt(propinvsigma2M[jsplit + 1]);

      for (i = 0; i < invrM[jsplit].length(); i++){
        obs = invrM[jsplit][i];
        probs[0] = winvsigmaM[0] * exp(-0.5 * propinvsigma2M[jsplit] * 
                      (regresResM[obs] - propmuM[jsplit] + intcptadd)*(regresResM[obs] - propmuM[jsplit] + intcptadd));
        probs[1] = winvsigmaM[1] * exp(-0.5 * propinvsigma2M[jsplit + 1] * 
                      (regresResM[obs] - propmuM[jsplit + 1] + intcptadd)*(regresResM[obs] - propmuM[jsplit + 1] + intcptadd));
        sumprobs = probs[0] + probs[1];
        unif = runif(0, sumprobs);
        if (unif < probs[0]){          // component jsplit
          probs[0] /= sumprobs;
          if (probs[0] <= 0.0) logPalloc = -FLT_MAX;
          else                 logPalloc += log(probs[0]);
          proprM[obs] = jsplit;
          propmixtureNM[jsplit]++;
          propinvrM[jsplit].addNode(obs);
        }
        else{                       // component jsplit + 1
          probs[1] /= sumprobs;
          if (probs[1] <= 0.0) logPalloc = -FLT_MAX;
          else                 logPalloc += log(probs[1]);
          proprM[obs] = jsplit + 1;          
          propmixtureNM[jsplit + 1]++;
          propinvrM[jsplit + 1].addNode(obs);
        }
      }
    }

    return logPalloc;
}   // end of function allocSplit


// ********** allocCombine **********
//
// Give allocations for observations belonging to the combined component
//   and compute a logarithm of the allocation probability
//
// RETURN:
//   logarithm of the allocation probability
//   * it returns -FLT_MAX if that probability is zero
//
// INPUT PARAMETERS:
//
// proprM
// propinvrM
// propmixtureNM
// jsplit ......... index of the first combined component
// propkP ..........number of components after the combined move
//
// OUTPUT PARAMETERS:
//
// rM
// invrM
// mixtureNM
//
double
allocCombine(int* proprM,               List<int>* propinvrM,    int* propmixtureNM,
             const int* rM,             const List<int>* invrM,  const int* mixtureNM,
             const double* wM,          const double* muM,       const double* invsigma2M,
             const int jsplit,          const int* propkP,
             const double* regresResM,  const double* Eb0,       const int* randomIntP)
{
    int i, j, obs;
    double logPalloc = 0.0;
    double intcptadd = (*randomIntP ? (*Eb0) : 0.0);    
 
      // Shift rM to proprM (only observations that do not belong to combined component)
    for (j = 0; j < jsplit; j++){
      propinvrM[j] = invrM[j]; 
      for (i = 0; i < invrM[j].length(); i++){
        proprM[invrM[j][i]] = j;
      }
      propmixtureNM[j] = mixtureNM[j];
    }
    for (j = jsplit + 1; j < *propkP; j++){
      propinvrM[j] = invrM[j + 1]; 
      for (i = 0; i < invrM[j + 1].length(); i++){
        proprM[invrM[j + 1][i]] = j;
      }
      propmixtureNM[j] = mixtureNM[j + 1];
    }    
    propinvrM[jsplit] = List<int>();
    propmixtureNM[jsplit] = 0;
    propinvrM[*propkP] = List<int>();
    propmixtureNM[*propkP] = 0;


      // Change allocation of the observations from the combined components
    if (invrM[jsplit].length() || invrM[jsplit + 1].length()){    // Will be the combined component non-empty?

      double unif;
      double probs[2];
      double sumprobs;
      double winvsigmaM[2];
      winvsigmaM[0] = wM[jsplit] * sqrt(invsigma2M[jsplit]);
      winvsigmaM[1] = wM[jsplit + 1] * sqrt(invsigma2M[jsplit + 1]);

      for (int js = 0; js <= 1; js++){
        for (i = 0; i < invrM[js + jsplit].length(); i++){
          obs = invrM[js + jsplit][i];
          probs[0] = winvsigmaM[0] * exp(-0.5 * invsigma2M[jsplit] * 
                        (regresResM[obs] - muM[jsplit] + intcptadd)*(regresResM[obs] - muM[jsplit] + intcptadd));
          probs[1] = winvsigmaM[1] * exp(-0.5 * invsigma2M[js + 1] * 
                        (regresResM[obs] - muM[jsplit + 1] + intcptadd)*(regresResM[obs] - muM[jsplit + 1] + intcptadd));
          sumprobs = probs[0] + probs[1];
          probs[js] /= sumprobs;
          if (probs[js] <= 0.0) logPalloc = -FLT_MAX;
          else                  logPalloc += log(probs[js]);
          proprM[obs] = jsplit;
          propmixtureNM[jsplit]++;
          propinvrM[jsplit].addNode(obs);
        }
      }
    }

    return logPalloc;
}   // end of function allocCombine


// ********** logPostRatioSplitCombine **********
//
// Compute a logarithm of the ratio of posterior that appear in the acceptance probability
//   of the split-combine move.
// The result corresponds to the split move, so that it must be multiplied by (-1) to be used
//   with the combine move.  
//
// INPUT PARAMETERS:
//
// jsplit .................. either index of the splitted component 
//                           or index of the first combined component
// shortkP ................. number of components in the shorter state, i.e.
//                           with split move, *shortkP = # components before the split move was proposed
//                           with combine move, *shortkP = # components after the combine move
//
// longwM
// shortwM     
// longmuM
// shortmuM ................ appropriate parameters with the same convention for 'short' and 'long' as above
// longinvsigma2M
// shortinvsigma2M
// longmixtureNM
// shortmixtureNM
//                           
// deltaP .................. hyperparameter
// xiP ..................... hyperparameter
// invkappaP ............... transformed hyperparameter
// halfl2pikappaP .......... (transformed) hyperparameter
// zetaP ................... hyperparameter
// etaP .................... hyperparameter
// lgammazetaP ............. (transformed) hyperparameter
// llambdaP ................ hyperparameter
// priorForkP .............. indicator of the prior for k (0 = Poisson, 1 = uniform)
//
double 
logPostRatioSplitCombine(const int jsplit,              const int* shortkP,
                         const double* longwM,          const double* shortwM,     
                         const double* longmuM,         const double* shortmuM,
                         const double* longinvsigma2M,  const double* shortinvsigma2M,
                         const int* longmixtureNM,      const int* shortmixtureNM,
                         const double* deltaP,          const double* xiP,              const double* invkappaP,                   
                         const double* halfl2pikappaP,  const double* zetaP,            const double* etaP,                        
                         const double* lgammazetaP,     const double* llambdaP,         const int* priorForkP)
{
  double logPostRatio = 0.0;

    // log(ratio of priors of allocations r_{i,l}: r_{i,l} = j, j1, j2  and weights w_j, w_{j1}, w_{j2})
  logPostRatio += (longmixtureNM[jsplit] + (*deltaP) - 1) * log(longwM[jsplit]) 
                      + (longmixtureNM[jsplit + 1] + (*deltaP) - 1) * log(longwM[jsplit + 1]) 
                      - (shortmixtureNM[jsplit] + (*deltaP) - 1) * log(shortwM[jsplit]) 
                      - lbeta(*deltaP, (*shortkP) * (*deltaP));
    

    // log(ratio of priors of means mu_j, mu_{j1}, mu_{j2})
    //  (without adding log(*shortkP + 1)!!!)
  logPostRatio += -(*halfl2pikappaP) 
                  - 0.5*(*invkappaP) * ((longmuM[jsplit] - (*xiP)) * (longmuM[jsplit] - (*xiP)) + 
                                       (longmuM[jsplit + 1] - (*xiP)) * (longmuM[jsplit + 1] - (*xiP)) -
                                       (shortmuM[jsplit] - (*xiP)) * (shortmuM[jsplit] - (*xiP)));

    // log(ratio of priors of variances sigma2_j, sigma2_{j1}, sigma2_{j2})
  if (*etaP <= 0.0) return -FLT_MAX;
  logPostRatio += (*zetaP) * log(*etaP) - (*lgammazetaP) 
                  + (*zetaP + 1) * (log(longinvsigma2M[jsplit]) + log(longinvsigma2M[jsplit + 1]) 
                                    - log(shortinvsigma2M[jsplit])) 
                  - (*etaP) * (longinvsigma2M[jsplit] + longinvsigma2M[jsplit + 1]
                               - shortinvsigma2M[jsplit]);

    // log(ratio of priors of the numbers of components k)
  switch (*priorForkP){
    case Poisson:                // Truncated Poisson prior for k
      logPostRatio += *llambdaP;     // I do not subtract log(*kP + 1)!!!
      break;                         // It cancels with the term I did not add few lines above
    case Uniform:                // Uniform prior for k
      logPostRatio += log(double(*shortkP + 1));  // This is actually the term I did not add in the command where
      break;                                      // log(ratio of priors of means) was computed
  }

  return logPostRatio;

}   // end of function logPostRatioSplitCombine


// ********** logJacobianSplitCombine **********
//
// Compute the log-Jacobian of the split-combine move
//
// INPUT PARAMETERS:
//
// w ............ mixture weight of either the splitted component
//                                  or     the component that originates from the combine move
// mu0 .......... mean of either the first component that originates from the split move
//                        or     the first combined component
// mu1 .......... mean of either the second component that originates from the split move
//                        or     the second combined component
// invsigma20 ... inverse variance of either the first component that originates from the split move
//                                    or     the first combined component
// invsigma21 ... inverse variance of either the second component that originates from the split move
//                                    or     the second combined component
// invsigma2 .... inverse variance of either the splitted component
//                                    or     the component that originates from the combine move
// vM ........... value of the three-dimensional proposal vector 
//
double
logJacobianSplitCombine(const double w,
                        const double mu0,           const double mu1,
                        const double invsigma20,    const double invsigma21,
                        const double invsigma2,     const double* vM)
{
  double logJacobian;

  if (fabs(vM[1]) > ZERO){
    logJacobian = log(w) + log(invsigma2) + log(fabs(mu1 - mu0)) 
                  - log(invsigma20) - log(invsigma21) - log(vM[1]) - log(1 - vM[1] * vM[1])
                  - log(vM[2]) - log(1 - vM[2]);
  }
  else{
    double hulp0, hulp1;
    if(vM[0] < 0.5){
      hulp0 = sqrt(vM[0] / (1 - vM[0]));
      hulp1 = 1 / hulp0;
    }
    else{
      hulp1 = sqrt((1 - vM[0]) / vM[0]);
      hulp0 = 1 / hulp1;
    }
    
    logJacobian = log(w) + 0.5 * log(invsigma2) + log(fabs(hulp0 - hulp1)) 
                  - log(invsigma20) - log(invsigma21) - log(1 - vM[1] * vM[1])
                  - log(vM[2]) - log(1 - vM[2]);  
  } 

  return logJacobian;

}   // end of function logJacobianSplitCombine
