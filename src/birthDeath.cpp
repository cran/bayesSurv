// Function to perform birth-death move

// 13/01/2004: start working on it
// 17/02/2004: rewritten to be used with auxiliary variables
// 13/03/2004: added possibility to have a mean of random intercept depended on a total mixture mean
// Remark: * birthDeath changes the overall mixture mean so also log-likelihood of random effects and log-likelihood of observations
//           must be taken into account in Hastings ratio
//         * when computing new value of logLikelihood, old values of all parameters (kP, wM, muM, etc. may be supplied since component pertinence does not change,
//           only a value of Eb0 must be changed)
//             
//
#include "bayessurvreg.h"

using namespace std;

// ********** birthDeath **********
//
// PARAMETERS:
//
// acceptedP ............... +1 if the move will be accepted, +0 if the move will be rejected
// birthP .................. +1 if the birth move was proposed, +0 if the death move was proposed
// kP ...................... on INPUT: current number of mixture components
//                           on OUTPUT: number of mixture components after the birth-death move
// wM ...................... on INPUT: current mixture weights
//                           on OUTPUT: mixture weights after the birth-death move
// muM ..................... on INPUT: current mixture means
//                           on OUTPUT: mixture means after the birth-death move
// invsigma2M .............. on INPUT: current mixture inverse-variances
//                           on OUTPUT: mixture inverse-variances after the birth-death move
// mixMomentM .............. on INPUT: current mixture moments
//                           on OUTPUT: mixture moments after the birth-death move
// rM ...................... on INPUT: current allocations of the observations
//                           on OUTPUT: allocations of the observations after the birth-death move
// invrM ................... on INPUT: current indeces of the observations belonging to the mixture components
//                           on OUTPUT: indeces of the observations belonging to the mixture components 
//                                      after the birth-death move
// mixtureNM ............... on INPUT: current numbers of observations at each mixture component
//                           on OUTPUT: numbers of observations at each mixture component after the birth-death move
// propkP .................. working space to store proposed values of k
// uM ...................... canonical proposal vector of length 3
//                           (pointer to an appropriate part of full uM)
//                           INPUT value from uM[0], uM[1], uM[2] used if the birth move is proposed
//                           OUTPUT value returns in uM[-3], uM[-2], uM[-1] computed canonical proposal vector 
//                                        if the death move is proposed and accepted
// (*logdu) ................ function to compute a log-density of a canonical proposal vector
// (*transu) ............... function used to transform a canonical proposal vector u -> v
// (*invtransu) ............ inverse to that transformation (v -> u
// (*logJtransu) ........... log(Jacobian) of transformation v(u)
// kmaxP ................... maximal number of mixture components
// piBirthM ................ probabilities of the birth move (given k) ((*kmaxP + 1) x 1)
// logpiBirthM ............. log-probabilities of the birth move (given k) ((*kmaxP + 1) x 1)
// logpiDeathM ............. log-probabilities of the death move (given k) ((*kmaxP + 1) x 1)
// deltaP .................. hyperparameter
// xiP ..................... hyperparameter
// invkappaP ............... (transformed) hyperparameter
// sqrtkappaP .............. (transformed) hyperparameter
// halfl2pikappaP .......... (transformed) hyperparameter
// zetaP ................... hyperparameter
// etaP .................... hyperparameter
// lgammazetaP ............. (transformed) hyperparameter
// llambdaP ................ (transformed) hyperparameter
// priorParmu .............. parameters of the prior distribution of the canonical proposal vector u
// transParmu .............. parameters for transformation of a canonical proposal vector u
// priorForkP .............. indicator of the prior for k (0 = Poisson, 1 = uniform)
// Eb0dependMix ............ 0/1 whether Eb0 depend on a mixture or not
// nP ...................... number of observations
//
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
           double (*logdtrans) (double))
{
  if (*priorForkP == Fixed) return;

  int errorType = Mixture;
  bool RichardsonGreen = true;       // Is the proposal for w, mu and invsigma2 same as their prior? I.e. that one
                                     //   used in Richardson and Green (1997)?

  int i, j;
  double unif;
  int accepted = 1;
  int birth = 0;
  double vM[3];             // transformed proposal vector, here it is directly w, mu, invsigma2
  double propuM[3];         // computed proposed canonical proposal vector (used in a combine move)
  double propmixMean, mixEX2;
  double proploglik, proprandomloglik;

  // Determine whether birth or death move will be proposed
  // ========================================================
  unif = runif(0, 1);
  if (unif < piBirthM[*kP]){
    birth = 1;
    *birthP += 1;
  }

  //  cout << (birth ? "BIRTH" : "DEATH") << endl;

  int* emptyCompM = new int[*kP];
  for (j = 0; j < *kP; j++) emptyCompM[j] = -1;
  int numE = numEmpty(emptyCompM, kP, mixtureNM);    // number of empty mixture components and their indeces

  int jdeath;    // index of the component which is to be removed (death move)
                 // OR index of the newly born component (birth move)

  double logduJtransu;               // log of the density of the proposal vector
  double logPostRatioJacobian;
  double logpiDeathBirth;            // log(pi_{k+1}(death) / pi_k(birth)) OR log(pi_{k-1}(birth) / pi_k(death))
  double loglikRatio = 0.0;
  double Paccept;                    // acceptance probability

  if (RichardsonGreen) transParmu[5] = *etaP;

  // Birth move
  // ==========
  if (birth){
    *propkP = *kP + 1;
    if (RichardsonGreen) transParmu[1] = *kP;

    logpiDeathBirth = logpiDeathM[*propkP] - logpiBirthM[*kP] - log(double(numE + 1));

    (*transu)(vM, uM, transParmu);

    logduJtransu = logdtransBirthDeath(vM, uM, priorParmu, transParmu, logdu, logJtransu, &RichardsonGreen);

    logPostRatioJacobian = logPostRatioJacobianBirthDeath(kP, vM, nP,
                                                          deltaP, xiP, invkappaP, halfl2pikappaP, zetaP, etaP, 
                                                          lgammazetaP, llambdaP, priorForkP, &RichardsonGreen);

    if (*Eb0dependMix && (*randomIntP)){
      propmixMean = (1 - vM[0])*mixMomentM[0] + vM[0]*vM[1];
      logLikelihood(&proploglik, *proploglikobs, nP, regresResM, YsM, kP, rM, wM, muM, invsigma2M, &propmixMean, &errorType, randomIntP, logdtrans);
      randomLogLikelihood(&proprandomloglik, *proprandomllcl, &ZERO_INT, nclusterP, nclusterP,
                          bM, betaM, Dcm, &propmixMean, indbinXA);
      loglikRatio = proploglik - (*loglikelhood);
      loglikRatio += proprandomloglik - (*randomloglik);
    }

    Paccept = exp(logPostRatioJacobian + logpiDeathBirth + logduJtransu + loglikRatio);
  }    // end of the birth move

  // Death move
  // ==========
  else{
    accepted = (numE != 0);           // death move is rejected if there are no empty components
    if (!accepted){
      delete [] emptyCompM;
      return;
    }

    *propkP = *kP - 1;
    if (RichardsonGreen) transParmu[1] = *propkP;

    logpiDeathBirth = logpiBirthM[*propkP] - logpiDeathM[*kP] + log(double(numE));

    proposeDeath(jdeath, vM, numE, emptyCompM, wM, muM, invsigma2M);
    
    (*invtransu)(propuM, vM, transParmu);

    logduJtransu = logdtransBirthDeath(vM, propuM, priorParmu, transParmu, logdu, logJtransu, &RichardsonGreen);

    logPostRatioJacobian = logPostRatioJacobianBirthDeath(propkP, vM, nP, 
                                                          deltaP, xiP, invkappaP, halfl2pikappaP, zetaP, etaP, 
                                                          lgammazetaP, llambdaP, priorForkP, &RichardsonGreen);

    if (*Eb0dependMix && (*randomIntP)){
      propmixMean = (mixMomentM[0] - vM[0]*vM[1])/(1 - vM[0]);
      logLikelihood(&proploglik, *proploglikobs, nP, regresResM, YsM, kP, rM, wM, muM, invsigma2M, &propmixMean, &errorType, randomIntP, logdtrans);
      randomLogLikelihood(&proprandomloglik, *proprandomllcl, &ZERO_INT, nclusterP, nclusterP,
                          bM, betaM, Dcm, &propmixMean, indbinXA);
      loglikRatio = (*loglikelhood) - proploglik;
      loglikRatio += (*randomloglik) - proprandomloglik;      
    }

    Paccept = exp(-logPostRatioJacobian + logpiDeathBirth - logduJtransu - loglikRatio);
  }    // end of the death move


  // Accept or reject the proposal
  // =============================
  if (Paccept < 1.0){
    unif = runif(0, 1);
    if (unif > Paccept){
      delete [] emptyCompM;
      return;
    }  
  }
  
  // Proposal was accepted
  // =====================
  *acceptedP += 1;
  moveParamsBirthDeath(jdeath, wM, muM, invsigma2M, rM, invrM, mixtureNM, propkP, vM, &birth);
  *kP = *propkP;
  if (!birth){
    for (i = 0; i < 3; i++) uM[i - 3] = propuM[i];
  }
  mixEX2 = mixMomentM[1]*mixMomentM[1] + mixMomentM[0]*mixMomentM[0];   // original E(epsilon^2) 
  if (*Eb0dependMix && (*randomIntP)){
    mixMomentM[0] = propmixMean;
    *loglikelhood = proploglik;
    *randomloglik = proprandomloglik;
    changePointers(loglikobs, proploglikobs);
    changePointers(randomllcl, proprandomllcl);
  }
  else{
    if (birth) mixMomentM[0] = (1 - vM[0])*mixMomentM[0] + vM[0]*vM[1];
    else       mixMomentM[0] = (mixMomentM[0] - vM[0]*vM[1])/(1 - vM[0]);
  }
  if (birth) mixEX2 = (1 - vM[0])*mixEX2 + vM[0]*(vM[1]*vM[1] + (1/vM[2]));
  else       mixEX2 = (mixEX2 - vM[0]*(vM[1]*vM[1] + (1/vM[2])))/(1 - vM[0]);
  mixMomentM[1] = mixEX2 - mixMomentM[0]*mixMomentM[0];
  if (mixMomentM[1] <= 0.0) mixMomentM[1] = 0.0;
  else                      mixMomentM[1] = sqrt(mixMomentM[1]);
 
  delete [] emptyCompM;
  return;
}   // end of function birthDeath


// ********** proposeDeath **********
//
// Propose which component is to be deleted within the death move
//
// INPUT PARAMETERS:
//
// numE .......... number of empty components
//                 (this number must be at least 1!!!)
// emptyCompM .... indeces of empty components on the first numE places
// wM ............ weights of the mixture components
// muM ........... means of the mixture components
// invsigma2M .... inverse-variances of the mixture components
//
// OUTPUT PARAMETERS:
//
// jdeath .................. index of the component to be killed
// vM[0] = propw ........... weight of the component to be killed
// vM[1] = propmu .......... mean of the component to be killed
// vM[2] = propinvsigma2 ... inv-variance of the component to be killed
//
void
proposeDeath(int& jdeath,        double* vM,
             const int& numE,    const int* emptyCompM,
             const double* wM,   const double* muM,        const double* invsigma2M)
{
  // Choose randomly the index of the component that is to be killed
  discreteUniformSampler(&jdeath, &numE, &ONE_INT, &ZERO_INT);
  jdeath = emptyCompM[jdeath];

  // Give weight, mean and inv-variance of the component to be killed
  vM[0] = wM[jdeath];
  vM[1] = muM[jdeath];
  vM[2] = invsigma2M[jdeath];

  return; 
}    // end of function proposeDeath



// ********** logdtransBirthDeath **********
//
// Compute minus log-density of the proposal or weight, mean and variance 
//   of the component to be removed or of a new component.
// Possibly do not add terms that cancels with log(posterior ratio) when Richardson & Green proposal is used.
//
// The result with (+) corresponds to the birth move, so that it must be multiplied by (-1) to be used
//   with the death move.
//
// RETURN: 
//    If RichardsonGreenP == true, log-proposal density of vM[0] = propw only (the rest cancels with posterior ratio)
//    else                       , log-proposal density evaluated in vM = (propw, propmu, propinvsigma2)
//
// INPUT PARAMETERS:
//
// vM ................. transformed canonical parameter (here: weight, mean, inv-variance
//                      of a component to be killed or born)
// uM ................. original canonical parameter
//
// priorParmu ......... prior parameters to generate u
// transParmu ......... parameters used to transform u -> v
//
double
logdtransBirthDeath(const double* vM,          const double* uM,  
                    const double* priorParmu,  const double* transParmu,          
                    double (*logdu) (const double*, const double*),
                    double (*logJtransu) (const double*, const double*, const double*),
                    const bool* RichardsonGreenP)
{
  double logDens;
  if (*RichardsonGreenP){
    logDens = lbeta(1, transParmu[1]) - (transParmu[1] - 1) * log(1 - vM[0]);     // - log-density of w
  }
  else{
    logDens = (-1) * (*logdu)(uM, priorParmu);
    logDens += (*logJtransu)(uM, vM, transParmu);
  }
 
  return logDens;
}


// ********** logPostRatioJacobianBirthDeath **********
//
// Compute a logarithm of the ratio of posterior that appear in the acceptance probability
//   of the birth-death move.
//
// Include also log(Jacobian) since it cancels with one term in the ratio.
//   (this Jacobian is equal to shortk * log(1 - propw))
//
// The result corresponds to the birth move, so that it must be multiplied by (-1) to be used
//   with the death move.
//
// INPUT PARAMETERS:
//
// shortkP ................. number of components in the shorter state, i.e.
//                           with birth move, *shortkP = # components before the birth move was proposed
//                           with death move, *shortkP = # components after the death move
//
// vM[0] = propw ........... proposed weight for newly born component (birth move)
//                           weight of the component to be removed (death move)
// vM[1] = propmu .......... proposed mean for newly born component (birth move)
//                           mean of the component to be removed (death move)
// vM[2] = propinvsigma2 ... proposed inv-variance for newly born component (birth move)
//                           inv-variance of the component to be removed (death move)
//
// nP ...................... total number of observations
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
// RichardsonGreenP ........ do we use the proposal for w, mu, invsigma2 as suggested by Richardson, Green (1997)?
//
double
logPostRatioJacobianBirthDeath(const int* shortkP,       const double* vM,
                               const int* nP,
                               const double* deltaP,
                               const double* xiP,        const double* invkappaP,  const double* halfl2pikappaP,
                               const double* zetaP,      const double* etaP,       const double* lgammazetaP,
                               const double* llambdaP,   const int* priorForkP,    const bool* RichardsonGreenP)
{
  double log1w = log(1 - vM[0]);
  double logPostRatio = 0.0;

    // log(ratio of priors of allocations r_{i,l})
  logPostRatio += (*nP) * log1w;

    // log(ratio of priors of weights) + log(Jacobian)
  logPostRatio += (*deltaP - 1) * log(vM[0]) + (*shortkP)*(*deltaP) * log1w - lbeta(*deltaP, (*shortkP) * (*deltaP));

  if (!(*RichardsonGreenP)){        // the following terms cancel in the case of Richardson & Green proposal
    // log(ratio of priors of means)
    //  (without adding log(*shortkP + 1)!!!)
    logPostRatio += -(*halfl2pikappaP) - 0.5*(*invkappaP) * (vM[1] - (*xiP)) * (vM[1] - (*xiP));

    // log(ratio of priors of variances)
    if (*etaP <= 0.0) return -FLT_MAX;
    logPostRatio += (*zetaP) * log(*etaP) - (*lgammazetaP) 
                    + (*zetaP + 1) * log(vM[2]) - (*etaP) * vM[2];
  }
                              
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
}   // end of the function logPostRatioBirthDeath


// ********** moveParamsBirthDeath **********
// 
// Move appropriately mixture parameters and pertinence indicators after the accepted birth/death move,
//  recalculate the weights
//
// INPUT PARAMETERS:
//
// jdeath ............... index of the component to be deleted (for death move)
// wM ................... mixture weights before the birth-death move was proposed
// muM .................. mixture means before the birth-death move was proposed
// invsigma2M ........... mixture inverse-variances before the birth-death move was proposed 
// rM ................... mixture pertinence before the birth-death move was proposed
// invrM ................ inverse-mixture pertinence before the birth-death move was proposed
// mixtureNM ............ numbers of observations at each mixture component before the birth-death move was proposed
// propkP ............... number of mixture components after the birth-death move
// v[0] = propw ......... proposed weight 
// v[1] = propmu ........ proposed mean (for birth move)
// v[2] = propinvsigma2 . proposed inverse-variance (for birth move)
// birthP ............... indicator of birth/death
//
// OUTPUT PARAMETERS:
// 
// jdeath ......... index of the newly born component (for birth move)
// wM ............. mixture weights after the birth-death move was performed
// muM ............ mixture means after the birth-death move was performed
// invsigma2M ..... mixture inverse-variances after the birth-death move was performed 
// rM ............. mixture pertinence after the birth-death move was performed
// invrM .......... inverse-mixture pertinence after the birth-death move was performed
// mixtureNM ...... numbers of observations at each mixture component afterthe birth-death move was performed
//
void
moveParamsBirthDeath(int& jdeath,                   
                     double* wM,          double* muM,        double* invsigma2M,
                     int* rM,             List<int>* invrM,   int* mixtureNM,
                     const int* propkP,   const double* vM,   const int* birthP)
{
  int j, obs;

  // Birth move - add the new component to the system and shift mixture parameters with higher indeces
  if (*birthP){
    jdeath = 0;
    while ((jdeath < *propkP - 1) && (muM[jdeath] < vM[1])) jdeath++;

    for (j = *propkP - 1; j > jdeath; j--){
      wM[j] = wM[j - 1] * (1 - vM[0]);
      muM[j] = muM[j - 1];
      invsigma2M[j] = invsigma2M[j - 1];
      mixtureNM[j] = mixtureNM[j - 1];
      invrM[j] = invrM[j - 1];
      for (obs = 0; obs < (invrM[j]).length(); obs++) rM[(invrM[j])[obs]] = j;        
    }    
    wM[jdeath] = vM[0];
    muM[jdeath] = vM[1];
    invsigma2M[jdeath] = vM[2];
    mixtureNM[jdeath] = 0;
    invrM[jdeath] = List<int>();
    for (j = jdeath - 1; j >= 0; j--)
      wM[j] = wM[j] * (1 - vM[0]);
  }

  // Death move - delete one component and shift the rest
  else{   
    for (j = 0; j < jdeath; j++)
      wM[j] = wM[j] / (1 - vM[0]);
    for (j = jdeath; j < *propkP; j++){
      wM[j] = wM[j + 1] / (1 - vM[0]);
      muM[j] = muM[j + 1];
      invsigma2M[j] = invsigma2M[j + 1];
      mixtureNM[j] = mixtureNM[j + 1];
      invrM[j] = invrM[j + 1];
      for (obs = 0; obs < (invrM[j]).length(); obs++) rM[(invrM[j])[obs]] = j;        
    }    
    wM[*propkP] = 0.0;
    muM[*propkP] = 0.0;
    invsigma2M[*propkP] = 0.0;
    mixtureNM[*propkP] = 0;
    invrM[j] = List<int>();
  }

  return;
}    // end of the function moveParamsBirthDeath


// ********** numEmpty **********
//
// Compute the number of empty components
//   and give the vector of their indeces
//
int
numEmpty(int* emptyCompM, const int *kP,  const int* mixtureNM)
{
  int numE = 0;
  for (int j = 0; j < *kP; j++)
    if (!(mixtureNM[j])){
      emptyCompM[numE] = j;
      numE++;
    }

  return numE;
}   // end of the function numEmpty
