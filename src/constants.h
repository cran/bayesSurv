#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

enum typeError {Mixture, Spline, PolyaTree, WhoKnows};
enum typeUpdate {RandomWalk, AdaptiveM, Gibbs, Fixed};
enum priorMatrix {InvWishart, SDUniform};
enum priorRandomEff {Normal_, Gspline_};
enum typePrior_k {Uniform, Poisson, Fixed_k};

const int AKINT_MAX = 999999;
const int ONE_INT = 1;
const int ZERO_INT = 0;

const double LOG_SQRT_2PI = 0.918938533204672741780329736406;	     // log(sqrt(2*pi)) 
const double NORM_ZERO = 1e-16;                                      // qnorm(1 - 1e-16) is still non infty in R
const double ZERO = 1e-50;
const double invFLT_MAX = 1e-50;
const double LOG_2 = 0.6931472;
const double SCALE_ZERO = 1e-20;
const double invSQRT_TWO_PI = 0.39894228040143270286;



#endif
