// Functions to find a maximum of a function using a Newton-Raphson algorithm
//  and solve an equation using the Newton-Raphson algorithm
//
// 04/11/2004: 'newton_raphson'
//
#include "newton_raphson.h"

#ifdef __cplusplus
extern "C" {
#endif

// =================================================================================
// ***** newton_raphson *****
// =================================================================================
//
// x ............. on INPUT: starting point for the Newton-Raphson
//                 on OUTPUT: point where the maximum is reached
// gx ............ on INPUT : value of g at the starting point
//                 on OUTPUT: value of g at the maximum
// dgx............ on INPUT : value of g' at the starting point
//                 on OUTPUT: value of g' at the maximum
// ddgx........... on INPUT : value of -g'' at the starting point
//                 on OUTPUT: value of -g'' at the maximum
// parmD ......... additional parameters to evaluate g
// parmI ......... additional parameters to evaluate g
// eval2 ......... routine to compute the function we want to maximize and its first and second derivatives
//                 const double* ..... point x where to evalulate the function g
//                 double* ........... g(x)
//                 double* ........... g'(x)
//                 double* ........... -g''(x)
//                 const double* ..... additional parameters to evaluate g
//                 const int* ........ additional parameters to evaluate g
//                 const int& ........ integer indicating what should be computed
//                                   0 = all g(x), g'(x), g''(x)
//                                   1 = only g(x)
//                                   2 = only g'(x) and g''(x)
// iter .......... on OUTPUT: number of iterations needed
// maxiter ....... maximum number of iterations
// max_stephalf... maximum number of step-halving steps
// toler ......... relative toleration to detect convergence
// zero .......... small number to detect singular second derivative
// err ........... error flag
//                 1 = maximum number of step-halving steps was performed without increase of the objective function
//                 2 = maximum number of iterations was used without convergence
//                 3 = bad initials
//                 4 = NaN appeared in g'(x) or g''(x)
//
void
newton_raphson(double* x,            double* gx,          double* dgx,  double* ddgx,  
               const double* parmD,  const int* parmI,
	       void (*eval2)(const double*,  double*, double*, double*, const double*, const int*, const int&),
               int* iter,            const int* maxiter,  const int* max_stephalf,  
               const double* toler,  const double* zero,  int* err)
{
  *err = 0;
  static double newgx, newx, NRstep, relat_diff;
  static int halfstep;

  if (!R_finite(*gx) || !R_finite(*dgx) || !R_finite(*ddgx)){
    *err = 3;
    return;
  }

  for (*iter = 0; *iter < *maxiter; (*iter)++){
    if (fabs(*ddgx) <= *zero) *ddgx = (*zero);
    NRstep = (*dgx)/(*ddgx);
    newx = (*x) + NRstep;
    for (halfstep = 0; halfstep < *max_stephalf; halfstep++){
      eval2(&newx, &newgx, dgx, ddgx, parmD, parmI, 1);
      relat_diff = fabs(1 - ((*gx)/newgx));
      if (newgx >= *gx || relat_diff <= *toler) break;
      newx = 0.5*((*x) + newx);
    }
    if (halfstep == *max_stephalf){
      *err = 1;
      return;
    }
    *x = newx;
    *gx = newgx;
    eval2(x, gx, dgx, ddgx, parmD, parmI, 2);
    if (!R_finite(*dgx) || !R_finite(*ddgx)){
      *err = 4;
      return;
    }    
    if (relat_diff <= *toler) break;
  }
  
  if (*iter == *maxiter) *err = 2;
  return;
}


// =======================================================================================================
// ***** solver_newton_raphson:  Find a solution of an equation g(x) = b using a Newton-Raphson algorithm
// =======================================================================================================
// 
// ASSUMPTIONS: * starting point is in correct region (i.e. only one solution is found)
//              * there exists a solution
//         typical usage: g(x) is concave and unimodal with no limits imposed on the range of the function 'g'
//
// x ............. on INPUT: starting point for the Newton-Raphson
//                 on OUTPUT: point where the equation is solved
// gx ............ on INPUT : value of g at the starting point
//                 on OUTPUT: value of g at the solution (this should be close to b)
// dgx............ on INPUT : value of g' at the starting point
//                 on OUTPUT: value of g' at the solution
// b ............. right hand side
// parmD ......... additional parameters to evaluate g
// parmI ......... additional parameters to evaluate g
// eval2 ......... routine to compute the function g(x)and its first derivative
//                 const double* ..... point x where to evalulate the function g
//                 double* ........... g(x)
//                 double* ........... g'(x)
//                 double* ........... -g''(x) (nowhere needed by this routine)
//                 const double* ..... additional parameters to evaluate g
//                 const int* ........ additional parameters to evaluate g
//                 const int& ........ integer indicating what should be computed
//                                   0 = all g(x), g'(x), g''(x)
//                                   1 = only g(x)
//                                   2 = only g'(x) and g''(x)
//                                   3 = only g(x) and g'(x)
// iter .......... on OUTPUT: number of iterations needed
// maxiter ....... maximum number of iterations
// toler ......... relative toleration to detect convergence
// zero .......... small number to detect singular first derivative
// err ........... error flag
//                 (1) = maximum number of step-halving steps was performed without increase of the objective function (not applicable here)
//                 2   = maximum number of iterations was used without convergence
//                 3   = bad initials
//                 4   = NaN appeared in g(x) or g'(x)
//
void
solver_newton_raphson(double* x,            double* gx,          double* dgx,  const double* b,
                      const double* parmD,  const int* parmI,
	              void (*eval2)(const double*,  double*, double*, double*, const double*, const int*, const int&),
                      int* iter,            const int* maxiter,
                      const double* toler,  const double* zero,  int* err)
{
  *err = 0;
  static double ddgx, NRstep, _diff;

  if (!R_finite(*gx) || !R_finite(*dgx) || !R_finite(*b)){
    *err = 3;
    return;
  }

  _diff = *b - (*gx);
  for (*iter = 0; *iter < *maxiter; (*iter)++){
    if (fabs(*dgx) <= *zero) *dgx = (*zero);
    NRstep = _diff/(*dgx);
    (*x) += NRstep;
    eval2(x, gx, dgx, &ddgx, parmD, parmI, 3);
    if (!R_finite(*gx) || !R_finite(*dgx)){
      *err = 4;
      return;
    }    
    _diff = *b - (*gx);
    if (fabs(_diff/(*b)) <= *toler) break;
  }
  
  if (*iter == *maxiter) *err = 2;
  return;
}


#ifdef __cplusplus
}                          /* end of extern "C" */
#endif
