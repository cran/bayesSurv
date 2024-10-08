/* adaptive rejection metropolis sampling */

/* *********************************************************************** */

// Original arms.c written by P. Wild and W. R. Gilks
//
// This is a very slight modification of the original program made by 
// Arnost Komarek
// 08/12/2006
// 19/06/2012  exit() replaced by error() to avoid RCMD check complains
// 16/09/2024  error() replaced by Rf_error() to avoid RCMD check complains
//             Calloc() and Free() replaced by R_Calloc() and R_Free() to avoid RCMD check complains
//
//
// It is supposed to be used in R
//  (e.g. R uniform random numbers generator is specified here)
//
// ORIGINAL DOCUMENTATION:
// =======================
//
// ARMS - Adaptive Rejection Metropolis Sampling
// _____________________________________________
//
// Arguments of the C function
// ___________________________
//
// Note: ARMS can be called either through the arms function, or though the
// arms_simple function, which has a simplified parameter list. A non-zero
// return code indicates an error.
//
// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
//
// int arms (double *xinit, int ninit, double *xl, double *xr,
//           double (*myfunc)(double x, void *mydata), void *mydata,
//           double *convex, int npoint, int dometrop, double *xprev,
//           double *xsamp, int nsamp, double *qcent, double *xcent, int ncent,
//           int *neval);
//
// int arms_simple (int ninit, double *xl, double *xr,
// 	            double (*myfunc)(double x, void *mydata), void *mydata,
//                  int dometrop, double *xprev, double *xsamp);
//
// -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
//
// xinit        : Pointer to first element of an array storing
//                starting values for x in ascending order. SEE
//                IMPORTANT NOTES ABOVE ON CHOICE OF STARTING VALUES
//                (INITIAL CONSTRUCTION POINTS).
// 
// ninit        : Number of starting values supplied. Usually four
//                will suffice.
// 
// xl           : Pointer to left bound. Bounds must be supplied even
//                if the support is truly unbounded. This should
//                present no problem if the bounds are set far enough out.
// 
// xr           : Pointer to right bound. Comments as for left bound.
// 
// myfunc       : A function to evaluate the log-density (user-supplied).
//                This function should have two arguments, the first a 
//                double holding the point x at which the log-density is
//                to be evaluated; the second a pointer to a user-defined
//                struct which must contain all the information needed
//                to calculate the log-density at x. The function must
//                return the value of the log-density (determined up to
//                an arbitrary additive constant) at x.
// 
// mydata       : A pointer to the struct holding the data required by myfunc.
//  
// convex       : A pointer to a double. The double should be non-negative.
//                We suggest setting the double to 1.0. If ARMS tends to keep
//                sampling the same value, try a larger setting.
//                Larger settings will make a larger envelope
//                in regions of non-log-concavity.
//  
// npoint       : Maximum number of envelope points to be allowed. Setting 
//                npoint=50 should be ample.
//
// dometrop     : Whether metropolis step is required. Set dometrop=0 if you
//                know the log-density is concave (see above); if you know it 
//                isn't, or you're not sure, set dometrop=1.
//
// xprev        : Pointer to a double holding the previous value from the Markov
//                chain (see above).
// 
// xsamp        : Pointer to the first element of an array which is to hold
//                sampled values from the target density.
//
// nsamp        : The number of sampled values to be obtained (normally
//                nsamp=1 in applications of Gibbs sampling).
//  
// qcent        : Pointer to start of an array holding requested percentages
//                for envelope centiles. This may be useful if you want to
//                generate some good initial values for the next updating
//                of the current parameter (BUT SEE CAUTION ABOVE).
//  
// xcent        : Pointer to start of an array to hold the requested centiles.
//  
// ncent        : Number of centiles requested. Set ncent=0 if you don't want 
//                centiles to be calculated.
//
// neval        : Pointer to an int which will hold, on exit, the number of
//               log-density evaluations that were performed.
//
// IMPLEMENTATION NOTE
// ___________________ 
//
// The program contains the statement:
// 
// #define RAND_MAX 2147483647      /* For Sun4 */
// 
// The constant RAND_MAX is used in function u_random, which returns a standard
// uniform random variate. You might need to comment out this statement, or
// reset it to be the largest integer that can be returned by the
// library function rand().
//    AK (08/12/2006):  Function u_random() has been removed and replaced by R unif_rand() in the code
//
//
// References
// __________
// 
// Gilks, W. R. (1992) Derivative-free adaptive rejection sampling
//   for Gibbs sampling. Bayesian Statistics 4, (eds. Bernardo, J.,
//   Berger, J., Dawid, A. P., and Smith, A. F. M.) Oxford 
//   University Press.
//
// Gilks, W. R., Best, N. G. and Tan, K. K. C. (1995) Adaptive
//   rejection Metropolis sampling. Applied Statistics, 44, 455-472.
//
// Gilks, W. R. and Wild, P. (1992) Adaptive rejection sampling
//   for Gibbs sampling. Applied Statistics 41, pp 337-348.
//
// Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H.
//   and Teller, E. (1953) Equations of state calculations by fast
//   computing machines. J. Chem. Phys., 21, 1087-1092.
//
// Ripley, B. (1987) Stochastic Simulation. New York, Wiley. 
//
//
#include "arms.h"


int arms_simple (int ninit, double *xl, double *xr,
                 double (*myfunc)(double x, void *mydata), void *mydata,
                 int dometrop, double *xprev, double *xsamp)

/* adaptive rejection metropolis sampling - simplified argument list */
/* ninit        : number of starting values to be used */
/* *xl          : left bound */
/* *xr          : right bound */
/* *myfunc      : function to evaluate log-density */
/* *mydata      : data required by *myfunc */
/* dometrop     : whether metropolis step is required */
/* *xprev       : current value from markov chain */
/* *xsamp       : to store sampled value */

{
  double *xinit = R_Calloc(ninit, double); 
  double convex = 1.0;
  double qcent, xcent;
  int err, i, npoint=100, nsamp=1, ncent=0, neval; 
 
  /* set up starting values */
  for(i=0; i<ninit; i++){
    xinit[i] = *xl + (i + 1.0) * (*xr - *xl)/(ninit + 1.0);
  }

  err = arms(xinit,ninit,xl,xr,myfunc,mydata,&convex,npoint,dometrop,xprev,xsamp,
             nsamp,&qcent,&xcent,ncent,&neval);

  R_Free(xinit);
  return err;
}

/* *********************************************************************** */

int arms (double *xinit, int ninit, double *xl, double *xr, 
	 double (*myfunc)(double x, void *mydata), void *mydata,
         double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
         int nsamp, double *qcent, double *xcent,
         int ncent, int *neval)

/* to perform derivative-free adaptive rejection sampling with metropolis step */
/* *xinit       : starting values for x in ascending order */
/* ninit        : number of starting values supplied */
/* *xl          : left bound */
/* *xr          : right bound */
/* *myfunc      : function to evaluate log-density */
/* *mydata      : data required by *myfunc */
/* *convex      : adjustment for convexity */
/* npoint       : maximum number of envelope points */
/* dometrop     : whether metropolis step is required */
/* *xprev       : previous value from markov chain */
/* *xsamp       : to store sampled values */
/* nsamp        : number of sampled values to be obtained */
/* *qcent       : percentages for envelope centiles */
/* *xcent       : to store requested centiles */
/* ncent        : number of centiles requested */
/* *neval       : on exit, the number of function evaluations performed */

{

  ENVELOPE *env;      /* rejection envelope */
  POINT pwork;        /* a working point, not yet incorporated in envelope */
  int msamp=0;        /* the number of x-values currently sampled */
  FUNBAG lpdf;        /* to hold density function and its data */
  METROPOLIS *metrop; /* to hold bits for metropolis step */
  int i,err;

  /* check requested envelope centiles */
  for(i=0; i<ncent; i++){
    if((qcent[i] < 0.0) || (qcent[i] > 100.0)){
      /* percentage requesting centile is out of range */
      return 1005;
    }
  }

  /* incorporate density function and its data into FUNBAG lpdf */
  lpdf.mydata = mydata;
  lpdf.myfunc = myfunc;

  /* set up space required for envelope */
  env = (ENVELOPE *)malloc(sizeof(ENVELOPE));
  if(env == NULL){
    /* insufficient space */
    return 1006;
  }

  /* start setting up metropolis struct */
  metrop = (METROPOLIS *)malloc(sizeof(METROPOLIS));
  if(metrop == NULL){
    /* insufficient space */
    return 1006;
  }
  metrop->on = dometrop;

  /* set up initial envelope */
  err = initial(xinit,ninit,*xl,*xr,npoint,&lpdf,env,convex,
        neval,metrop);
  if(err)return err;

  /* finish setting up metropolis struct (can only do this after */
  /* setting up env) */
  if(metrop->on){
    if((*xprev < *xl) || (*xprev > *xr)){
      /* previous markov chain iterate out of range */
      return 1007;
    }
    metrop->xprev = *xprev;
    metrop->yprev = perfunc(&lpdf,env,*xprev);
  }

  /* now do adaptive rejection */
  do {
    /* sample a new point */
    sample (env,&pwork);

    /* perform rejection (and perhaps metropolis) tests */
    i = test(env,&pwork,&lpdf,metrop);
    if(i == 1){
      /* point accepted */
      xsamp[msamp++] = pwork.x;
    } else if (i != 0) {
      /* envelope error - violation without metropolis */
      return 2000;
    }  
  } while (msamp < nsamp);

  /* nsamp points now sampled */
  /* calculate requested envelope centiles */
  for (i=0; i<ncent; i++){
    invert(qcent[i]/100.0,env,&pwork);
    xcent[i] = pwork.x;
  }

  /* free space */
  free(env->p);
  free(env);
  free(metrop);

  return 0;
}

/* *********************************************************************** */

int initial (double *xinit, int ninit, double xl, double xr, int npoint,
	     FUNBAG *lpdf, ENVELOPE *env, double *convex, int *neval,
             METROPOLIS *metrop)

/* to set up initial envelope */
/* xinit        : initial x-values */
/* ninit        : number of initial x-values */
/* xl,xr        : lower and upper x-bounds */
/* npoint       : maximum number of POINTs allowed in envelope */
/* *lpdf        : to evaluate log density */
/* *env         : rejection envelope attributes */
/* *convex      : adjustment for convexity */
/* *neval       : current number of function evaluations */
/* *metrop      : for metropolis step */

{
  int i,j,k,mpoint;
  POINT *q;

  if(ninit<3){
    /* too few initial points */
    return 1001;
  }

  mpoint = 2*ninit + 1;
  if(npoint < mpoint){
    /* too many initial points */
    return 1002;
  }

  if((xinit[0] <= xl) || (xinit[ninit-1] >= xr)){
    /* initial points do not satisfy bounds */
    return 1003;
  }

  for(i=1; i<ninit; i++){
    if(xinit[i] <= xinit[i-1]){
      /* data not ordered */
      return 1004;
    }
  }

  if(*convex < 0.0){
    /* negative convexity parameter */
    return 1008;
  }

  /* copy convexity address to env */
  env->convex = convex;

  /* copy address for current number of function evaluations */
  env->neval = neval;
  /* initialise current number of function evaluations */
  *(env->neval) = 0;

  /* set up space for envelope POINTs */
  env->npoint = npoint;
  env->p = (POINT *)malloc(npoint*sizeof(POINT));
  if(env->p == NULL){
    /* insufficient space */
    return 1006;
  }

  /* set up envelope POINTs */
  q = env->p;
  /* left bound */
  q->x = xl;
  q->f = 0;
  q->pl = NULL;
  q->pr = q+1;
  for(j=1, k=0; j<mpoint-1; j++){
    q++;
    if(j%2){
      /* point on log density */
      q->x = xinit[k++];
      q->y = perfunc(lpdf,env,q->x);
      q->f = 1;
    } else {
      /* intersection point */
      q->f = 0;
    }
    q->pl = q-1;
    q->pr = q+1;
  }
  /* right bound */
  q++;
  q->x = xr;
  q->f = 0;
  q->pl = q-1;
  q->pr = NULL;

  /* calculate intersection points */
  q = env->p;
  for (j=0; j<mpoint; j=j+2, q=q+2){
    if(meet(q,env,metrop)){
      /* envelope violation without metropolis */
      return 2000;
    }
  }

  /* exponentiate and integrate envelope */
  cumulate(env);

  /* note number of POINTs currently in envelope */
  env->cpoint = mpoint;

  return 0;
}

/* *********************************************************************** */

void sample(ENVELOPE *env, POINT *p)

/* To sample from piecewise exponential envelope */
/* *env    : envelope attributes */
/* *p      : a working POINT to hold the sampled value */

{
  double prob;

  /* sample a uniform */
  prob = unif_rand();
  /* get x-value correponding to a cumulative probability prob */
  invert(prob,env,p);

  return;
}

/* *********************************************************************** */

void invert(double prob, ENVELOPE *env, POINT *p)

/* to obtain a point corresponding to a qiven cumulative probability */
/* prob    : cumulative probability under envelope */
/* *env    : envelope attributes */
/* *p      : a working POINT to hold the sampled value */

{
  double u,xl,xr,yl,yr,eyl,eyr,prop;
  POINT *q;

  /* find rightmost point in envelope */
  q = env->p;
  while(q->pr != NULL)q = q->pr;

  /* find exponential piece containing point implied by prob */
  u = prob * q->cum;
  while(q->pl->cum > u)q = q->pl;

  /* piece found: set left and right POINTs of p, etc. */
  p->pl = q->pl;
  p->pr = q;
  p->f = 0;
  p->cum = u;

  /* calculate proportion of way through integral within this piece */
  prop = (u - q->pl->cum) / (q->cum - q->pl->cum);

  /* get the required x-value */
  if (q->pl->x == q->x){
    /* interval is of zero length */
    p->x = q->x;
    p->y = q->y;
    p->ey = q->ey;

    xl = p->x;      /* These lines added by AK on 20130530 to avoid warning: variable 'xl' is used uninitialized whenever 'if' condition is true */
    xr = p->x;      /* which relates to use of xl and xr in if ((p->x < xl) || (p->x > xr)) below. This definition of xl, xr should never lead   */
                    /* to error("arms error 1\n").                                                                                               */
  } else {
    xl = q->pl->x;
    xr = q->x;
    yl = q->pl->y;
    yr = q->y;
    eyl = q->pl->ey;
    eyr = q->ey;
    if(fabs(yr - yl) < YEPS){
      /* linear approximation was used in integration in function cumulate */
      if(fabs(eyr - eyl) > EYEPS*fabs(eyr + eyl)){
	p->x = xl + ((xr - xl)/(eyr - eyl))
	       * (-eyl + sqrt((1. - prop)*eyl*eyl + prop*eyr*eyr));
      } else {
	p->x = xl + (xr - xl)*prop;
      }
      p->ey = ((p->x - xl)/(xr - xl)) * (eyr - eyl) + eyl;
      p->y = logshift(p->ey, env->ymax);
    } else {
      /* piece was integrated exactly in function cumulate */
      p->x = xl + ((xr - xl)/(yr - yl))
	      * (-yl + logshift(((1.-prop)*eyl + prop*eyr), env->ymax));
      p->y = ((p->x - xl)/(xr - xl)) * (yr - yl) + yl;
      p->ey = expshift(p->y, env->ymax);
    }
  }

  /* guard against imprecision yielding point outside interval */
  //if ((p->x < xl) || (p->x > xr))exit(1);
  if ((p->x < xl) || (p->x > xr)) Rf_error("arms error 1\n");

  return;
}

/* *********************************************************************** */

int test(ENVELOPE *env, POINT *p, FUNBAG *lpdf, METROPOLIS *metrop)

/* to perform rejection, squeezing, and metropolis tests */
/* *env          : envelope attributes */
/* *p            : point to be tested */
/* *lpdf         : to evaluate log-density */
/* *metrop       : data required for metropolis step */

{
  double u,y,ysqueez,ynew,yold,znew,zold,w;
  POINT *ql,*qr;
  
  /* for rejection test */
  u = unif_rand() * p->ey;
  y = logshift(u,env->ymax);

  if(!(metrop->on) && (p->pl->pl != NULL) && (p->pr->pr != NULL)){
    /* perform squeezing test */
    if(p->pl->f){
      ql = p->pl;
    } else {
      ql = p->pl->pl;
    }
    if(p->pr->f){
      qr = p->pr;
    } else {
      qr = p->pr->pr;
    }
    ysqueez = (qr->y * (p->x - ql->x) + ql->y * (qr->x - p->x))
               /(qr->x - ql->x);
    if(y <= ysqueez){
      /* accept point at squeezing step */
      return 1;
    }
  }

  /* evaluate log density at point to be tested */
  ynew = perfunc(lpdf,env,p->x);
  
  /* perform rejection test */
  if(!(metrop->on) || ((metrop->on) && (y >= ynew))){
    /* update envelope */
    p->y = ynew;
    p->ey = expshift(p->y,env->ymax);
    p->f = 1;
    if(update(env,p,lpdf,metrop)){
      /* envelope violation without metropolis */
      return -1;
    }
    /* perform rejection test */
    if(y >= ynew){
      /* reject point at rejection step */
      return 0;
    } else {
      /* accept point at rejection step */
      return 1;
    }
  }

  /* continue with metropolis step */
  yold = metrop->yprev;
  /* find envelope piece containing metrop->xprev */
  ql = env->p;
  while(ql->pl != NULL)ql = ql->pl;
  while(ql->pr->x < metrop->xprev)ql = ql->pr;
  qr = ql->pr;
  /* calculate height of envelope at metrop->xprev */
  w = (metrop->xprev - ql->x)/(qr->x - ql->x);
  zold = ql->y + w*(qr->y - ql->y);
  znew = p->y;
  if(yold < zold)zold = yold;
  if(ynew < znew)znew = ynew;
  w = ynew-znew-yold+zold;
  if(w > 0.0)w = 0.0;

  if(w > -YCEIL){
    w = exp(w);
  } else {
    w = 0.0;
  }
  u = unif_rand();
  if(u > w){
    /* metropolis says dont move, so replace current point with previous */
    /* markov chain iterate */
    p->x = metrop->xprev;
    p->y = metrop->yprev;
    p->ey = expshift(p->y,env->ymax);
    p->f = 1;
    p->pl = ql;
    p->pr = qr;
  } else {
    /* trial point accepted by metropolis, so update previous markov */
    /* chain iterate */
    metrop->xprev = p->x;
    metrop->yprev = ynew;
  }
  return 1;
}

/* *********************************************************************** */

int update(ENVELOPE *env, POINT *p, FUNBAG *lpdf, METROPOLIS *metrop)

/* to update envelope to incorporate new point on log density*/
/* *env          : envelope attributes */
/* *p            : point to be incorporated */
/* *lpdf         : to evaluate log-density */
/* *metrop       : for metropolis step */

{
  POINT *m,*ql,*qr,*q;

  if(!(p->f) || (env->cpoint > env->npoint - 2)){
    /* y-value has not been evaluated or no room for further points */
    /* ignore this point */
    return 0;
  }

  /* copy working POINT p to a new POINT q */
  q = env->p + env->cpoint++;
  q->x = p->x;
  q->y = p->y;
  q->f = 1;

  /* allocate an unused POINT for a new intersection */
  m = env->p + env->cpoint++;
  m->f = 0;
  if((p->pl->f) && !(p->pr->f)){
    /* left end of piece is on log density; right end is not */
    /* set up new intersection in interval between p->pl and p */
    m->pl = p->pl;
    m->pr = q;
    q->pl = m;
    q->pr = p->pr;
    m->pl->pr = m;
    q->pr->pl = q;
  } else if (!(p->pl->f) && (p->pr->f)){
    /* left end of interval is not on log density; right end is */
    /* set up new intersection in interval between p and p->pr */
    m->pr = p->pr;
    m->pl = q;
    q->pr = m;
    q->pl = p->pl;
    m->pr->pl = m;
    q->pl->pr = q;
  } else {
    /* this should be impossible */
    //exit(10);
    Rf_error("arms error 10");
  }

  /* now adjust position of q within interval if too close to an endpoint */
  if(q->pl->pl != NULL){
    ql = q->pl->pl;
  } else {
    ql = q->pl;
  }
  if(q->pr->pr != NULL){
    qr = q->pr->pr;
  } else {
    qr = q->pr;
  }
  if (q->x < (1. - XEPS) * ql->x + XEPS * qr->x){
    /* q too close to left end of interval */
    q->x = (1. - XEPS) * ql->x + XEPS * qr->x;
    q->y = perfunc(lpdf,env,q->x);
  } else if (q->x > XEPS * ql->x + (1. - XEPS) * qr->x){
    /* q too close to right end of interval */
    q->x = XEPS * ql->x + (1. - XEPS) * qr->x;
    q->y = perfunc(lpdf,env,q->x);
  }

  /* revise intersection points */
  if(meet(q->pl,env,metrop)){
    /* envelope violation without metropolis */
    return 1;
  }
  if(meet(q->pr,env,metrop)){
    /* envelope violation without metropolis */
    return 1;
  }
  if(q->pl->pl != NULL){
    if(meet(q->pl->pl->pl,env,metrop)){
      /* envelope violation without metropolis */
      return 1;
    }
  }
  if(q->pr->pr != NULL){
    if(meet(q->pr->pr->pr,env,metrop)){
      /* envelope violation without metropolis */
      return 1;
    }
  }

  /* exponentiate and integrate new envelope */
  cumulate(env);

  return 0;
}

/* *********************************************************************** */

void cumulate(ENVELOPE *env)

/* to exponentiate and integrate envelope */
/* *env     : envelope attributes */

{
  POINT *q,*qlmost;

  qlmost = env->p;
  /* find left end of envelope */
  while(qlmost->pl != NULL)qlmost = qlmost->pl;

  /* find maximum y-value: search envelope */
  env->ymax = qlmost->y;
  for(q = qlmost->pr; q != NULL; q = q->pr){
    if(q->y > env->ymax)env->ymax = q->y;
  }

  /* exponentiate envelope */
  for(q = qlmost; q != NULL; q = q->pr){
    q->ey = expshift(q->y,env->ymax);
  }

  /* integrate exponentiated envelope */
  qlmost->cum = 0.;
  for(q = qlmost->pr; q != NULL; q = q->pr){
    q->cum = q->pl->cum + area(q);
  }

  return;
}

/* *********************************************************************** */

int meet (POINT *q, ENVELOPE *env, METROPOLIS *metrop)
/* To find where two chords intersect */
/* q         : to store point of intersection */
/* *env      : envelope attributes */
/* *metrop   : for metropolis step */

{
  double gl,gr,grl,dl,dr;
  int il,ir,irl;

  if(q->f){
    /* this is not an intersection point */
    //exit(30);
    Rf_error("arms error 30");
  }

  /* calculate coordinates of point of intersection */
  if ((q->pl != NULL) && (q->pl->pl->pl != NULL)){
    /* chord gradient can be calculated at left end of interval */
    gl = (q->pl->y - q->pl->pl->pl->y)/(q->pl->x - q->pl->pl->pl->x);
    il = 1;
  } else {
    /* no chord gradient on left */
    il = 0;
  }
  if ((q->pr != NULL) && (q->pr->pr->pr != NULL)){
    /* chord gradient can be calculated at right end of interval */
    gr = (q->pr->y - q->pr->pr->pr->y)/(q->pr->x - q->pr->pr->pr->x);
    ir = 1;
  } else {
    /* no chord gradient on right */
    ir = 0;
  }
  if ((q->pl != NULL) && (q->pr != NULL)){
    /* chord gradient can be calculated across interval */
    grl = (q->pr->y - q->pl->y)/(q->pr->x - q->pl->x);
    irl = 1;
  } else {
    irl = 0;
  }

  if(irl && il && (gl<grl)){
    /* convexity on left exceeds current threshold */
    if(!(metrop->on)){
      /* envelope violation without metropolis */
      return 1;
    }
    /* adjust left gradient */
    gl = gl + (1.0 + *(env->convex)) * (grl - gl);
  }

  if(irl && ir && (gr>grl)){
    /* convexity on right exceeds current threshold */
    if(!(metrop->on)){
      /* envelope violation without metropolis */
      return 1;
    }
    /* adjust right gradient */
    gr = gr + (1.0 + *(env->convex)) * (grl - gr);
  }

  if(il && irl){
    dr = (gl - grl) * (q->pr->x - q->pl->x);
    if(dr < YEPS){
      /* adjust dr to avoid numerical problems */
      dr = YEPS;
    }
  }

  if(ir && irl){
    dl = (grl - gr) * (q->pr->x - q->pl->x);
    if(dl < YEPS){
      /* adjust dl to avoid numerical problems */
      dl = YEPS;
    }
  }

  if(il && ir && irl){
    /* gradients on both sides */
    q->x = (dl * q->pr->x + dr * q->pl->x)/(dl + dr);
    q->y = (dl * q->pr->y + dr * q->pl->y + dl * dr)/(dl + dr);
  } else if (il && irl){
    /* gradient only on left side, but not right hand bound */
    q->x = q->pr->x;
    q->y = q->pr->y + dr;
  } else if (ir && irl){
    /* gradient only on right side, but not left hand bound */
    q->x = q->pl->x;
    q->y = q->pl->y + dl;
  } else if (il){
    /* right hand bound */
    q->y = q->pl->y + gl * (q->x - q->pl->x);
  } else if (ir){
    /* left hand bound */
    q->y = q->pr->y - gr * (q->pr->x - q->x);
  } else {
    /* gradient on neither side - should be impossible */
    //exit(31);
    Rf_error("arms error 31");
  }
  if(((q->pl != NULL) && (q->x < q->pl->x)) ||
     ((q->pr != NULL) && (q->x > q->pr->x))){
    /* intersection point outside interval (through imprecision) */
    //exit(32);
    Rf_error("arms error 32");
  }
  /* successful exit : intersection has been calculated */
  return 0;
}

/* *********************************************************************** */

double area(POINT *q)

/* To integrate piece of exponentiated envelope to left of POINT q */

{
  double a;

  if(q->pl == NULL){
    /* this is leftmost point in envelope */
    //exit(1);
    Rf_error("arms error 1");
  } else if(q->pl->x == q->x){
    /* interval is zero length */
    a = 0.;
  } else if (fabs(q->y - q->pl->y) < YEPS){
    /* integrate straight line piece */
    a = 0.5*(q->ey + q->pl->ey)*(q->x - q->pl->x);
  } else {
    /* integrate exponential piece */
    a = ((q->ey - q->pl->ey)/(q->y - q->pl->y))*(q->x - q->pl->x);
  }
  return a;
}

/* *********************************************************************** */

double expshift(double y, double y0)

/* to exponentiate shifted y without underflow */
{
  if(y - y0 > -2.0 * YCEIL){
    return exp(y - y0 + YCEIL);
  } else {
    return 0.0;
  }
}

/* *********************************************************************** */

double logshift(double y, double y0)

/* inverse of function expshift */
{
  return (log(y) + y0 - YCEIL);
}

/* *********************************************************************** */

double perfunc(FUNBAG *lpdf, ENVELOPE *env, double x)

/* to evaluate log density and increment count of evaluations */

/* *lpdf   : structure containing pointers to log-density function and data */
/* *env    : envelope attributes */
/* x       : point at which to evaluate log density */

{
  double y;

  /* evaluate density function */
  y = (lpdf->myfunc)(x,lpdf->mydata);

  /* increment count of function evaluations */
  (*(env->neval))++;

  return y;
}

/* *********************************************************************** */

void display(FILE *f, ENVELOPE *env)

/* to display envelope - for debugging only */
/* all fprintf commands commented by AK on 05/10/2022 to avoid CRAN warnings:    */
/* arms.cpp:917:45: warning: format specifies type 'unsigned int' but the argument has type 'POINT *' (aka 'point *') [-Wformat] */
/* arms.cpp:927:56: warning: format specifies type 'unsigned int' but the argument has type 'POINT *' (aka 'point *') [-Wformat] */
/* arms.cpp:927:58: warning: format specifies type 'unsigned int' but the argument has type 'struct point *' [-Wformat]          */
/* arms.cpp:927:64: warning: format specifies type 'unsigned int' but the argument has type 'struct point *' [-Wformat]          */
    
{
  POINT *q;

  /* print envelope attributes */
  //fprintf(f,"========================================================\n");
  //fprintf(f,"envelope attributes:\n");
  //fprintf(f,"points in use = %d, points available = %d\n",
  //        env->cpoint,env->npoint);
  //fprintf(f,"function evaluations = %d\n",*(env->neval));
  //fprintf(f,"ymax = %f, p = %x\n",env->ymax,env->p);
  //fprintf(f,"convexity adjustment = %f\n",*(env->convex));
  //fprintf(f,"--------------------------------------------------------\n");

  /* find leftmost POINT */
  q = env->p;
  while(q->pl != NULL)q = q->pl;

  /* now print each POINT from left to right */
  for(q = env->p; q != NULL; q = q->pr){
    //fprintf(f,"point at %x, left at %x, right at %x\n",q,q->pl,q->pr);
    //fprintf(f,"x = %f, y = %f, ey = %f, cum = %f, f = %d\n",
    //        q->x,q->y,q->ey,q->cum,q->f);
  }
  //fprintf(f,"========================================================\n");

  return;
}

