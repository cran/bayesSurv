// Function to update possibly censored data distributed as a G-spline
//
// 22/10/2004: 'update_Data_GS'
// 10/12/2004: 2nd specification of the G-spline implemented in 'update_Data_GS'
// 12/01/2005: 'update_Data_GS_doubly'
//             'update_Data_GS_regres'
// 20/01/2005: 'update_Data_GS_doubly' rewritten
// 06/06/2013: 'update_Data_GS_regres_misclass' added
//
#include "update_Data_GS.h"

// ****** update_Data_GS ***********************
//
// Version without regression
// =============================================
//
// YsM[nP x gg->dim()] .... (augmented) log-event times
// rM[nP] ................. component labels taking values 0, 1, ..., gg->total_length()-1
//
void
update_Data_GS(double* YsM,
               const double*  y_left,  
               const double*  y_right,   
               const int*     status,
               const int*     rM,         
               const Gspline* gg,
               const int*     nP,         
               const int*     n_censored)
{
  if (!(*n_censored)) return;

  int obs, j;
  double mu_jk = 0;
  double PhiL = 0;
  double PhiU = 0;  
  double u = 0;
  double PhiInv = 0;
  double stres = 0;

  double invsigma[_max_dim];
  double invscale[_max_dim];
  for (j = 0; j < gg->dim(); j++){
    invsigma[j] = 1/gg->sigma(j);
    invscale[j] = 1/gg->scale(j);
  }

  double* y_obs = YsM;
  const double* y1 = y_left;
  const double* y2 = y_right;
  const int* stat = status;
  const int* rp = rM;
  for (obs = 0; obs < *nP; obs++){
    for (j = 0; j < gg->dim(); j++){

      switch (*stat){
      case 1:   /* exactly observed */
        break;

      case 0:   /* right censored */
        mu_jk = gg->mu_component(j, *rp);
        stres = (*y1 - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
        PhiL = pnorm(stres, 0, 1, 1, 0);
        if (PhiL >= 1 - NORM_ZERO){        // censored time irrealistic large (out of the prob. scale)
          *y_obs = *y1;
        }
        else{
          if (PhiL <= NORM_ZERO){         // censoring time equal to "zero", generate an exact time from N(mean, variance), 
                                          //   i.e. from the full  not-truncated distribution
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            *y_obs = gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
          }
          else{
            u = runif(0, 1) * (1 - PhiL) + PhiL;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            if (PhiInv == R_PosInf){    // u was equal to 1, additional check added 16/12/2004
              *y_obs = *y1;
            }
            else{
              *y_obs = gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv;
            }
          }  
        }            
        break;

      case 2:   /* left censored */
        mu_jk = gg->mu_component(j, *rp);
        stres = (*y1 - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
        PhiU = pnorm(stres, 0, 1, 1, 0);
        if (PhiU <= NORM_ZERO){           // left censoring time irrealistic low (equal to "zero")
          *y_obs = *y1;
        }
        else{
          if (PhiU >= 1 - NORM_ZERO){    // left censoring time equal to "infty", generate an exact time from N(mean, variance), 
                                           //   i.e. from the full  not-truncated distribution
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            *y_obs = gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
          }
          else{
            u = runif(0, 1) * PhiU;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            if (PhiInv == R_NegInf){  // u was equal to 0,  additional check added 16/12/2004
              *y_obs = *y1;
            }
            else{
              *y_obs = gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
            }
          }
        }
        break;

      case 3:   /* interval censored */
        mu_jk = gg->mu_component(j, *rp);
        stres = (*y1 - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
        PhiL = pnorm(stres, 0, 1, 1, 0);
        stres = (*y2 - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
        PhiU = pnorm(stres, 0, 1, 1, 0);
        PhiInv = PhiU - PhiL;
        if (PhiInv <= NORM_ZERO){       // too narrow interval, or the interval out of the probability scale
                                        //   (both limits in "zero" probability region)
                                        //   generate something inbetween
          u = runif(0, 1);
          *y_obs = *y1 + u*((*y2) - (*y1)); 
        }
        else{
          if (PhiInv >= 1 - NORM_ZERO){ // too large interval, practically (-infty, +infty), generate an exact time from N(mean, variance)
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            *y_obs = gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
          }
          else{
            u = runif(0, 1) * PhiInv + PhiL;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            if (!R_finite(PhiInv)){    // u was either zero or one,  additional check added 16/12/2004
              u = runif(0, 1);
              *y_obs = *y1 + u*((*y2) - (*y1)); 
            }
            else{
              *y_obs = gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
            }
          }  
        }      
        break;
      }  /** end of switch (status) **/

      /*** This section just performs additional checks to prevent simulations with NaN's ***/
      if (!R_finite(*y_obs)){
        int condit;
        REprintf("\nY[%d,%d]=%e, r[%d,%d]=%d,  status[%d,%d]=%d,  stres=%e", 
		 obs, j, *y_obs, obs, j, *rp, obs, j, *stat, stres);
        REprintf(";  mean=%e", mu_jk); 
        REprintf(";  invvar=%e", gg->invsigma2(j)); 
        REprintf("\nu=%3.20e,  PhiL=%3.20e,  PhiU=%3.20e,  PhiInv=%3.20e", u, PhiL, PhiU, PhiInv);
        REprintf("NORM_ZERO=%3.20e,  1-NORM_ZERO=%3.20e", NORM_ZERO, 1-NORM_ZERO);
        switch (*stat){
        case 0:
          condit = 1*(PhiL >= 1 - NORM_ZERO);
          REprintf("\nPhiL >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiL <= NORM_ZERO);
          REprintf("\nPhiL <= NORM_ZERO: %d", condit);
          break;
        case 2:
          condit = 1*(PhiU >= 1 - NORM_ZERO);
          REprintf("\nPhiU >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiU <= NORM_ZERO);
          REprintf("\nPhiU <= NORM_ZERO: %d", condit);
          break;
        case 3:
          condit = 1*(PhiU-PhiL >= 1 - NORM_ZERO);
          REprintf("\nPhiU-PhiL >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiU-PhiL <= NORM_ZERO);
          REprintf("\nPhiU-PhiL <= NORM_ZERO: %d", condit);
          break;
        }        
        REprintf("\n");
        throw returnR("Trap in update_Data_GS: NaN generated.", 1);
      }

      y_obs++;
      y1++;
      y2++;
      stat++;
    }
    rp++;
  }
  
  return;
}    /*** end of function update_Data_GS ***/


// ****** update_Data_GS_regres ***********************
//
// Version with possible regression
// ================================
//
// YsM[nP x gg->dim()] ........... on INPUT:  current vector of (imputed) log(event times)
//                                 on OUTPUT: updated vector of (augmented) log(event times)
// regresResM[nP x gg->dim()] .... on INPUT:  current vector of regression residuals (y - x'beta - z'b))
//                                 on OUTPUT: updated vector of regression residuals
//
// rM[nP] ........................ component labels taking values 0, 1, ..., gg->total_length()-1
//
void
update_Data_GS_regres(double* YsM,           
                      double* regresResM,
                      const double*  y_left,  
                      const double*  y_right,   
                      const int*     status,
                      const int*     rM,         
                      const Gspline* gg,       
                      const int* nP)
{
  int obs, j;
  double mu_jk = 0;
  double PhiL = 0;
  double PhiU = 0;  
  double u = 0;
  double PhiInv = 0;
  double stres = 0;

  double invsigma[_max_dim];
  double invscale[_max_dim];
  for (j = 0; j < gg->dim(); j++){
    invsigma[j] = 1/gg->sigma(j);
    invscale[j] = 1/gg->scale(j);
  }

  //Rprintf("\nG-spline dim: %d\n", gg->dim());
  //Rprintf("mu[0, 0]  = %g\n", gg->mu_component(0, 0));
  //Rprintf("sigma[0]  = %g\n", gg->sigma(0));
  //Rprintf("intcpt[0] = %g\n", gg->intcpt(0));
  //Rprintf("scale[0]  = %g\n", gg->scale(0));      

  double* y_obs = YsM;
  double* regRes = regresResM;
  const double* y1 = y_left;
  const double* y2 = y_right;
  const int* stat = status;
  const int* rp = rM;
  for (obs = 0; obs < *nP; obs++){
    for (j = 0; j < gg->dim(); j++){

      switch (*stat){
      case 1:   /* exactly observed */
        break;

      case 0:   /* right censored */
        mu_jk = gg->mu_component(j, *rp);
        *regRes -= *y_obs;
        stres = (*y1 + (*regRes) - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
        PhiL = pnorm(stres, 0, 1, 1, 0);
        if (PhiL >= 1 - NORM_ZERO){        // censored time irrealistic large (out of the prob. scale)
          *y_obs = *y1;
        }
        else{
          if (PhiL <= NORM_ZERO){         // censoring time equal to "zero", generate an exact time from N(mean, variance), 
                                          //   i.e. from the full  not-truncated distribution
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            *y_obs = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
          }
          else{
            u = runif(0, 1) * (1 - PhiL) + PhiL;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            if (PhiInv == R_PosInf){    // u was equal to 1, additional check added 16/12/2004
              *y_obs = *y1;
            }
            else{
              *y_obs = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv;
            }
          }  
        }
        *regRes += (*y_obs);
        break;

      case 2:   /* left censored */
        mu_jk = gg->mu_component(j, *rp);
        *regRes -= *y_obs;
        stres = (*y1 + (*regRes) - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
        PhiU = pnorm(stres, 0, 1, 1, 0);
        if (PhiU <= NORM_ZERO){           // left censoring time irrealistic low (equal to "zero")
          *y_obs = *y1;
        }
        else{
          if (PhiU >= 1 - NORM_ZERO){    // left censoring time equal to "infty", generate an exact time from N(mean, variance), 
                                           //   i.e. from the full  not-truncated distribution
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            *y_obs = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
          }
          else{
            u = runif(0, 1) * PhiU;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            if (PhiInv == R_NegInf){  // u was equal to 0,  additional check added 16/12/2004
              *y_obs = *y1;
            }
            else{
              *y_obs = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
            }
          }
        }
        *regRes += *y_obs;
        break;

      case 3:   /* interval censored */
        mu_jk = gg->mu_component(j, *rp);
        *regRes -= *y_obs;
        stres = (*y1 + (*regRes) - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
        PhiL = pnorm(stres, 0, 1, 1, 0);
        stres = (*y2 + (*regRes) - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
        PhiU = pnorm(stres, 0, 1, 1, 0);
        PhiInv = PhiU - PhiL;
        if (PhiInv <= NORM_ZERO){       // too narrow interval, or the interval out of the probability scale
                                        //   (both limits in "zero" probability region)
                                        //   generate something inbetween
          u = runif(0, 1);
          *y_obs = *y1 + u*((*y2) - (*y1)); 
        }
        else{
          if (PhiInv >= 1 - NORM_ZERO){ // too large interval, practically (-infty, +infty), generate an exact time from N(mean, variance)
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            *y_obs = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
          }
          else{
            u = runif(0, 1) * PhiInv + PhiL;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            if (!R_finite(PhiInv)){    // u was either zero or one,  additional check added 16/12/2004
              u = runif(0, 1);
              *y_obs = *y1 + u*((*y2) - (*y1)); 
            }
            else{
              *y_obs = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
            }
          }  
        }      
        *regRes += *y_obs;
        break;
      }  /** end of switch (status) **/

      /*** This section just performs additional checks to prevent simulations with NaN's ***/
      if (!R_finite(*y_obs) || !R_finite(*regRes)){
        int condit;
        REprintf("\nY[%d,%d]=%e,  regRes[%d,%d]=%e,  r[%d,%d]=%d,  status[%d,%d]=%d,  stres=%e", 
		 obs, j, *y_obs, obs, j, *regRes, obs, j, *rp, obs, j, *stat, stres);
        REprintf(";  mean=%e", mu_jk); 
        REprintf(";  invvar=%e", gg->invsigma2(j)); 
        REprintf("\nu=%3.20e,  PhiL=%3.20e,  PhiU=%3.20e,  PhiInv=%3.20e", u, PhiL, PhiU, PhiInv);
        REprintf("NORM_ZERO=%3.20e,  1-NORM_ZERO=%3.20e", NORM_ZERO, 1-NORM_ZERO);
        switch (*stat){
        case 0:
          condit = 1*(PhiL >= 1 - NORM_ZERO);
          REprintf("\nPhiL >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiL <= NORM_ZERO);
          REprintf("\nPhiL <= NORM_ZERO: %d", condit);
          break;
        case 2:
          condit = 1*(PhiU >= 1 - NORM_ZERO);
          REprintf("\nPhiU >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiU <= NORM_ZERO);
          REprintf("\nPhiU <= NORM_ZERO: %d", condit);
          break;
        case 3:
          condit = 1*(PhiU-PhiL >= 1 - NORM_ZERO);
          REprintf("\nPhiU-PhiL >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiU-PhiL <= NORM_ZERO);
          REprintf("\nPhiU-PhiL <= NORM_ZERO: %d", condit);
          break;
        }        
        REprintf("\n");
        throw returnR("Trap in update_Data_GS_regres: NaN generated.", 1);
      }

      y_obs++;
      regRes++;
      y1++;
      y2++;
      stat++;
    }
    rp++;
  }
  
  return;
}    /*** end of function update_Data_GS_regres ***/


// ****** update_Data_GS_doubly ***********************
// *** Update of the event-time in the case of doubly censored data
//
// Yevent[nP x gg->dim()] ........ on INPUT:  current vector of (imputed) log(event times)
//                                 on OUTPUT: updated vector of (augmented) log(event times)
//                                 i.e. augmented log(T2 - T1), where T1 = onset time, T2 = event time (on a study scale)
// regresResM[nP x gg->dim()] .... on INPUT:  current vector of regression residuals (y - x'beta - z'b))
//                                 on OUTPUT: updated vector of regression residuals
// Yonset[nP x gg->dim()] .... log-onset times 
//                             i.e. log(T1)
// t_left[nP x gg->dim()] .... 
// t_right[nP x gg->dim()].... observed event times (on a study scale)
// status[nP x gg->dim()] .... censoring status for event
// rM[nP] .................... component labels taking values 0, 1, ..., gg->total_length()-1
// gg ........................ G-spline defining the distribution of the log-time-to-event (log(T2 - T1))
// nP ........................ number of observational vectors
// n_censored ................ number of censored event times
//
void
update_Data_GS_doubly(double* Yevent,        
                      double* regresResM,
                      const double*  Yonset, 
  	              const double*  t_left,  
                      const double*  t_right,    
                      const int*     status,
                      const int*     rM,         
                      const Gspline* gg,       
                      const int*     nP)
{
  int obs, j;
  double t_onset, yL, yU, help;
  double mu_jk = 0; 
  double PhiL = 0;
  double PhiU = 0;  
  double u = 0;
  double PhiInv = 0;
  double stres = 0;

  double invsigma[_max_dim];
  double invscale[_max_dim];
  for (j = 0; j < gg->dim(); j++){
    invsigma[j] = 1/gg->sigma(j);
    invscale[j] = 1/gg->scale(j);
  }

  double* y_event = Yevent;
  double* regRes = regresResM;
  const double* y_onset = Yonset;
  const double* t1 = t_left;
  const double* t2 = t_right;
  const int* stat = status;
  const int* rp = rM;
  for (obs = 0; obs < *nP; obs++){
    for (j = 0; j < gg->dim(); j++){

      t_onset = (*y_onset > -_emax ? exp(*y_onset) : 0.0);
      if (!R_finite(t_onset)) throw returnR("Trap: t_onset equal to NaN in 'update_Data_GS_doubly'", 1);      

      *regRes -= *y_event;     
      switch (*stat){
      case 1:   /* exactly observed, but the onset time might not be observed exactly */
        help = (*t1) - t_onset;
        if (help <= _ZERO_TIME_) *y_event = _LOG_ZERO_TIME_;
        else                     *y_event = log(help);
        break;

      case 0:   /* right censored */
        mu_jk = gg->mu_component(j, *rp);
        help = (*t1) - t_onset;
        if (help <= _ZERO_TIME_){      // time-to-event right censored at 0, generate an exact time from N(mean, variance)
          u = runif(0, 1);
          PhiInv = qnorm(u, 0, 1, 1, 0);
          *y_event = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
        }
        else{
          yL = log(help);
          stres = (yL + (*regRes) - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
          PhiL = pnorm(stres, 0, 1, 1, 0);
          if (PhiL >= 1 - NORM_ZERO){        // censored time irrealistic large (out of the prob. scale)
            *y_event = yL;
          }
          else{
            if (PhiL <= NORM_ZERO){         // censoring time equal to "zero", generate an exact time from N(mean, variance), 
                                            //   i.e. from the full  not-truncated distribution
              u = runif(0, 1);
              PhiInv = qnorm(u, 0, 1, 1, 0);
              *y_event = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
            }  
            else{
              u = runif(0, 1) * (1 - PhiL) + PhiL;
              PhiInv = qnorm(u, 0, 1, 1, 0);
              if (PhiInv == R_PosInf){    // u was equal to 1, additional check added 16/12/2004
                *y_event = yL;
              }
              else{
                *y_event = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv;
              }
            }  
          }
        }
        break;

      case 2:   /* left censored event => onset had to be left censored as well at the same time */
        mu_jk = gg->mu_component(j, *rp);
        help = (*t1) - t_onset;
        if (help <= _ZERO_TIME_) *y_event = _LOG_ZERO_TIME_;        // time-to-event left censored at 0 => time-to-event = 0
        else{
          yL = log(help);
          stres = (yL + (*regRes) - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
          PhiU = pnorm(stres, 0, 1, 1, 0);
          if (PhiU <= NORM_ZERO){           // left censoring time irrealistic low (equal to "zero")
            *y_event = _LOG_ZERO_TIME_;
          }
          else{
            if (PhiU >= 1 - NORM_ZERO){      // left censoring time equal to "infty", generate an exact time from N(mean, variance), 
                                             //   i.e. from the full  not-truncated distribution
              u = runif(0, 1);
              PhiInv = qnorm(u, 0, 1, 1, 0);
              *y_event = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
            }
            else{
              u = runif(0, 1) * PhiU;
              PhiInv = qnorm(u, 0, 1, 1, 0);
              if (PhiInv == R_NegInf){  // u was equal to 0,  additional check added 16/12/2004
                *y_event = yL;
              }
              else{
                *y_event = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
              }
            }
          }
        }
        break;

      case 3:   /* interval censored */
        mu_jk = gg->mu_component(j, *rp);

        help = (*t1) - t_onset;
        if (help <= _ZERO_TIME_){        // time-to-event will be left censored
          help = (*t2) - t_onset;
          if (help <= _ZERO_TIME_){      // too narrow interval located close to zero
            *y_event = _LOG_ZERO_TIME_;
          }
          else{                          // code for left censored observations
            yL = log(help);
            stres = (yL + (*regRes) - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
            PhiU = pnorm(stres, 0, 1, 1, 0);
            if (PhiU <= NORM_ZERO){           // left censoring time irrealistic low (equal to "zero")
              *y_event = _LOG_ZERO_TIME_;
            }
            else{
              if (PhiU >= 1 - NORM_ZERO){      // left censoring time equal to "infty", generate an exact time from N(mean, variance), 
                                               //   i.e. from the full  not-truncated distribution
                u = runif(0, 1);
                PhiInv = qnorm(u, 0, 1, 1, 0);
                *y_event = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
              }
              else{
                u = runif(0, 1) * PhiU;
                PhiInv = qnorm(u, 0, 1, 1, 0);
                if (PhiInv == R_NegInf){  // u was equal to 0,  additional check added 16/12/2004
                  *y_event = yL;
                }
                else{
                  *y_event = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
                }
              }
            }            
          }
        }
        else{
          yL = log(help);

          help = (*t2) - t_onset;
          if (help <= _ZERO_TIME_){      // too narrow interval located close to zero
            *y_event = _LOG_ZERO_TIME_;
          }
          else{
            yU = log(help);

            stres = (yL + (*regRes) - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
            PhiL = pnorm(stres, 0, 1, 1, 0);
            stres = (yU + (*regRes) - gg->intcpt(j) - gg->scale(j)*mu_jk) * invsigma[j]*invscale[j];
            PhiU = pnorm(stres, 0, 1, 1, 0);
            PhiInv = PhiU - PhiL;
            if (PhiInv <= NORM_ZERO){       // too narrow interval, or the interval out of the probability scale
                                            //   (both limits in "zero" probability region)
                                            //   generate something inbetween
              u = runif(0, 1);
              *y_event = yL + u*(yU - yL); 
            }
            else{
              if (PhiInv >= 1 - NORM_ZERO){ // too large interval, practically (-infty, +infty), generate an exact time from N(mean, variance)
                u = runif(0, 1);
                PhiInv = qnorm(u, 0, 1, 1, 0);
                *y_event = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
              }
              else{
                u = runif(0, 1) * PhiInv + PhiL;
                PhiInv = qnorm(u, 0, 1, 1, 0);
                if (!R_finite(PhiInv)){    // u was either zero or one,  additional check added 16/12/2004
                  u = runif(0, 1);
                  *y_event = yL + u*(yU - yL); 
                }
                else{
                  *y_event = -(*regRes) + gg->intcpt(j) + gg->scale(j)*mu_jk + gg->sigma(j)*gg->scale(j)*PhiInv; 
                }
              }  
            }      
          }
        }
        break;
      }  /** end of switch (status) **/
      *regRes += (*y_event);


      /*** This section just performs additional checks to prevent simulations with NaN's ***/
      if (!R_finite(*y_event) || !R_finite(*regRes)){
        int condit;
        REprintf("\nY[%d,%d]=%e,  regRes[%d,%d]=%e,  r[%d,%d]=%d,  status[%d,%d]=%d,  stres=%e", 
		 obs, j, *y_event, obs, j, *regRes, obs, j, *rp, obs, j, *stat, stres);
        REprintf(";  mean=%e", mu_jk); 
        REprintf(";  invvar=%e", gg->invsigma2(j)); 
        REprintf("\nu=%3.20e,  PhiL=%3.20e,  PhiU=%3.20e,  PhiInv=%3.20e", u, PhiL, PhiU, PhiInv);
        REprintf("NORM_ZERO=%3.20e,  1-NORM_ZERO=%3.20e", NORM_ZERO, 1-NORM_ZERO);
        switch (*stat){
        case 0:
          condit = 1*(PhiL >= 1 - NORM_ZERO);
          REprintf("\nPhiL >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiL <= NORM_ZERO);
          REprintf("\nPhiL <= NORM_ZERO: %d", condit);
          break;
        case 2:
          condit = 1*(PhiU >= 1 - NORM_ZERO);
          REprintf("\nPhiU >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiU <= NORM_ZERO);
          REprintf("\nPhiU <= NORM_ZERO: %d", condit);
          break;
        case 3:
          condit = 1*(PhiU-PhiL >= 1 - NORM_ZERO);
          REprintf("\nPhiU-PhiL >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiU-PhiL <= NORM_ZERO);
          REprintf("\nPhiU-PhiL <= NORM_ZERO: %d", condit);
          break;
        }        
        REprintf("\n");
        throw returnR("Trap in update_Data_GS_doubly: NaN generated.", 1);
      }

      y_event++;
      regRes++;
      y_onset++;
      t1++;
      t2++;
      stat++;
    }
    rp++;
  }
  
  return;
}    /*** end of function update_Data_GS_doubly ***/


// ****** update_Data_GS_regres_misclass ***********************
//
// Version with possible regression and misclassification of the event status
// ============================================================================
//
// REMARK to calculation of iPML: iPML calculated here conditions also by component allocations.
//                                This is different to 'iPML_misclass_GJK' function below
//                                which integrates the component allocations out by using the
//                                whole mixture in the expression for iPML.
//
// Created in 201305 by modification of 'update_Data_GS_regres' function.
// -----------------------------------------------------------------------------
//
// This function assumes that gg->dim() = 1.
//
// YsM[nP x gg->dim()] ........... on INPUT:  current vector of (imputed) log(event times)
//                                 on OUTPUT: updated vector of (augmented) log(event times)
// regresResM[nP x gg->dim()] .... on INPUT:  current vector of regression residuals (y - x'beta - z'b))
//                                 on OUTPUT: updated vector of regression residuals
// n00[nExaminer * nFactor] ...... INPUT:  whatsever
//                                 OUTPUT: numbers of (0-0) correctly classified events for each examiner:factor
// n10[nExaminer * nFactor] ...... INPUT:  whatsever
//                                 OUTPUT: numbers of (Classification = 1 | True = 0) incorrectly classified events for each examiner:factor
// n01[nExaminer * nFactor] ...... INPUT:  whatsever
//                                 OUTPUT: numbers of (Classification = 0 | True = 1) incorrectly classified events for each examiner:factor
// n11[nExaminer * nFactor] ...... INPUT:  whatsever
//                                 OUTPUT: numbers of (1-1) correctly classified events for each examiner:factor
//
// iPML[nExaminer * nFactor] ..... INPUT:  whatsever
//                                 OUTPUT: individual contributions needed to calculate the pseudo marginal likelihood and also the deviance
//
// dwork[(1 + max(nvisit)) * 6] .. working array
//
// sens[nExaminer * nFactor]...... sensitivities for each examiner:factor
// spec[nExaminer * nFactor]...... specificities for each examiner:factor
// logvtime[nP * sum(nvisit)] .... logarithms of visit times for each observation
// status[nP * sum(nvisit)] ...... classified event status for each visit
// 
// nvisit[nP] .................... numbers of visits for each observation
// Examiner[nP * sum(nvisit)] .... examiner (0, 1, ..., nExaminer - 1) identification at each visit
// Factor[nP * sum(nvisit)] ...... factor (0, 1, ..., nFactor - 1) identification at each visit
//
// rM[nP] ........................ component labels taking values 0, 1, ..., gg->total_length()-1
//
void
update_Data_GS_regres_misclass(double* YsM,           
                               double* regresResM,
                               int*    n00,
                               int*    n10,
                               int*    n01,
                               int*    n11,
                               double* iPML,
                               double* dwork,
                               const double*  sens,
                               const double*  spec,
                               const double*  logvtime,
                               const int*     status,
                               const int*     nExaminer,
                               const int*     nFactor,
                               const int*     nvisit,
                               const int*     maxnvisit,
                               const int*     Examiner,
                               const int*     Factor,
                               const int*     rM,         
                               const Gspline* gg,       
                               const int*     nP)
{ 
  if (gg->dim() > 1) REprintf("update_Data_GS_regres_misclass: Error, not implemented for gg->dim() > 1.\n");

  /*** Some general variables ***/
  int    obs, m, k, L;
  double mu_i = 0;
  double u    = 0;
  double Phi  = 0;
  double stres_sampled = 0;

  double invsigma_invscale = 1 / (gg->sigma(0) * gg->scale(0));

  /*** Working arrays and related variables ***/
  double *A      = dwork;                             /* A numbers                                                                    */
  double *cumInt = A + (1 + *maxnvisit);              /* cumsum(A * int_{y_{k-1}}^{y_k} f(s)ds), the last is the normalizing constant */
                                                      /* and also the value of iPML                                                   */
  double *cprod_sens = cumInt + (1 + *maxnvisit);     /* cumulative product needed for 'A's based on sensitivities                    */
  double *cprod_spec = cprod_sens + (1 + *maxnvisit); /* cumulative product needed for 'A's based on specificities                    */
  double *stres_cut  = cprod_spec + (1 + *maxnvisit); /* limits of intervals on the scale of standardized residuals                   */ 
  double *Phi_cut    = stres_cut + (1 + *maxnvisit);  /* Phi(stres_cut)                                                               */

  double *A_k;
  double *cumInt_k; 
  double *cprod_sens_k;
  double *cprod_spec_k;
  double *stres_cut_k;
  double *Phi_cut_k;

  /*** Reset classification matrices ***/
  int* n00P = n00;
  int* n10P = n10;
  int* n01P = n01;
  int* n11P = n11;
  for (m = 0; m < *nExaminer * *nFactor; m++){
    *n00P = 0;
    *n10P = 0;
    *n01P = 0;
    *n11P = 0;

    n00P++;
    n10P++;
    n01P++;
    n11P++;    
  }

  /*** Main loop over observations ***/
  double* y_i      = YsM;
  double* regRes_i = regresResM;
  double* iPML_i   = iPML;
  
  const int*    nvisit_i    = nvisit;

  const double* logvtime_i = logvtime;
  const double* logvtime_ik;

  const int* status_i = status;
  const int* status_ik;

  const int* Examiner_i = Examiner;
  const int* Examiner_ik;

  const int* Factor_i = Factor;
  const int* Factor_ik;
  
  const int* r_i = rM;

  for (obs = 0; obs < *nP; obs++){

    mu_i = gg->mu_component(0, *r_i);
    *regRes_i -= *y_i;

    /*** Calculate cumulative products based on specificities needed for 'A' numbers ***/
    cprod_spec_k = cprod_spec;
    *cprod_spec_k = 1.0;            /* k = 0*/
    cprod_spec_k++;

    status_ik   = status_i;
    Examiner_ik = Examiner_i;
    Factor_ik   = Factor_i;
    for (k = 1; k <= *nvisit_i; k++){
      *cprod_spec_k = *(cprod_spec_k - 1) * (*status_ik == 1 ? (1 - spec[*nFactor * *Examiner_ik + *Factor_ik]) : spec[*nFactor * *Examiner_ik + *Factor_ik]);
      cprod_spec_k++;
      status_ik++;
      Examiner_ik++;
      Factor_ik++;
    }
    
    /*** Calculate cumulative products based on sensitivities needed for 'A' numbers ***/
    cprod_sens_k = cprod_sens + *nvisit_i;
    *cprod_sens_k = 1.0;                    /* k = nvisit */
    cprod_sens_k--;

    status_ik--;
    Examiner_ik--;
    Factor_ik--;
    for (k = *nvisit_i - 1; k >= 0; k--){
      *cprod_sens_k = *(cprod_sens_k + 1) * (*status_ik == 1 ? sens[*nFactor * *Examiner_ik + *Factor_ik] : (1 - sens[*nFactor * *Examiner_ik + *Factor_ik]));
      cprod_sens_k--;
      status_ik--;
      Examiner_ik--;
      Factor_ik--;
    }    

    /*** Calculate the 'A' numbers and 'cumInt' for this observation ***/
    A_k          = A;
    cprod_sens_k = cprod_sens;
    cprod_spec_k = cprod_spec;
    cumInt_k     = cumInt;
    stres_cut_k  = stres_cut;
    Phi_cut_k    = Phi_cut;
    logvtime_ik  = logvtime_i;

    /** k = 0: first visit - like left-censored) **/
    *A_k = *cprod_sens_k * *cprod_spec_k;

    *stres_cut_k = (*logvtime_ik + (*regRes_i) - gg->intcpt(0) - gg->scale(0) * mu_i) * invsigma_invscale;
    *Phi_cut_k   = pnorm(*stres_cut_k, 0, 1, 1, 0);

    *cumInt_k = *A_k * *Phi_cut_k;

    A_k++;
    cprod_sens_k++;
    cprod_spec_k++;
    cumInt_k++;
    stres_cut_k++;
    Phi_cut_k++;
    logvtime_ik++;

    /** k = 1, ..., *nvisit_i - 1: like interval-censored **/
    for (k = 1; k < *nvisit_i; k++){      
      *A_k = *cprod_sens_k * *cprod_spec_k;
   
      *stres_cut_k = (*logvtime_ik + (*regRes_i) - gg->intcpt(0) - gg->scale(0) * mu_i) * invsigma_invscale;
      *Phi_cut_k   = pnorm(*stres_cut_k, 0, 1, 1, 0);

      *cumInt_k = *(cumInt_k - 1) + *A_k * (*Phi_cut_k - *(Phi_cut_k - 1));

      A_k++;
      cprod_sens_k++;
      cprod_spec_k++;
      cumInt_k++;
      stres_cut_k++;
      Phi_cut_k++;
      logvtime_ik++;
    }

    /** k = *nvisit_i: like right-censored **/
    *A_k = *cprod_sens_k * *cprod_spec_k;
    *cumInt_k = *(cumInt_k - 1) + *A_k * (1 - *(Phi_cut_k - 1));

    /** Normalizing constant and also iPML **/
    *iPML_i = *cumInt_k;

    /** Debuging section **/
    //if (obs == 5){
    //  Rprintf("alpha <- c("); for (k = 0; k < *nFactor * *nExaminer; k++) Rprintf("%g, ", sens[k]); Rprintf(")\n");
    //  Rprintf("eta <- c(");   for (k = 0; k < *nFactor * *nExaminer; k++) Rprintf("%g, ", spec[k]); Rprintf(")\n");
    //  Rprintf("nvisit = %d\n", *nvisit_i);
    //  Rprintf("   logv <- c("); for (k = 0; k < *nvisit_i; k++) Rprintf("%g, ", logvtime_i[k]); Rprintf(")\n");
    //  Rprintf("   stres <- c("); for (k = 0; k < *nvisit_i; k++) Rprintf("%g, ", stres_cut[k]); Rprintf(")\n");
    //  Rprintf("   Phi  <- c("); for (k = 0; k < *nvisit_i; k++) Rprintf("%g, ", Phi_cut[k]); Rprintf(")\n");
    //  Rprintf("   Y <- c("); for (k = 0; k < *nvisit_i; k++) Rprintf("%d, ", status_i[k]); Rprintf(")\n");
    //  Rprintf("   A <- c("); for (k = 0; k <= *nvisit_i; k++) Rprintf("%g, ", A[k]); Rprintf(")\n");
    //  Rprintf("   cumInt <- c("); for (k = 0; k <= *nvisit_i; k++) Rprintf("%g, ", cumInt[k]); Rprintf(")\n\n");
    //}

    /** Sample a uniform random variable **/
    u = runif(0, 1);

    /** Find out to which piece the 'u' value points out **/
    cumInt_k    = cumInt;
    A_k         = A;
    //stres_cut_k = stres_cut;
    Phi_cut_k   = Phi_cut;
    for (L = 0; L < *nvisit_i; L++){
      if (u <=  *cumInt_k / *iPML_i) break;
      cumInt_k++;
      A_k++;
      //stres_cut_k++;
      Phi_cut_k++;
    }
    /*** Now: L = 0:                  u belongs to piece (-infty, vtime[0]],       A_k = A[0],   stres_cut_k = stres[0]      ***/
    /***      L = 1:                  u belongs to piece (vtime[0], vtime[1]],     A_k = A[1],   stres_cut_k = stres[1]      ***/
    /***      ...                                                                                                            ***/
    /***      L = nvisit - 1 = K - 1: u belongs to piece (vtime[K-2], vtime[K-1]], A_k = A[K-1], stres_cut_k = stres[K-1]    ***/
    /***      L = nvisit = K        : u belongs to piece (vtime[K-1], infty),      A_k = A[K],   stres_cut_k = N.A.          ***/

    /*** Get the sampled value of the standardized residual ***/
    if (L == 0){                     /*** Like LEFT-CENSORED observation     ***/
      Phi = (*iPML_i * u) / *A_k;
    }else{                           /*** L = 1, ..., nvisit_i: Like INTERVAL or RIGHT-CENSORED observation ***/
      Phi = (*iPML_i * u - *(cumInt_k - 1)) / *A_k + *(Phi_cut_k - 1);
    }
    if (Phi <= NORM_ZERO){        // PROBLEM1: stres_sampled = -infty
      stres_sampled = -QNORM_ONE;
    }else{
      if (Phi >= 1 - NORM_ZERO){  // PROBLEM2: stres_sampled = infty
        stres_sampled = QNORM_ONE;
      }else{                      // NO PROBLEMS
        stres_sampled = qnorm(Phi, 0, 1, 1, 0);
      }
    }

    /*** Calculate the sampled value of the log event time and the regression residual ***/
    *y_i = gg->sigma(0) * gg->scale(0) * stres_sampled - *regRes_i + gg->intcpt(0) + gg->scale(0) * mu_i; 
    *regRes_i += *y_i;

    /*** Update the classification matrices                                           ***/
    /*** Shift pointers logvtime_i, status_i, Examiner_i, Factor_i at the same time.  ***/
    for (k = 0; k < *nvisit_i; k++){
      if (*y_i <= *logvtime_i){    /*** True status is 1. ***/
        if (*status_i == 1){          /** Correct (1, 1)   **/
          n11[*nFactor * *Examiner_i + *Factor_i] += 1;
        }else{                        /** Incorrect (0, 1) **/
          n01[*nFactor * *Examiner_i + *Factor_i] += 1;
        }
      }else{                       /*** True status is 0. ***/
        if (*status_i == 1){          /** Incorrect (1, 0)   **/
          n10[*nFactor * *Examiner_i + *Factor_i] += 1;
        }else{                        /** Correct (0, 0) **/
          n00[*nFactor * *Examiner_i + *Factor_i] += 1;
        } 
      }
      logvtime_i++;
      status_i++;
      Examiner_i++;
      Factor_i++;
    }

    /*** Shift remaining pointers ***/
    y_i++;
    regRes_i++;
    r_i++;
    nvisit_i++;
    iPML_i++;
  }
  
  return;
}    /*** end of function update_Data_GS_regres_misclass ***/


// ****** iPML_misclass_GJK ***********************
//
// This is more or less a subest of update_Data_GS_regres_misclass which only calculates iPML.
// Its main purpose is to use it to calculate "D.in.bar" by a simple .C call from R.
//
// Difference in calculation of iPML compared to 'update_Data_GS_regres_misclass' function:
//   * component allocations are integrated out here by using the whole mixture in the definition
//     of iPML
// 
// It assumes only dim = 1 G-spline and all needed G-spline parameters are supplied as simple arrays.
//
// min_etaM[nP] ................ minus linear predictors (including the random effects)
//                               = regresResM - YsM from the 'update_...' function
//
// dwork[] ................. this might be shorter than the same argument of 'update_Data_GS_regres_misclass'
//
extern "C"{

void
iPML_misclass_GJK(double* iPML,                  
                  double* dwork,
                  const double* min_etaM,
                  const double* sens,
                  const double* spec,
                  const double* logvtime,
                  const int*    status,
                  const int*    nExaminer,
                  const int*    nFactor,
                  const int*    nvisit,
                  const int*    maxnvisit,
                  const int*    Examiner,
                  const int*    Factor,
                  const int*    gg_K,
                  const double* gg_gamma,
                  const double* gg_delta,         
                  const double* gg_sigma,        
                  const double* gg_intcpt,        
                  const double* gg_scale,
                  const double* gg_w,        
                  const int*    nP)
{ 
  /*** Some general variables ***/
  int    obs, k, r;
  double mu_r = 0;

  double invsigma_invscale = 1 / (*gg_sigma * *gg_scale);

  /*** Working arrays and related variables ***/
  double *A      = dwork;                             /* A numbers                                                                    */
  double *cprod_sens = A + (1 + *maxnvisit);          /* cumulative product needed for 'A's based on sensitivities                    */
  double *cprod_spec = cprod_sens + (1 + *maxnvisit); /* cumulative product needed for 'A's based on specificities                    */
  //double *next  = cprod_spec + (1 + *maxnvisit);    /* limits of intervals on the scale of standardized residuals                   */ 

  double *A_k;
  double *cprod_sens_k;
  double *cprod_spec_k;
  double stres_cut_k, stres_cut_common_r;
  double Phi_cut_k, Phi_cut_k_min_1;

  /*** Main loop over observations ***/
  const double* min_eta_i  = min_etaM;
  double* iPML_i = iPML;
  
  const int*    nvisit_i    = nvisit;

  const double* logvtime_i = logvtime;
  const double* logvtime_ik;

  const int* status_i = status;
  const int* status_ik;

  const int* Examiner_i = Examiner;
  const int* Examiner_ik;

  const int* Factor_i = Factor;
  const int* Factor_ik;
  
  const double* w_r;

  for (obs = 0; obs < *nP; obs++){

    /*** Calculate cumulative products based on specificities needed for 'A' numbers ***/
    cprod_spec_k = cprod_spec;
    *cprod_spec_k = 1.0;            /* k = 0*/
    cprod_spec_k++;

    status_ik   = status_i;
    Examiner_ik = Examiner_i;
    Factor_ik   = Factor_i;
    for (k = 1; k <= *nvisit_i; k++){
      *cprod_spec_k = *(cprod_spec_k - 1) * (*status_ik == 1 ? (1 - spec[*nFactor * *Examiner_ik + *Factor_ik]) : spec[*nFactor * *Examiner_ik + *Factor_ik]);
      cprod_spec_k++;
      status_ik++;
      Examiner_ik++;
      Factor_ik++;
    }
    
    /*** Calculate cumulative products based on sensitivities needed for 'A' numbers ***/
    cprod_sens_k = cprod_sens + *nvisit_i;
    *cprod_sens_k = 1.0;                    /* k = nvisit */
    cprod_sens_k--;

    status_ik--;
    Examiner_ik--;
    Factor_ik--;
    for (k = *nvisit_i - 1; k >= 0; k--){
      *cprod_sens_k = *(cprod_sens_k + 1) * (*status_ik == 1 ? sens[*nFactor * *Examiner_ik + *Factor_ik] : (1 - sens[*nFactor * *Examiner_ik + *Factor_ik]));
      cprod_sens_k--;
      status_ik--;
      Examiner_ik--;
      Factor_ik--;
    }    

    /*** Calculate the 'A' numbers for this observation ***/
    A_k          = A;
    cprod_sens_k = cprod_sens;
    cprod_spec_k = cprod_spec;
    for (k = 0; k <= *nvisit_i; k++){
      *A_k = *cprod_sens_k * *cprod_spec_k;
      A_k++;
      cprod_sens_k++;
      cprod_spec_k++;
    }

    /*** Calculate iPML for this observation ***/
    *iPML_i = 0.0;
    w_r  = gg_w;
    mu_r = *gg_gamma - *gg_K * *gg_delta;                     // = mu(0), see Gspline.h, function mu_component.
    for (r = 0; r <= 2 * *gg_K; r++){
      A_k = A;
      logvtime_ik  = logvtime_i;
      stres_cut_common_r = ((*min_eta_i) - *gg_intcpt - *gg_scale * mu_r) * invsigma_invscale;

      /** k = 0: first visit - like left-censored) **/
      stres_cut_k = stres_cut_common_r + *logvtime_ik * invsigma_invscale;
      Phi_cut_k = pnorm(stres_cut_k, 0, 1, 1, 0);
      *iPML_i += *A_k * *w_r * Phi_cut_k;

      A_k++;
      logvtime_ik++;

      /** k = 1, ..., *nvisit_i - 1: like interval-censored **/
      for (k = 1; k < *nvisit_i; k++){  
        Phi_cut_k_min_1 = Phi_cut_k;

        stres_cut_k = stres_cut_common_r + *logvtime_ik * invsigma_invscale;
        Phi_cut_k   = pnorm(stres_cut_k, 0, 1, 1, 0);
        *iPML_i += *A_k * *w_r * (Phi_cut_k - Phi_cut_k_min_1);

        A_k++;
        logvtime_ik++;
      }

      /** k = *nvisit_i: like right-censored **/
      *iPML_i += *A_k * *w_r * (1 - Phi_cut_k);
 
      /** Shift mixture pointers **/
      w_r++;
      mu_r += *gg_delta;                                     // = mu(r + 1), see Gspline.h, function mu_component.
    }

    /*** Shift pointers logvtime_i, status_i, Examiner_i, Factor_i at the same time.  ***/
    logvtime_i += *nvisit_i;
    status_i   += *nvisit_i;
    Examiner_i += *nvisit_i;
    Factor_i   += *nvisit_i;

    /*** Shift remaining pointers ***/
    min_eta_i++;
    nvisit_i++;
    iPML_i++;
  }
  
  return;
}    /*** end of function iPML_misclass_GJK ***/

}    // end of extern "C" {
