// Function to update possibly censored data distributed as a G-spline
//
// 22/10/2004: 'update_Data_GS'
// 10/12/2004: 2nd specification of the G-spline implemented in 'update_Data_GS'
// 12/01/2005: 'update_Data_GS_doubly'
//             'update_Data_GS_regres'
// 20/01/2005: 'update_Data_GS_doubly' rewritten
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
               const double* y_left,  const double* y_right,   const int* status,
               const int* rM,         const Gspline* gg,
               const int* nP,         const int* n_censored)
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


// ****** update_Data_GS ***********************
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
update_Data_GS_regres(double* YsM,           double* regresResM,
                      const double* y_left,  const double* y_right,   const int* status,
                      const int* rM,         const Gspline* gg,       const int* nP)
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
update_Data_GS_doubly(double* Yevent,        double* regresResM,
                      const double* Yonset, 
  	              const double* t_left,  const double* t_right,   const int* status,
                      const int* rM,         const Gspline* gg,       const int* nP)
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
