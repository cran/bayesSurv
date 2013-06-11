\name{bayesBisurvreg}
\alias{bayesBisurvreg}
\title{Population-averaged accelerated failure time model for bivariate,
  possibly doubly-interval-censored data. The error distribution is
  expressed as a~penalized bivariate normal mixture with high number
  of components (bivariate G-spline).
}
\description{A function to estimate a regression model with bivariate
  (possibly right-, left-, interval- or doubly-interval-censored) data.
  In the case of doubly interval censoring, different regression models
  can be specified for the onset and event times.
  
  The error density of the regression model is specified as a mixture of Bayesian G-splines
  (normal densities with equidistant means and constant variance matrices). This function performs an
  MCMC sampling from the posterior distribution of unknown quantities.

  For details, see Komárek (2006) and
  Komárek and Lesaffre (2006).
  
  We explain first in more detail a model without doubly censoring.
  Let \eqn{T_{i,l},\; i=1,\dots, N,\; l=1, 2}{T[i,l], i=1,..., N, l=1, 2}
  be event times for \eqn{i}{i}th cluster and the first and the second
  unit. The following regression model is assumed:
  \deqn{\log(T_{i,l}) = \beta'x_{i,l} + \varepsilon_{i,l},\quad i=1,\dots,N,\;l=1,2}{%
        log(T[i,l]) = beta'x[i,l] + epsilon[i,l], i=1,..., N, l=1,2}
  where \eqn{\beta}{beta} is unknown regression parameter vector and
  \eqn{x_{i,l}}{x[i,l]} is a vector of covariates. The bivariate error terms
  \eqn{\varepsilon_i=(\varepsilon_{i,1},\,\varepsilon_{i,2})',\;i=1,\dots,N}{epsilon[i] = (epsilon[i,1], epsilon[i,2])', i=1,..., N}
  are assumed to be i.i.d. with a~bivariate density
  \eqn{g_{\varepsilon}(e_1,\,e_2)}{g[epsilon](e[1], e[2])}. This density is expressed as
  a~mixture of Bayesian G-splines (normal densities with equidistant
  means and constant variance matrices). We distinguish two,
  theoretically equivalent, specifications.

  \describe{
    \item{Specification 1}{
      \deqn{(\varepsilon_1,\,\varepsilon_2)' \sim
	\sum_{j_1=-K_1}^{K_1}\sum_{j_2=-K_2}^{K_2} w_{j_1,j_2} N_2(\mu_{(j_1,j_2)},\,\mbox{diag}(\sigma_1^2,\,\sigma_2^2))}{%
        (epsilon[1],\,epsilon[2])' is distributed as
        sum[j[1]=-K[1]][K[1]] sum[j[2]=-K[2]][K[2]] w[j[1],j[2]]
	N(mu[(j[1],j[2])], diag(sigma[1]^2, sigma[2]^2))}
      where \eqn{\sigma_1^2,\,\sigma_2^2}{sigma[1]^2, sigma[2]^2} are
      \bold{unknown} basis variances and
      \eqn{\mu_{(j_1,j_2)} = (\mu_{1,j_1},\,\mu_{2,j_2})'}{%
	         mu[(j[1],j[2])] = (mu[1,j[1]],	mu[2,j[2]])'}
      is an~equidistant grid of knots symmetric around the
      \bold{unknown} point \eqn{(\gamma_1,\,\gamma_2)'}{(gamma[1], gamma[2])'} 
      and related to the unknown basis variances through the
      relationship
      \deqn{\mu_{1,j_1} = \gamma_1 + j_1\delta_1\sigma_1,\quad j_1=-K_1,\dots,K_1,}{%
	mu[1,j[1]] = gamma[1] + j[1]*delta[1]*sigma[1], j[1]=-K[1],...,	K[1]}
      \deqn{\mu_{2,j_2} = \gamma_2 + j_2\delta_2\sigma_2,\quad j_2=-K_2,\dots,K_2,}{%
	    mu[2,j[2]] = gamma[2] + j[2]*delta[2]*sigma[2], j[2]=-K[2],..., K[2]}
      where \eqn{\delta_1,\,\delta_2}{delta[1], delta[2]} are fixed
      constants, e.g. \eqn{\delta_1=\delta_2=2/3}{delta[1]=delta[2]=2/3}
      (which has a~justification of being close to cubic B-splines).
      \deqn{}{}
    }
    \item{Specification 2}{
      \deqn{(\varepsilon_1,\,\varepsilon_2)' \sim (\alpha_1,\,\alpha_2)'+ \bold{S}\,(V_1,\,V_2)'}{%
	(epsilon[1],\,epsilon[2])' is distributed as (alpha[1],	alpha[2])' + S (V[1], V[2])'}
      where \eqn{(\alpha_1,\,\alpha_2)'}{(alpha[1], alpha[2])'} is an
      \bold{unknown} intercept term and
      \eqn{\bold{S} \mbox{ is a diagonal matrix with } \tau_1 \mbox{ and }\tau_2 \mbox{ on a diagonal,}}{%
	S is a diagonal matrix with tau[1] and tau[2] on a diagonal,}
      i.e. \eqn{\tau_1,\,\tau_2}{tau[1], tau[2]} are \bold{unknown} scale
      parameters. \eqn{(V_1,\,V_2)')}{(V[1], V[2])'} is then
      standardized bivariate error term which is distributed according
      to the bivariate normal mixture, i.e.
      \deqn{(V_1,\,V_2)'\sim \sum_{j_1=-K_1}^{K_1}\sum_{j_2=-K_2}^{K_2}
	w_{j_1,j_2} N_2(\mu_{(j_1,j_2)},\,\mbox{diag}(\sigma_1^2, \sigma_2^2))}{%
	(V[1], V[2])' is distributed as sum[j[1]=-K[1]][K[1]]
	sum[j[2]=-K[2]][K[2]] w[j[1],j[2]] N(mu[(j[1],j[2])],
	diag(sigma[1]^2, sigma[2]^2))}
      where \eqn{\mu_{(j_1,j_2)} =
	(\mu_{1,j_1},\,\mu_{2,j_2})'}{mu[(j[1],j[2])] = (mu[1,j[1]],
	mu[2,j[2]])'} is an~equidistant grid of \bold{fixed} knots (means), usually
      symmetric about the \bold{fixed} point \eqn{(\gamma_1,\,\gamma_2)'=(0,
	0)'}{(gamma[1], gamma[2])' = (0, 0)'} and
      \eqn{\sigma_1^2,\,\sigma_2^2}{sigma[1]^2, sigma[2]^2} are
      \bold{fixed} basis variances. Reasonable values for the numbers of grid
      points \eqn{K_1}{K[1]} and \eqn{K_2}{K[2]} are
      \eqn{K_1=K_2=15}{K[1]=K[2]=15} with the distance between the two
      knots equal to \eqn{\delta=0.3}{delta=0.3} and for the basis
      variances
      \eqn{\sigma_1^2\sigma_2^2=0.2^2.}{sigma[1]^2=sigma[2]^2=0.2^2.}
    }  
  }
  Personally, I found Specification 2 performing better. In the paper
  Komárek and Lesaffre (2006) only
  Specification 2 is described.

  The mixture weights
  \eqn{w_{j_1,j_2},\;j_1=-K_1,\dots, K_1,\;j_2=-K_2,\dots,
    K_2}{w[j[1],j[2]], j[1]=-K[1],..., K[1], j[2]=-K[2],..., K[2]} are
  not estimated directly. To avoid the constraints
  \eqn{0 < w_{j_1,j_2} < 1}{0 < w[j[1],j[2]] < 1} and
  \eqn{\sum_{j_1=-K_1}^{K_1}\sum_{j_2=-K_2}^{K_2}w_{j_1,j_2} =
    1}{sum[j[1]=-K[1]][K[1]]sum[j[2]=-K[2]][K[2]]w[j[1],j[2]]=1}
  transformed weights  \eqn{a_{j_1,j_2},\;j_1=-K_1,\dots, K_1,\;j_2=-K_2,\dots,
    K_2}{a[j[1],j[2]], j[1]=-K[1],..., K[1], j[2]=-K[2],..., K[2]}
  related to the original weights by the logistic transformation:
  \deqn{a_{j_1,j_2} =
    \frac{\exp(w_{j_1,j_2})}{\sum_{m_1}\sum_{m_2}\exp(w_{m_1,m_2})}}{%
  a[j[1],j[2]] = exp(w[j[1],j[2]])/sum[m[1]]sum[m[2]] exp(w[m[1],m[2]])}
  are estimated instead.

  A~Bayesian model is set up for all unknown parameters. For more
  details I refer to Komárek and Lesaffre (2006)
  and to Komárek (2006).

  If there are doubly-censored data the model of the same type as above
  can be specified for both the onset time and the time-to-event. 
}
\usage{
bayesBisurvreg(formula, formula2, data = parent.frame(),
   na.action = na.fail, onlyX = FALSE,
   nsimul = list(niter = 10, nthin = 1, nburn = 0, nwrite = 10),
   prior, prior.beta, init = list(iter = 0),
   mcmc.par = list(type.update.a = "slice", k.overrelax.a = 1,
                   k.overrelax.sigma = 1, k.overrelax.scale = 1),
   prior2, prior.beta2, init2,
   mcmc.par2 = list(type.update.a = "slice", k.overrelax.a = 1,
                    k.overrelax.sigma = 1, k.overrelax.scale = 1),
   store = list(a = FALSE, a2 = FALSE, y = FALSE, y2 = FALSE,
                r = FALSE, r2 = FALSE),
   dir = getwd())
}
\arguments{
  \item{formula}{model formula for the regression. In the case of
    doubly-censored data, this is the model formula for the onset
    time. Data are assumed to be sorted according to subjects and within
    subjects according to the types of the events that determine the
    bivariate survival distribution, i.e. the response vector must be
    \eqn{t_{1,1}, t_{1,2}, t_{2,1}, t_{2,2}, t_{3,1}, t_{3,2}, \dots,
      t_{n,1},
      t_{n,2}}{t[1,1],t[1,2],t[2,1],t[2,2],t[3,1],t[3,2],...,t[n,1],t[n,2]}.
    The rows of the design matrix with covariates must be sorted
    analogically.

    The left-hand side of the formula must be an object created using
    \code{\link[survival]{Surv}}.
  }
  \item{formula2}{model formula for the regression of the time-to-event in
    the case of doubly-censored data. Ignored otherwise. The same remark as
    for \code{formula} concerning the sort order applies here.
  }
  \item{data}{optional data frame in which to interpret the variables
    occuring in the formulas.
  }
  \item{na.action}{the user is discouraged from changing the default
    value \code{na.fail}.}
  \item{onlyX}{if \code{TRUE} no MCMC sampling is performed and only the
      design matrix (matrices) are returned. This can be useful to set up
    correctly priors for regression parameters in the presence of
    \code{factor} covariates.}
  \item{nsimul}{a list giving the number of iterations of the MCMC and
    other parameters of the simulation.
    \describe{
      \item{niter}{total number of sampled values after discarding
	thinned ones, burn-up included;}
      \item{nthin}{thinning interval;}
      \item{nburn}{number of sampled values in a burn-up period after
	discarding thinned values. This value should be smaller than
	\code{niter}. If not, \code{nburn} is set to \code{niter - 1}. It can be set to zero;}
      \item{nwrite}{an interval at which information about the number of
	performed iterations is print on the screen and during the
	burn-up period an interval with which the sampled values are
	writen to files;}
    }    
  }    
  \item{prior}{a~list specifying the prior distribution of the G-spline
    defining the distribution of the error term in the regression model
    given by \code{formula}. See \code{prior} argument of
    \code{\link{bayesHistogram}} function for more detail. In this list
    also \sQuote{Specification} as described above is specified.
  }
  \item{prior.beta}{prior specification for the regression parameters,
    in the case of doubly censored data for the regression parameters of
    the onset time. I.e. it is related to \code{formula}.
    
    This should be a~list with the following components:
    \describe{
    \item{mean.prior}{a~vector specifying a~prior mean for each
      \code{beta} parameter in the model.}
    \item{var.prior}{a~vector specifying a~prior variance for each
      \code{beta} parameter.}
    }
    It is recommended to run the function
    bayesBisurvreg first with its argument \code{onlyX} set to \code{TRUE}
    to find out how the betas are sorted. They must correspond to a
    design matrix X taken from \code{formula}.
  }  %% end of item{prior.beta}
  \item{init}{an~optional list with initial values for the MCMC related
    to the model given by \code{formula}. The list can have the following components:
    \describe{
    \item{iter}{the number of the iteration to which the initial values
      correspond, usually zero.}
    \item{beta}{a~vector of initial values for the regression
      parameters. It must be sorted in the same way as are the columns
      in the design matrix. Use \code{onlyX=TRUE} if you do not know how
      the columns in the design matrix are created.}
    \item{a}{a~matrix of size \eqn{(2K_1+1)\times(2K_2+1)}{(2*K[1]+1) x
	(2*K[2]+1)} with the initial values of transformed mixture weights.}
    \item{lambda}{initial values for the Markov random fields precision
      parameters. According to the chosen prior for the transformed
      mixture weights, this is either a~number or a~vector of length 2.
    }
    \item{gamma}{a~vector of length 2 of initial values for the middle
      knots \eqn{\gamma_1, \gamma_2}{gamma[1], gamma[2]} in each
      dimension.
      
      If \sQuote{Specification} is 2, this value will not be changed
      by the MCMC and it is recommended (for easier
      interpretation of the results) to set \code{init$gamma} to zero
      for all dimensions (default behavior).
      
      If \sQuote{Specification} is 1 \code{init$gamma} should be
      approximately equal to the mean value of the residuals in each
      margin.
    }
    \item{sigma}{a~vector of length~2 of initial values of the basis
      standard deviations \eqn{\sigma_1, \sigma_2}{sigma[1], sigma[2]}.
      
      If \sQuote{Specification} is 2 this value will not be changed
        by the MCMC and it is recommended to set it
	approximately equal to the range of standardized data (let say 4
	+ 4) divided by the number of knots in each margin and
	multiplied by something like 2/3.

	If \sQuote{Specification} is 1
	this should be approximately equal to the range of the residuals
	divided by the number of knots in each margin and
	multiplied again by something like 2/3.
      }  
    \item{intercept}{a~vector of length~2 of initial values of the
      intercept terms \eqn{\alpha_1, \alpha_2}{alpha[1], alpha[2]}.
      
      If \sQuote{Specification} is 1 this value is not changed by the
      MCMC and the initial value is always changed to zero for both dimensions.}
    \item{scale}{a~vector of length~2 of initial values of the scale
      parameters \eqn{\tau_1, \tau_2}{tau[1], tau[2]}.

      If \sQuote{Specification} is 1 this value is not changed by the MCMC
      and the initial value is always changed to one for both dimensions.}
    \item{y}{a~matrix with 2 columns and \eqn{N}{N} rows with initial
      values of log-event-times for each cluster in rows.}
    \item{r}{a~matrix with 2 columns and \eqn{N}{N} rows with initial
      component labels for each bivariate residual in rows. All values
      in the first column must be between \eqn{-K_1}{-K[1]}} and
    \eqn{K_1}{K[1]} and all values in the second column must be between
    \eqn{-K_2}{-K[2]} and \eqn{K_2}{K[2]}. See argument \code{init} of
    the function \code{\link{bayesHistogram}} for more details.}
  }  
  \item{mcmc.par}{a~list specifying how some of the G-spline parameters
    related to \code{formula} are to be updated. The list can have the following components (all
    of them have their default values):
    \describe{
      \item{type.update.a}{G-spline transformed weights \eqn{a}{a} can
	be updated using one of the following algorithms:
        \describe{
	  \item{slice}{slice sampler of Neal (2003)}
	  \item{ars.quantile}{adaptive rejection sampling of Gilks and
	    Wild (1992) with starting abscissae being quantiles of the
	    envelop at the previous iteration}
	  \item{ars.mode}{adaptive rejection sampling of Gilks and Wild
	    (1992) with starting abscissae being the mode plus/minus 3
	    times estimated standard deviation of the full conditional
	    distribution}
        }
	Default is \code{slice}.
      }
      \item{k.overrelax.a}{if \code{type.update.a == "slice"} some
	updates are overrelaxed. Then every \code{k.overrelax.a}th
	iteration is not overrelaxed. Default is \code{k.overrelax.a =
	  1}, i.e. no overrelaxation}
      \item{k.overrelax.sigma}{G-spline basis standard deviations are
	updated using the slice sampler of Neal (2003). At the same
	time, overrelaxation can be used. Then every
	k.overrelax.sigma th update is not overrelaxed. Default is
	\code{k.overrelax.sigma = 1}, i.e. no overrelaxation}
      \item{k.overrelax.scale}{G-spline scales are
	updated using the slice sampler of Neal (2003). At the same
	time, overrelaxation can be used. Then every
	k.overrelax.scale th update is not overrelaxed. Default is
      \code{k.overrelax.scale = 1}, i.e. no overrelaxation}
    }  %% end of describe  
  }  %% end of item{mcmc.par}
  \item{prior2}{a~list specifying the prior distribution of the G-spline
    defining the distribution of the error term in the regression model
    given by \code{formula2}. See \code{prior} argument of
    \code{\link{bayesHistogram}} function for more detail.  }
  \item{prior.beta2}{prior specification for the regression parameters
    of time-to-event in the case of doubly censored data (related to
    \code{formula2}).
    This should be a~list with the same structure as \code{prior.beta}.
  }  %% end of item{prior.beta2}
  \item{init2}{an~optional list with initial values for the MCMC related
    to the model given by \code{formula2}. The list has the same
    structure as \code{init}.
  }
  \item{mcmc.par2}{a~list specifying how some of the G-spline parameters
    related to \code{formula2} are to be updated. The list has the same
    structure as \code{mcmc.par}.
  }
  \item{store}{a~list of logical values specifying which chains that are
    not stored by default are to be stored. The list can have the
    following components.
    \describe{
      \item{a}{if \code{TRUE} then all the transformed mixture weights
	\eqn{a_{k_1,\,k_2},}{a[k[1],k[2]],}
	\eqn{k_1=-K_1,\dots,K_1,}{k[1]=-K[1],..., K[1],}
	\eqn{k_2=-K_2,\dots,K_2,}{k[2]=-K[2],..., K[2],}
	related to the G-spline of \code{formula} are stored.}
      \item{a2}{if \code{TRUE} and there are doubly-censored data then all the transformed mixture weights
	\eqn{a_{k_1,\,k_2},}{a[k[1],k[2]],}
	\eqn{k_1=-K_1,\dots,K_1,}{k[1]=-K[1],..., K[1],}
	\eqn{k_2=-K_2,\dots,K_2,}{k[2]=-K[2],..., K[2],}
	related to the G-spline of \code{formula2} are stored.}
      \item{y}{if \code{TRUE} then augmented log-event times for all
	observations related to the \code{formula} are stored.}
      \item{y2}{if \code{TRUE} then augmented log-event times for all
	observations related to \code{formula2} are stored.}
      \item{r}{if \code{TRUE} then labels of mixture components for
	residuals related to \code{formula} are stored.}
      \item{r2}{if \code{TRUE} then labels of mixture components for
	residuals related to \code{formula2} are stored.}
    }  %% end of describe  
  }  %% end of item{store}
  \item{dir}{a string that specifies a directory where all sampled
    values are to be stored.
  }  
}
\value{
  A list of class \code{bayesBisurvreg} containing an information
  concerning the initial values and prior choices.
}
\section{Files created}{  
  Additionally, the following files with sampled values
  are stored in a directory specified by \code{dir} argument of this
  function (some of them are created only on request, see \code{store}
  parameter of this function).

  Headers are written to all files created by default and to files asked
  by the user via the argument \code{store}. During the burn-in, only
  every \code{nsimul$nwrite} value is written. After the burn-in, all
  sampled values are written in files created by default and to files
  asked by the user via the argument \code{store}. In the files for
  which the corresponding \code{store} component is \code{FALSE}, every
  \code{nsimul$nwrite} value is written during the whole MCMC (this
  might be useful to restart the MCMC from some specific point).
  
  The following files are created:
  \describe{
    \item{iteration.sim}{one column labeled \code{iteration} with
      indeces of MCMC iterations to which the stored sampled values
      correspond.
    }
    \item{mixmoment.sim}{columns labeled \code{k}, \code{Mean.1}, \code{Mean.2},
      \code{D.1.1}, \code{D.2.1}, \code{D.2.2}, where
      
      \bold{k} = number of mixture components that had probability
      numerically higher than zero;
      
      \bold{Mean.1} =
      \eqn{\mbox{E}(\varepsilon_{i,1})}{E(epsilon[i,1])};
      
      \bold{Mean.2} =
      \eqn{\mbox{E}(\varepsilon_{i,2})}{E(epsilon[i,2])};
      
      \bold{D.1.1} =
      \eqn{\mbox{var}(\varepsilon_{i,1})}{var(epsilon[i,1])};
      
      \bold{D.2.1} =
      \eqn{\mbox{cov}(\varepsilon_{i,1},\,\varepsilon_{i,2})}{cov(epsilon[i,1], epsilon[i,2])};

      \bold{D.2.2} =
      \eqn{\mbox{var}(\varepsilon_{i,2})}{var(epsilon[i,2])};
      
      all related to the distribution of the error term from the model given by \code{formula}.
    }
    \item{mixmoment\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{mixmoment.sim}, however related to the model
      given by \code{formula2}.      
    }
    \item{mweight.sim}{sampled mixture weights
      \eqn{w_{k_1,\,k_2}}{w[k[1],k[2]]} of mixture components that had
      probabilities numerically higher than zero. Related to the model
      given by \code{formula}.
    }
    \item{mweight\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{mweight.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{mmean.sim}{indeces \eqn{k_1,\;k_2,}{k[1], k[2],}
      \eqn{k_1 \in\{-K_1, \dots, K_1\},}{k[1] in {-K[1], ..., K[1]},}
      \eqn{k_2 \in\{-K_2, \dots, K_2\}}{k[2] in {-K[2], ..., K[2]}}
      of mixture components that had probabilities numerically higher
      than zero. It corresponds to the weights in
      \code{mweight.sim}. Related to the model given by \code{formula}.
    }
    \item{mmean\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{mmean.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{gspline.sim}{characteristics of the sampled G-spline
      (distribution of
      \eqn{(\varepsilon_{i,1},\,\varepsilon_{i,2})'}{(epsilon[i,1],
	epsilon[i,2])'}) related to the model given by
      \code{formula}. This file together with \code{mixmoment.sim},
      \code{mweight.sim} and \code{mmean.sim} can be used to reconstruct
      the G-spline in each MCMC iteration.
      
      The file has columns labeled \code{gamma1},
      \code{gamma2}, \code{sigma1}, \code{sigma2}, \code{delta1},
      \code{delta2}, \code{intercept1}, \code{intercept2},
      \code{scale1}, \code{scale2}. The meaning of the values in these
      columns is the following:
      
      \bold{gamma1} = the middle knot \eqn{\gamma_1}{gamma[1]} in the
      first dimension. If \sQuote{Specification} is 2, this column
      usually contains zeros;
      
      \bold{gamma2} = the middle knot \eqn{\gamma_2}{gamma[2]} in the
      second dimension. If \sQuote{Specification} is 2, this column
      usually contains zeros;
      
      \bold{sigma1} = basis standard deviation \eqn{\sigma_1}{sigma[1]}
      of the G-spline in the first dimension. This column contains
      a~fixed value if \sQuote{Specification} is 2;
      
      \bold{sigma2} = basis standard deviation \eqn{\sigma_2}{sigma[2]}
      of the G-spline in the second dimension. This column contains
      a~fixed value if \sQuote{Specification} is 2;
      
      \bold{delta1} = distance \eqn{delta_1}{delta[1]} between the two knots of the G-spline in
      the first dimension. This column contains
      a~fixed value if \sQuote{Specification} is 2;
      
      \bold{delta2} = distance \eqn{\delta_2}{delta[2]} between the two knots of the G-spline in
      the second dimension. This column contains a~fixed value if
      \sQuote{Specification} is 2;
      
      \bold{intercept1} = the intercept term \eqn{\alpha_1}{alpha[1]} of
      the G-spline in the first dimension. If \sQuote{Specification} is 1, this column
      usually contains zeros;
      
      \bold{intercept2} = the intercept term \eqn{\alpha_2}{alpha[2]} of
      the G-spline in the second dimension. If \sQuote{Specification} is 1, this column
      usually contains zeros;
      
      \bold{scale1} = the scale parameter \eqn{\tau_1}{tau[1]} of the
      G-spline in the first dimension. If \sQuote{Specification} is 1, this column
      usually contains ones;
      
      \bold{scale2} = the scale parameter \eqn{\tau_2}{tau[2]} of the
      G-spline in the second dimension. \sQuote{Specification} is 1, this column
      usually contains ones.      
    }
    \item{gspline\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{gspline.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{mlogweight.sim}{fully created only if \code{store$a = TRUE}. The
      file contains the transformed weights
      \eqn{a_{k_1,\,k_2},}{a[k[1],k[2]],}
      \eqn{k_1=-K_1,\dots,K_1,}{k[1]=-K[1],..., K[1],}
      \eqn{k_2=-K_2,\dots,K_2}{k[2]=-K[2],..., K[2]} of all mixture
      components, i.e. also of components that had numerically zero
      probabilities.
      This file is related to the model given by \code{formula}.
    }
    \item{mlogweight\_{}2.sim}{fully created only if \code{store$a2 =
      TRUE} and in the case of doubly-censored data, the same
      structure as \code{mlogweight.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{r.sim}{fully created only if \code{store$r = TRUE}. The file
      contains the labels of the mixture components into which the
      residuals are intrinsically assigned. Instead of double indeces
      \eqn{(k_1,\,k_2)}{(k[1], k[2])}, values from 1 to \eqn{(2\,K_1+1)\times
	(2\,K_2+1)}{(2*K[1]+1)*(2*K[2]+1)} are stored here. Function
      \code{\link{vecr2matr}} can be used to transform it back to double
      indeces.
    }
    \item{r\_{}2.sim}{fully created only if \code{store$r2 =
      TRUE} and in the case of doubly-censored data, the same
      structure as \code{r.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{lambda.sim}{either one column labeled \code{lambda} or two
      columns labeled \code{lambda1} and \code{lambda2}. These are the
      values of the smoothing parameter(s) \eqn{\lambda}{lambda}
      (hyperparameters of the prior distribution of the transformed
      mixture weights \eqn{a_{k_1,\,k_2}}{a[k[1],k[2]]}). This file is
      related to the model given by \code{formula}.
    }
    \item{lambda\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{lambda.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{beta.sim}{sampled values of the regression parameters
      \eqn{\beta}{beta} related to the model given by
      \code{formula}. The columns are labeled according to the
      \code{colnames} of the design matrix.
    }
    \item{beta\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{beta.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{Y.sim}{fully created only if \code{store$y = TRUE}. It
      contains sampled (augmented) log-event times for all observations
      in the data set.
    }
    \item{Y\_{}2.sim}{fully created only if \code{store$y2 =
      TRUE} and in the case of doubly-censored data, the same
      structure as \code{Y.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{logposter.sim}{columns labeled \code{loglik}, \code{penalty}
      or \code{penalty1} and
      \code{penalty2}, \code{logprw}. This file is related to the model
      given by \code{formula}. The columns have the following meaning.

      \bold{loglik}
      \eqn{=}{=} \eqn{%
	-N\Bigl\{\log(2\pi) + \log(\sigma_1) + \log(\sigma_2)\Bigr\}-
          0.5\sum_{i=1}^N\Bigl\{
	  (\sigma_1^2\,\tau_1^2)^{-1}\; (y_{i,1} - x_{i,1}'\beta - \alpha_1 - \tau_1\mu_{1,\,r_{i,1}})^2 +
	  (\sigma_2^2\,\tau_2^2)^{-1}\; (y_{i,2} - x_{i,2}'\beta - \alpha_2 - \tau_2\mu_{2,\,r_{i,2}})^2
          \Bigr\}
      }{%
	-N(log(2*pi) + log(sigma[1]) + log(sigma[2]))
          -0.5*sum[i=1][N](
          (sigma[1]^2*tau[1]^2)^(-1) * (y[i,1] - x[i,1]'beta - alpha[1] - tau[1]*mu[1,r[i,1]])^2 +
          (sigma[2]^2*tau[2]^2)^(-1) * (y[i,2] - x[i,2]'beta - alpha[2] - tau[2]*mu[2,r[i,2]])^2 	  
	)}
      
      where \eqn{y_{i,l}}{y[i,l]} denotes (augmented) \emph{(i,l)}th
      true log-event time. In other words, \code{loglik} is equal to the
      conditional log-density
      \eqn{\sum_{i=1}^N\,\log\Bigl\{p\bigl((y_{i,1},\,y_{i,2})\;\big|\;r_{i},\,\beta,\,\mbox{G-spline}\bigr)\Bigr\};}{%
      sum[i=1][N] log(p((y[i,1], y[i,2]) | r[i], beta, G-spline));
      }

      \bold{penalty1:} If \code{prior$neighbor.system} = \code{"uniCAR"}:
        the penalty term for the first dimension not multiplied by
	\code{lambda1};
      
      \bold{penalty2:} If \code{prior$neighbor.system} = \code{"uniCAR"}:
        the penalty term for the second dimension not multiplied by
	\code{lambda2};

      \bold{penalty:} If \code{prior$neighbor.system} is different from \code{"uniCAR"}:
        the penalty term not multiplied by \code{lambda};

      \bold{logprw} \eqn{=}{=}
      \eqn{-2\,N\,\log\bigl\{\sum_{k_1}\sum_{k_2}a_{k_1,\,k_2}\bigr\} +
	\sum_{k_1}\sum_{k_2}N_{k_1,\,k_2}\,a_{k_1,\,k_2},}{%
	-2*N*log(sum[k[1]]sum[k[2]] exp(a[k[1],k[2]])) +
	sum[k[1]]sum[k[2]] N[k[1],k[2]]*a[k[1],k[2]],}
      where \eqn{N_{k_1,\,k_2}}{N[k[1],k[2]]} is the number of residuals
      assigned intrinsincally to the \eqn{(k_1,\,k_2)}{(k[1], k[2])}th
      mixture component.

      In other words, \code{logprw} is equal to the conditional
      log-density
      \eqn{\sum_{i=1}^N \log\bigl\{p(r_i\;|\;\mbox{G-spline
	  weights})\bigr\}.}{%
      sum[i=1][N] log(p(r[i] | G-spline weights)).}
    }
    \item{logposter\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{lambda.sim}, however related to the model
      given by \code{formula2}. 
    }
  }
}
\references{
  Gilks, W. R. and Wild, P. (1992).
  Adaptive rejection sampling for Gibbs sampling.
  \emph{Applied Statistics,} \bold{41}, 337 - 348.

  Komárek, A. (2006).
  \emph{Accelerated Failure Time Models for Multivariate
    Interval-Censored Data with Flexible Distributional Assumptions}.
  PhD. Thesis, Katholieke Universiteit Leuven, Faculteit Wetenschappen.

  Komárek, A. and Lesaffre, E. (2006).
  Bayesian semi-parametric accelerated failure time model for paired
  doubly interval-censored data.
  \emph{Statistical Modelling}, \bold{6}, 3--22.
  
  Neal, R. M. (2003).
  Slice sampling (with Discussion).
  \emph{The Annals of Statistics,} \bold{31}, 705 - 767.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\examples{
## See the description of R commands for
## the population averaged AFT model
## with the Signal Tandmobiel data,
## analysis described in Komarek and Lesaffre (2006),
##
## R commands available in the documentation
## directory of this package as
## - see ex-tandmobPA.R and
##   http://www.karlin.mff.cuni.cz/~komarek/software/bayesSurv/ex-tandmobPA.pdf
##
}
\keyword{survival}
\keyword{regression}
\keyword{multivariate}
