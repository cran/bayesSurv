\name{bayessurvreg2}
\alias{bayessurvreg2}
\title{
  Cluster-specific accelerated failure time model for multivariate,
  possibly doubly-interval-censored data. The error distribution is
  expressed as a~penalized univariate normal mixture with high number of
  components (G-spline). The distribution of the vector of random
  effects is multivariate normal.
}
\description{
  A function to estimate a regression model with possibly clustered
  (possibly right, left, interval or doubly-interval censored) data.
  In the case of doubly-interval censoring, different regression models
  can be specified for the onset and event times.

  (Multivariate) random effects, normally distributed and  acting as
  in the linear mixed model, normally distributed, can be included to
  adjust for clusters.

  The error density of the regression model is specified as a mixture of
  Bayesian G-splines (normal densities with equidistant means and
  constant variances). This function performs an MCMC sampling from the
  posterior distribution of unknown quantities.

  For details, see \eqn{\mbox{Kom\'{a}rek}}{Komarek} (2006),
  and \eqn{\mbox{Kom\'{a}rek},}{Komarek,} Lesaffre and Legrand (2006).

  We explain first in more detail a model without doubly censoring.
  Let \eqn{T_{i,l},\; i=1,\dots, N,\; l=1,\dots, n_i}{%
           T[i,l], i=1,..., N, l=1,..., n[i]}
  be event times for \eqn{i}{i}th cluster and the units within that cluster
  The following regression model is assumed:
  \deqn{\log(T_{i,l}) = \beta'x_{i,l} + b_i'z_{i,l} + \varepsilon_{i,l},\quad i=1,\dots, N,\;l=1,\dots, n_i}{%
    log(T[i,l]) = beta'x[i,l] + b[i]'z[i,l] + epsilon[i,l], i=1,..., N, l=1,..., n[i]}
  where \eqn{\beta}{beta} is unknown regression parameter vector,
  \eqn{x_{i,l}}{x[i,l]} is a vector of covariates.
  \eqn{b_i}{b[i]} is a (multivariate) cluster-specific random effect
  vector and \eqn{z_{i,l}}{z[i,l]} is a vector of covariates for random
  effects.

  The random effect vectors \eqn{b_i,\;i=1,\dots, N}{b[i], i=1,..., N}
  are assumed to be i.i.d. with a (multivariate) normal distribution
  with the mean \eqn{\beta_b}{beta[b]} and a~covariance matrix
  \eqn{D}{D}. Hierarchical centring (see Gelfand, Sahu, Carlin, 1995) is
  used. I.e. \eqn{\beta_b}{beta[b]} expresses the average effect of the
  covariates included in \eqn{z_{i,l}}{z[i,l]}. Note that covariates
  included in \eqn{z_{i,l}}{z[i,l]} may not be included in the covariate
  vector \eqn{x_{i,l}}{x[i,l]}. The covariance matrix \eqn{D}{D} is
  assigned an inverse Wishart prior distribution in the next level of hierarchy.
  
  The error terms
  \eqn{\varepsilon_{i,l},\;i=1,\dots, N, l=1,\dots, n_i}{%
       epsilon[i,l], i=1,..., N, l=1,..., n[i]}
  are assumed to be i.i.d. with a~univariate density
  \eqn{g_{\varepsilon}(e)}{g[epsilon](e)}. This density is expressed as
  a~mixture of Bayesian G-splines (normal densities with equidistant
  means and constant variances). We distinguish two,
  theoretically equivalent, specifications.

  \describe{
    \item{Specification 1}{
      \deqn{\varepsilon \sim
	\sum_{j=-K}^{K} w_{j} N(\mu_{j},\,\sigma^2)}{%
        epsilon is distributed as
        sum[j=-K][K] w[j]
	N(mu[j], sigma^2)}
      where \eqn{\sigma^2}{sigma^2} is the
      \bold{unknown} basis variance and
      \eqn{\mu_{j},\;j=-K,\dots, K}{mu[j],\;j=-K,..., K}
      is an~equidistant grid of knots symmetric around the
      \bold{unknown} point \eqn{\gamma}{gamma} 
      and related to the unknown basis variance through the
      relationship
      \deqn{\mu_{j} = \gamma + j\delta\sigma,\quad j=-K,\dots,K,}{%
	mu[j] = gamma + j*delta*sigma, j=K,...,	K}
      where \eqn{\delta}{delta} is fixed
      constants, e.g. \eqn{\delta=2/3}{delta=2/3}
      (which has a~justification of being close to cubic B-splines).
      \deqn{}{}
    }
    \item{Specification 2}{
      \deqn{\varepsilon \sim \alpha + \tau\,V}{%
	epsilon[1] is distributed as alpha + tau * V}
      where \eqn{\alpha}{alpha} is an
      \bold{unknown} intercept term and
      \eqn{\tau}{tau} is an \bold{unknown} scale parameter.
      \eqn{V}{V} is then
      standardized error term which is distributed according
      to the univariate normal mixture, i.e.
      \deqn{V\sim \sum_{j=-K}^{K}
	w_{j} N(\mu_{j},\,\sigma^2)}{%
	V is distributed as sum[j=-K][K]
	w[j] N(mu[j], sigma^2)}
      where \eqn{\mu_{j},\;j=-K,\dots, K}{mu[j], j=-K,..., K}
      is an~equidistant grid of \bold{fixed} knots (means), usually
      symmetric about the \bold{fixed} point \eqn{\gamma=0}{gamma = 0} and
      \eqn{\sigma^2}{sigma^2} is \bold{fixed} basis variance.
      Reasonable values for the numbers of grid
      points \eqn{K}{K} is
      \eqn{K=15}{K=15} with the distance between the two
      knots equal to \eqn{\delta=0.3}{delta=0.3} and for the basis
      variance
      \eqn{\sigma^2=0.2^2.}{sigma^2=0.2^2.}
    }  
  }
  Personally, I found Specification 2 performing better. In the paper
  \eqn{\mbox{Kom{\'a}rek}}{Komarek} and Lesaffre (2006) only
  Specification 2 is described.

  The mixture weights
  \eqn{w_{j},\;j=-K,\dots, K}{w[j], j=-K,..., K} are
  not estimated directly. To avoid the constraints
  \eqn{0 < w_{j} < 1}{0 < w[j] < 1} and
  \eqn{\sum_{j=-K}^{K}\,w_j = 1}{sum[j=-K][K] w[j] = 1}
  transformed weights \eqn{a_{j},\;j=-K,\dots, K}{a[j], j=-K,..., K}
  related to the original weights by the logistic transformation:
  \deqn{a_{j} = \frac{\exp(w_{j})}{\sum_{m}\exp(w_{m})}}{%
        a[j] = exp(w[j])/sum[m] exp(w[m])}
  are estimated instead.

  A~Bayesian model is set up for all unknown parameters. For more
  details I refer to \eqn{\mbox{Kom\'{a}rek}}{Komarek} (2006)
  and to
  \eqn{\mbox{Kom\'{a}rek},}{Komarek,} Lesafre, and Legrand (2006).
      
  If there are doubly-censored data the model of the same type as above
  can be specified for both the onset time and the time-to-event. 
}
\usage{
bayessurvreg2(formula, random, formula2, random2,
   data = parent.frame(),
   na.action = na.fail, onlyX = FALSE,
   nsimul = list(niter = 10, nthin = 1, nburn = 0, nwrite = 10),   
   prior, prior.beta, prior.b, init = list(iter = 0),
   mcmc.par = list(type.update.a = "slice", k.overrelax.a = 1,
                   k.overrelax.sigma = 1, k.overrelax.scale = 1),
   prior2, prior.beta2, prior.b2, init2,
   mcmc.par2 = list(type.update.a = "slice", k.overrelax.a = 1,
                    k.overrelax.sigma = 1, k.overrelax.scale = 1),                           
   store = list(a = FALSE, a2 = FALSE, y = FALSE, y2 = FALSE,
                r = FALSE, r2 = FALSE, b = FALSE, b2 = FALSE), 
   dir = getwd())
}
\arguments{
  \item{formula}{model formula for the regression. In the case of
    doubly-censored data, this is the model formula for the onset
    time.

    The left-hand side of the \code{formula} must be an~object created
    using \code{\link[survival]{Surv}}.

    In the formula all covariates appearing both in the vector
    \eqn{x_{i,l}}{x[i,l]} and \eqn{z_{i,l}}{z[i,l]} must be
    mentioned. Intercept is implicitely included in the model by the
    estimation of the error distribution. As a~consequence \code{-1} in
    the model formula does not have any effect on the model
    specification.

    If \code{random} is used then the formula must contain
    an identification of clusters in the form \code{cluster(id)}, where
    \code{id} is a name of the variable that determines clusters, e.g.
    \tabular{c}{
      \code{Surv(time, event)~gender + cluster(id)}.
    }     
  }  %% end of item{formula}
  \item{random}{formula for the `random' part of the model, i.e. the
    part that specifies the covariates \eqn{z_{i,l}}{z[i,l]}. In the
    case of doubly-censored data, this is the \code{random} formula for
    the onset time.

    If omitted, no random part is included in the model. E.g. to specify the model with a
    random intercept, say \code{random=~1}. All effects mentioned in
    \code{random} should also be mentioned on the right-hand side of
    \code{formula}.

    When some random effects are included the random intercept is added
    by default. It can be removed using e.g. \code{random=~-1 + gender}.
  }  %% end of item{random}
  \item{formula2}{model formula for the regression of the time-to-event in
    the case of doubly-censored data. Ignored otherwise. The same structure as
    for \code{formula} applies here.
  }  %% end of item{formula2}
  \item{random2}{specification of the `random' part of the model for
    time-to-event in the case of doubly-censored data. Ignored
    otherwise. The same structure as for \code{random} applies here.
  }  %% end of item{random2}
  \item{data}{optional data frame in which to interpret the variables
    occuring in the \code{formula}, \code{formula2}, \code{random},
    \code{random2} statements.
  }  %% end of item{data}
  \item{na.action}{the user is discouraged from changing the default
    value \code{na.fail}.
  }  %% end of item{na.action}
  \item{onlyX}{if \code{TRUE} no MCMC sampling is performed and only the
    design matrix (matrices) are returned. This can be useful to set up
    correctly priors for regression parameters in the presence of
    \code{factor} covariates.
  }  %% end of item{onlyX}
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
  }  %% end of item{nsimul}
  \item{prior}{a~list specifying the prior distribution of the G-spline
    defining the distribution of the error term in the regression model
    given by \code{formula} and \code{random}. See \code{prior} argument of
    \code{\link{bayesHistogram}} function for more detail. In this list
    also \sQuote{Specification} as described above is specified.

    The item \code{prior$neighbor.system} can only be equal to
    \code{uniCAR} here.
  }  %% end of item{prior}
  \item{prior.b}{a list defining the way in which the random effects
    involved in \code{formula} and \code{random}
    are to be updated and the specification
    of priors for parameters related to these random effects. The list
    is assumed to have the following components.
    \describe{
    \item{prior.D}{a string defining the prior distribution for the
      covariance matrix of random effects \eqn{D}{D}. It can be either
      ``inv.wishart'' or ``sduniform''.
      \describe{
      \item{inv.wishart}{in that case is assumed that the prior distribution
        of the matrix \eqn{D}{D} is Inverse-Wishart with degrees of freedom
        equal to \eqn{\tau}{tau} and a scale matrix equal to
        \eqn{S}{S}. When \eqn{D}{D} is a matrix \eqn{q\times q}{q x q} a
        prior expectation of \eqn{D}{D} is equal to
        \eqn{(\tau - q - 1)^{-1}S}{(1/(tau - q - 1))S} if
        \eqn{\tau > q + 1}{tau > q + 1}. For
        \eqn{q - 1 < \tau \leq q + 1}{q - 1 < tau <= q + 1} a prior
        expectation is not finite.	
        Degrees of freedom parameter \eqn{\tau}{tau} does not have to be an
        integer. It has to only satisfy a condition
        \eqn{\tau > q - 1}{tau > q - 1}. \code{prior.b$df.D} gives a prior
        degrees of freedom parameter \eqn{\tau}{tau} and
        \code{prior.b$scale.D} determines the scale matrix \eqn{D}{D}.
        Inverse-Wishart is also the default choice.
      }
      \item{sduniform}{this can be used only when the random effect is
	univariate (e.g. only random intercept in the model).
	Then the matrix \eqn{D}{D} is just a scalar and the
	prior of \eqn{\sqrt{D}}{sqrt(D)} (standard deviation of the
	univariate random effect) is assumed to be uniform on interval
	\eqn{(0, S)}{(0, S)}. The upper limit \eqn{S}{S} is given by \code{prior.b$scale.D}.
      }	
      }  % end of describe
    }  % end of item(Prior.D)
    \item{df.D}{degrees of freedom parameter \eqn{\tau}{tau} in the case
      that the prior of the matrix \eqn{D}{D} is inverse-Wishart.
      }
    \item{scale.D}{a lower triangle of the scale matrix \eqn{S}{S} in
      the case that the prior of the matrix \eqn{D}{D} is
      inverse-Wishart or the upper limit \eqn{S}{S} of the uniform distribution in
      the case that \eqn{\sqrt{D} \sim \mbox{Unif}(0, S)}{sqrt(D) ~ Unif(0, S)}.
      }    
    }  %% end of describe
  }  %% end of prior.b  
  \item{prior.beta}{prior specification for the regression parameters,
    in the case of doubly-censored data for the regression parameters of
    the onset time, i.e. it is related to \code{formula} and
    \code{random}.
    Note that the \code{beta} vector contains both the fixed effects
    \eqn{\beta}{beta} and the means of the random effects (except the
    random intercept) \eqn{\beta_b}{beta[b]}.
    
    This should be a~list with the following components:
    \describe{
    \item{mean.prior}{a~vector specifying a~prior mean for each
      \code{beta} parameter in the model.}
    \item{var.prior}{a~vector specifying a~prior variance for each
      \code{beta} parameter.}
    }
    It is recommended to run the function
    bayessurvreg2 first with its argument \code{onlyX} set to \code{TRUE}
    to find out how the betas are sorted. They must correspond to a
    design matrix X taken from \code{formula}.
  }  %% end of item{prior.beta}
  \item{init}{an~optional list with initial values for the MCMC related
    to the model given by \code{formula} and \code{random}. The list can have the following components:
    \describe{
    \item{iter}{the number of the iteration to which the initial values
      correspond, usually zero.}
    \item{beta}{a~vector of initial values for the regression
      parameters (both the fixed effects and means of the random
      effects).
      It must be sorted in the same way as are the columns
      in the design matrix. Use \code{onlyX=TRUE} if you do not know how
      the columns in the design matrix are created.}
    \item{a}{a~vector of length \eqn{2K+1}{2*K+1} with the initial values of transformed mixture weights.}
    \item{lambda}{initial values for the Markov random fields precision parameter. 
    }
    \item{gamma}{an~initial values for the middle
      knot \eqn{\gamma.}{gamma.}
      
      If \sQuote{Specification} is 2, this value will not be changed
      by the MCMC and it is recommended (for easier
      interpretation of the results) to set \code{init$gamma} to zero
      (default behavior).
      
      If \sQuote{Specification} is 1 \code{init$gamma} should be
      approximately equal to the mean value of the residuals.
    }
    \item{sigma}{an~initial values of the basis
      standard deviation \eqn{\sigma}{sigma}.
      
      If \sQuote{Specification} is 2, this value will not be changed
      by the MCMC and it is recommended to set it
      approximately equal to the range of standardized data
      (let say 4 + 4)
      divided by the number of knots and
      multiplied by something like 2/3.

      If \sQuote{Specification} is 1
      this should be approximately equal to the range of the residuals
      divided by the number of knots \eqn{(2K+1)}{(2*K+1)} and
      multiplied again by something like 2/3.
      }  
    \item{intercept}{an~initial values of the
      intercept term \eqn{\alpha}{alpha}.
      
      If \sQuote{Specification} is 1 this value is not changed by the
      MCMC and the initial value is always changed to zero.}
    \item{scale}{an~initial value of the scale
      parameter \eqn{\tau}{tau}.

      If \sQuote{Specification} is 1 this value is not changed by the MCMC
      and the initial value is always changed to one.}
    \item{D}{initial value for the covariance matrix of random effects
      \eqn{D}{D}. Only its lower triangle must be given in a~vector,
      e.g. \code{c(d[1,1], d[2,1], d[3,1], d[2,2], d[3,2], d[3,3])}
      for a matrix \eqn{3 \times 3}{3 x 3}.}
    \item{b}{a~vector or matrix of the initial values of random effects
      \eqn{b_i,\;i=1,\dots,N}{b[i],\;i=1,..., N}
      for each cluster. The matrix should be of size
      \eqn{q\times N}{q x N}, where \eqn{q}{q} is the number of random
      effects. I.e. each column of the matrix contains the initial
      values for one cluster.}
    \item{y}{a~vector of length \eqn{\sum_{i=1}^N\,n_i}{sum[i=1][N] n[i]}
      with initial values of log-event-times.}
    \item{r}{a~vector of length \eqn{\sum_{i=1}^N\,n_i}{sum[i=1][N] n[i]}
      with initial
      component labels for each residual. All values must be between
      \eqn{-K}{-K} and \eqn{K.}{K.} See argument \code{init} of
      the function \code{\link{bayesHistogram}} for more details.}
    }
  }  %% end of item{init}  
  \item{mcmc.par}{a~list specifying how some of the G-spline parameters
    related to the distribution of the error term from \code{formula}
    are to be updated. See \code{\link{bayesBisurvreg}} for more
    details.

    In contrast to \code{\link{bayesBisurvreg}} function argument
    \code{mcmc.par$type.update.a} can also be equal to
    \code{"block"} in which case all \eqn{a}{a} coefficients are updated
    in 1 block using the Metropolis-Hastings algorithm. 
  }  %% end of item{mcmc.par}
  \item{prior2}{a~list specifying the prior distribution of the G-spline
    defining the distribution of the error term in the regression model
    given by \code{formula2} and \code{random2}. See \code{prior} argument of
    \code{\link{bayesHistogram}} function for more detail.
  }  %% end of item{prior2}
  \item{prior.b2}{prior specification for the parameters related to the
    random effects from \code{formula2} and \code{random2}. This should
    be a~list with the same structure as \code{prior.b}.
  }  %% end of item{prior.b2}  
  \item{prior.beta2}{prior specification for the regression parameters
    of time-to-event in the case of doubly censored data (related to
    \code{formula2} and \code{random2}).
    This should be a~list with the same structure as \code{prior.beta}.
  }  %% end of item{prior.beta2}
  \item{init2}{an~optional list with initial values for the MCMC related
    to the model given by \code{formula2} and \code{random2}.
    The list has the same structure as \code{init}.
  }  %% end of item{init2}
  \item{mcmc.par2}{a~list specifying how some of the G-spline parameters
    related to \code{formula2} are to be updated. The list has the same
    structure as \code{mcmc.par}.
  }  %% end of item{mcmc.par2}
  \item{store}{a~list of logical values specifying which chains that are
    not stored by default are to be stored. The list can have the
    following components.
    \describe{
      \item{a}{if \code{TRUE} then all the transformed mixture weights
	\eqn{a_{k},}{a[k],}
	\eqn{k=-K,\dots,K,}{k=-K,..., K,}
	related to the G-spline (error distribution) of \code{formula} are stored.}
      \item{a2}{if \code{TRUE} and there are doubly-censored data then
	all the transformed mixture weights
        \eqn{a_{k},}{a[k],}
	\eqn{k=-K,\dots,K,}{k=-K,..., K,}
	related to the G-spline (error distribution) of \code{formula2} are stored.}
      \item{y}{if \code{TRUE} then augmented log-event times for all
	observations related to the \code{formula} are stored.}
      \item{y2}{if \code{TRUE} then augmented log-event times for all
	observations related to \code{formula2} are stored.}
      \item{r}{if \code{TRUE} then labels of mixture components for
	residuals related to \code{formula} are stored.}
      \item{r2}{if \code{TRUE} then labels of mixture components for
	residuals related to \code{formula2} are stored.}
      \item{b}{if \code{TRUE} then the sampled values of the random
	effects related to \code{formula} and \code{random} are stored.}
      \item{b2}{if \code{TRUE} then the sampled values of the random
	effects related to \code{formula2} and \code{random2} are stored.}
    }  %% end of describe  
  }  %% end of item{store}
  \item{dir}{a string that specifies a directory where all sampled
    values are to be stored.
  }  %% end of item{dir}
}  %% end of arguments
\value{
  A list of class \code{bayessurvreg2} containing an information
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
    \item{mixmoment.sim}{columns labeled \code{k}, \code{Mean.1}, 
      \code{D.1.1}, where
      
      \bold{k} = number of mixture components that had probability
      numerically higher than zero;
      
      \bold{Mean.1} =
      \eqn{\mbox{E}(\varepsilon_{i,l})}{E(epsilon[i,l])};
            
      \bold{D.1.1} =
      \eqn{\mbox{var}(\varepsilon_{i,l})}{var(epsilon[i,l])};
            
      all related to the distribution of the error term from the
      model given by \code{formula}.
    }
    \item{mixmoment\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{mixmoment.sim}, however related to the model
      given by \code{formula2}.      
    }
    \item{mweight.sim}{sampled mixture weights
      \eqn{w_{k}}{w[k]} of mixture components that had
      probabilities numerically higher than zero. Related to the model
      given by \code{formula}.
    }
    \item{mweight\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{mweight.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{mmean.sim}{indeces \eqn{k,}{k,}
      \eqn{k \in\{-K, \dots, K\}}{k in {-K, ..., K}}
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
      \eqn{\varepsilon_{i,l}}{epsilon[i,l]})
      related to the model given by
      \code{formula}. This file together with \code{mixmoment.sim},
      \code{mweight.sim} and \code{mmean.sim} can be used to reconstruct
      the G-spline in each MCMC iteration.
      
      The file has columns labeled
      \code{gamma1},
      \code{sigma1},
      \code{delta1},
      \code{intercept1}, 
      \code{scale1},
      The meaning of the values in these columns is the following:
      
      \bold{gamma1} = the middle knot \eqn{\gamma}{gamma} 
      If \sQuote{Specification} is 2, this column usually contains zeros;
            
      \bold{sigma1} = basis standard deviation \eqn{\sigma}{sigma}
      of the G-spline. This column contains a~fixed value
      if \sQuote{Specification} is 2;
            
      \bold{delta1} = distance \eqn{delta}{delta} between the two knots of the G-spline.
      This column contains a~fixed value if \sQuote{Specification} is 2;
      
      \bold{intercept1} = the intercept term \eqn{\alpha}{alpha} of the G-spline.
      If \sQuote{Specification} is 1, this column usually contains zeros;
      
      \bold{scale1} = the scale parameter \eqn{\tau}{tau} of the G-spline.
      If \sQuote{Specification} is 1, this column usually contains ones;
    }
    \item{gspline\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{gspline.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{mlogweight.sim}{fully created only if \code{store$a = TRUE}. The
      file contains the transformed weights
      \eqn{a_{k},}{a[k],}
      \eqn{k=-K,\dots,K}{k=-K,..., K}
      of all mixture components, i.e. also of components that had numerically zero
      probabilities. This file is related to the error distribution of
      the model given by \code{formula}.
    }
    \item{mlogweight\_{}2.sim}{fully created only if \code{store$a2 =
      TRUE} and in the case of doubly-censored data, the same
      structure as \code{mlogweight.sim}, however related to the error
      distribution of the model given by \code{formula2}. 
    }
    \item{r.sim}{fully created only if \code{store$r = TRUE}. The file
      contains the labels of the mixture components into which the
      residuals are intrinsically assigned. Instead of indeces on the
      scale \eqn{\{-K,\dots, K\}}{{-K,..., K}}
      values from 1 to \eqn{(2\,K+1)}{(2*K+1)} are stored here. Function
      \code{\link{vecr2matr}} can be used to transform it back to
      indices from \eqn{-K}{-K} to \eqn{K}{K}.
    }
    \item{r\_{}2.sim}{fully created only if \code{store$r2 =
      TRUE} and in the case of doubly-censored data, the same
      structure as \code{r.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{lambda.sim}{one column labeled \code{lambda}. These are the
      values of the smoothing parameter\eqn{\lambda}{lambda}
      (hyperparameters of the prior distribution of the transformed
      mixture weights \eqn{a_{k}}{a[k]}). This file is
      related to the model given by \code{formula}.
    }
    \item{lambda\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{lambda.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{beta.sim}{sampled values of the regression parameters, both
      the fixed effects \eqn{\beta}{beta} and means of the random
      effects \eqn{\beta_b}{beta_b} (except the random intercept which
      has always the mean equal to zero).
      This file is related to the model given by \code{formula}.
      The columns are labeled according to the
      \code{colnames} of the design matrix.
    }
    \item{beta\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{beta.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{D.sim}{sampled values of the covariance matrix \eqn{D}{D} of
      the random effects. The file has \eqn{1 + 0.5\,q\,(q+1)}{1 +
	0.5*q*(q+1)} columns (\eqn{q}{q} is the dimension of the random
      effect vector \eqn{b_i}{b_i}). The first column labeled \code{det}
      contains the determinant of the sampled matrix, additional columns
      labeled \code{D.1.1}, \code{D.2.1}, \dots, \code{D.q.1}, \dots
      \code{D.q.q} contain the lower triangle of the sampled
      matrix. This file is related to the model specified by
      \code{formula} and \code{random}.
    }
    \item{D\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{D.sim}, however related to the model given by
      \code{formula2} and \code{random2}.
    }
    \item{b.sim}{fully created only if \code{store$b = TRUE}. It
      contains sampled values of random effects for all clusters in
      the data set. The file has \eqn{q\times N}{q*N} columns sorted as
      \eqn{b_{1,1},\dots,b_{1,q},\dots, b_{N,1},\dots,b_{N,q}}{%
           b[1,1],..., b[1,q],..., b[N,1],..., b[N,q]}. This file is
	 related to the model given by \code{formula} and \code{random}.
    }
    \item{b\_{}2.sim}{fully created only if \code{store$b2 =
      TRUE} and in the case of doubly-censored data, the same
      structure as \code{b.sim}, however related to the model
      given by \code{formula2} and \code{random2}. 
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
    \item{logposter.sim}{columns labeled \code{loglik}, \code{penalty},
      and \code{logprw}. This file is related to the model
      given by \code{formula}. The columns have the following meaning.
      
      \bold{loglik}
      \eqn{=}{=} \eqn{%
	- (\sum_{i=1}^N\,n_i)\,\Bigl\{\log(\sqrt{2\pi}) + \log(\sigma) \Bigr\}-
          0.5\sum_{i=1}^N\sum_{l=1}^{n_i}
	  \Bigl\{
	  (\sigma^2\,\tau^2)^{-1}\; (y_{i,l} - x_{i,l}'\beta - z_{i,l}'b_i -
	  \alpha - \tau\mu_{r_{i,l}})^2
          \Bigr\}
      }{%
	-(sum[i=1][N] n[i]) * (log(sqrt(2*pi)) + log(sigma))
          -0.5*sum[i=1][N] sum[l=1][n[i]](
          (sigma^2*tau^2)^(-1) * (y[i,l] - x[i,l]'beta - z[i,l]'b[i] - alpha - tau*mu[r[i,l]])^2)}
      
      where \eqn{y_{i,l}}{y[i,l]} denotes (augmented) \emph{(i,l)}th
      true log-event time.

      In other words, \code{loglik} is equal to the
      conditional log-density
      \deqn{\sum_{i=1}^N \sum_{l=1}^{n_i}\,\log\Bigl\{p\bigl(y_{i,l}\;\big|\;r_{i,l},\,\beta,\,b_i,\,\mbox{G-spline}\bigr)\Bigr\};}{%
      sum[i=1][N] sum[l=1][n[i]] log(p(y[i,l] | r[i,l], beta, b[i], G-spline));
      }
      
      \bold{penalty:}
      the penalty term
      \deqn{-\frac{1}{2}\sum_{k}\Bigl(\Delta\, a_k\Bigr)^2}{%
            -0.5*sum[k] (Delta a[k])^2}
      (not multiplied by \eqn{\lambda}{lambda});
      
      \bold{logprw} \eqn{=}{=}
      \eqn{-2\,(\sum_i n_i)\,\log\bigl\{\sum_{k}a_{k}\bigr\} +
	\sum_{k}N_{k}\,a_{k},}{%
	-2*(sum[i] n[i])*log(sum[k] exp(a[k])) +
	sum[k[1]] N[k]*a[k],}
      where \eqn{N_{k}}{N[k]} is the number of residuals
      assigned intrinsincally to the \eqn{k}{k}th
      mixture component.

      In other words, \code{logprw} is equal to the conditional
      log-density
      \eqn{\sum_{i=1}^N\sum_{l=1}^{n_i} \log\bigl\{p(r_{i,l}\;|\;\mbox{G-spline
	  weights})\bigr\}.}{%
      sum[i=1][N] sum[l=1][n[i]] log(p(r[i,l] | G-spline weights)).}
    }
    \item{logposter\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{logposter.sim}, however related to the model
      given by \code{formula2}. 
    }
  }
}  
\references{
  Gelfand, A. E., Sahu, S. K., and Carlin, B. P. (1995).
  Efficient parametrisations for normal linear mixed models.
  \emph{Biometrika,} \bold{82,} 479-488.

  \eqn{\mbox{Kom\'{a}rek, A.}}{Komarek, A.} (2006).
  \emph{Accelerated Failure Time Models for Multivariate
    Interval-Censored Data with Flexible Distributional Assumptions}.
  PhD. Thesis, Katholieke Universiteit Leuven, Faculteit Wetenschappen.

  \eqn{\mbox{Kom\'{a}rek}}{Komarek}, A. Lesaffre, E., and Legrand, C. (2006).
  EORTC data article.
  \emph{To be written.}    
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\examples{
## See the description of R commands for
## the model with EORTC data,
## analysis described in Komarek, Lesaffre and Legrand (2006),
## R commands will be soon available in the documentation
## directory of this package.
}
\keyword{survival}
\keyword{regression}
\keyword{multivariate}
