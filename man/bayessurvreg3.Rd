\name{bayessurvreg3}
\alias{bayessurvreg3}
\title{
  Cluster-specific accelerated failure time model for multivariate,
  possibly doubly-interval-censored data with flexibly specified random effects
  and/or error distribution.
}
\description{
  A function to estimate a regression model with possibly clustered
  (possibly right, left, interval or doubly-interval censored) data.
  In the case of doubly-interval censoring, different regression models
  can be specified for the onset and event times.

  A~univariate random effect (random intercept)
  with the distribution expressed as a~penalized normal
  mixture can be included in the model to adjust for clusters.

  The error density of the regression model is specified as a mixture of
  Bayesian G-splines (normal densities with equidistant means and
  constant variances). This function performs an MCMC sampling from the
  posterior distribution of unknown quantities.

  For details, see Komárek (2006)
  and Komárek and Lesaffre (2008).

  We explain first in more detail a model without doubly censoring.
  Let \eqn{T_{i,l},\; i=1,\dots, N,\; l=1,\dots, n_i}{%
           T[i,l], i=1,..., N, l=1,..., n[i]}
  be event times for \eqn{i}{i}th cluster and the units within that cluster
  The following regression model is assumed:
  \deqn{\log(T_{i,l}) = \beta'x_{i,l} + b_i + \varepsilon_{i,l},\quad i=1,\dots, N,\;l=1,\dots, n_i}{%
    log(T[i,l]) = beta'x[i,l] + b[i] + epsilon[i,l], i=1,..., N, l=1,..., n[i]}
  where \eqn{\beta}{beta} is unknown regression parameter vector,
  \eqn{x_{i,l}}{x[i,l]} is a vector of covariates.
  \eqn{b_i}{b[i]} is a cluster-specific random effect (random intercept).

  The random effects \eqn{b_i,\;i=1,\dots, N}{b[i], i=1,..., N}
  are assumed to be i.i.d. with a~univariate density \eqn{g_{b}(b)}{g[b](b)}.
  The error terms
  \eqn{\varepsilon_{i,l},\;i=1,\dots, N, l=1,\dots, n_i}{%
       epsilon[i,l], i=1,..., N, l=1,..., n[i]}
  are assumed to be i.i.d. with a~univariate density
  \eqn{g_{\varepsilon}(e)}{g[epsilon](e)}.
  
  Densities \eqn{g_{b}}{g[b]} and \eqn{g_{\varepsilon}}{g[epsilon]} are
  both expressed as
  a~mixture of Bayesian G-splines (normal densities with equidistant
  means and constant variances). We distinguish two,
  theoretically equivalent, specifications.

  In the following, the density for \eqn{\varepsilon}{epsilon}
  is explicitely described. The density for \eqn{b}{b} is obtained in
  an analogous manner.  

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
  Komárek and Lesaffre (2008) only
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
  details I refer to Komárek and Lesaffre (2008).
      
  If there are doubly-censored data the model of the same type as above
  can be specified for both the onset time and the time-to-event.

  In the case one wishes to link the random intercept of the onset model
  and the random intercept of the time-to-event model, there are the
  following possibilities.

  \bold{Bivariate normal distribution} \cr
  It is assumed that the pair of random intercepts from the onset and
  time-to-event part of the model are normally distributed with zero
  mean and an unknown covariance matrix \eqn{D}{D}.

  A priori, the inverse covariance matrix \eqn{D^{-1}}{D^(-1)} is
  addumed to follow a Wishart distribution.

  
  \bold{Unknown correlation between the basis G-splines} \cr
  Each pair of basis G-splines describing the distribution of the random
  intercept in the onset part and the time-to-event part of the model is
  assumed to be correlated with an unknown correlation coefficient
  \eqn{\varrho}{rho}. Note that this is just an experiment and is no
  more further supported.

  Prior distribution on \eqn{\varrho}{rho} is assumed to be
  uniform. In the MCMC, the Fisher Z transform of the \eqn{\varrho}{rho}
  given by
  \deqn{Z =
    -\frac{1}{2}\log\Bigl(\frac{1-\varrho}{1+\varrho}\Bigr)=\mbox{atanh}(\varrho)}{
  Z = -0.5*log((1-rho)/(1+rho)) = atanh(rho)}
  is sampled. Its prior is derived from the uniform prior
  \eqn{\mbox{Unif}(-1,\;1)}{Unif(-1, 1)} put on \eqn{\varrho.}{rho.}

  The Fisher Z transform is updated using the Metropolis-Hastings
  alhorithm. The proposal distribution is given either by a normal
  approximation obtained using the Taylor expansion of the full
  conditional distribution or by a Langevin proposal (see Robert and
  Casella, 2004, p. 318).  
}
\usage{
bayessurvreg3(formula, random, formula2, random2,
   data = parent.frame(),
   na.action = na.fail, onlyX = FALSE,
   nsimul = list(niter = 10, nthin = 1, nburn = 0, nwrite = 10),   
   prior, prior.beta, prior.b, init = list(iter = 0),
   mcmc.par = list(type.update.a = "slice", k.overrelax.a = 1,
                   k.overrelax.sigma = 1, k.overrelax.scale = 1,
                   type.update.a.b = "slice", k.overrelax.a.b = 1,
                   k.overrelax.sigma.b = 1, k.overrelax.scale.b = 1),
   prior2, prior.beta2, prior.b2, init2,
   mcmc.par2 = list(type.update.a = "slice", k.overrelax.a = 1,
                    k.overrelax.sigma = 1, k.overrelax.scale = 1,
                   type.update.a.b = "slice", k.overrelax.a.b = 1,
                   k.overrelax.sigma.b = 1, k.overrelax.scale.b = 1),
   priorinit.Nb,
   rho = list(type.update = "fixed.zero", init=0, sigmaL=0.1),
   store = list(a = FALSE, a2 = FALSE, y = FALSE, y2 = FALSE,
                r = FALSE, r2 = FALSE, b = FALSE, b2 = FALSE,
                a.b = FALSE, a.b2 = FALSE, r.b = FALSE, r.b2 = FALSE), 
   dir = getwd())
}
\arguments{
  \item{formula}{model formula for the regression. In the case of
    doubly-censored data, this is the model formula for the onset
    time.

    The left-hand side of the \code{formula} must be an~object created
    using \code{\link[survival]{Surv}}.

    Intercept is implicitely included in the model by the
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
  \item{random}{formula for the `random' part of the model.
    In the case of doubly-censored data, this is the \code{random} formula for
    the onset time. With this version of the function only
    \tabular{c}{
      \code{random = ~1}
    }     
    is allowed. If omitted, no random part is included in the model. 
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
  \item{prior.b}{a~list specifying the prior distribution of the
    G-spline defining the distribution of the random intercept in the
    regression model given by \code{formula} and \code{random}. See
    \code{prior} argument of \code{\link{bayesHistogram}} function for
    more detail. In this list
    also \sQuote{Specification} as described above is specified.

    It is ignored if the argument \code{priorinit.Nb} is given.

    The item \code{prior.b$neighbor.system} can only be equal to
    \code{uniCAR} here.
  }  %% end of item{prior.b}  
  \item{prior.beta}{prior specification for the regression parameters,
    in the case of doubly-censored data for the regression parameters of
    the onset time, i.e. it is related to \code{formula} and
    \code{random}.
    
    This should be a~list with the following components:
    \describe{
    \item{mean.prior}{a~vector specifying a~prior mean for each
      \code{beta} parameter in the model.}
    \item{var.prior}{a~vector specifying a~prior variance for each
      \code{beta} parameter.}
    }
    It is recommended to run the function
    bayessurvreg3 first with its argument \code{onlyX} set to \code{TRUE}
    to find out how the betas are sorted. They must correspond to a
    design matrix X taken from \code{formula}.
  }  %% end of item{prior.beta}
  \item{init}{an~optional list with initial values for the MCMC related
    to the model given by \code{formula} and \code{random}. The list can have the following components:
    \describe{
    \item{iter}{the number of the iteration to which the initial values
      correspond, usually zero.}
    \item{beta}{a~vector of initial values for the regression
      parameters.
      It must be sorted in the same way as are the columns
      in the design matrix. Use \code{onlyX=TRUE} if you do not know how
      the columns in the design matrix are created.}
    \item{a}{a~vector of length \eqn{2K+1}{2*K+1} with the initial
      values of transformed mixture weights for the G-spline defining
      the distribution of the error term \eqn{\varepsilon}{epsilon}.}
    \item{lambda}{initial values for the Markov random fields precision
      parameter for the G-spline defining
      the distribution of the error term \eqn{\varepsilon}{epsilon}. 
    }
    \item{gamma}{an~initial values for the middle
      knot \eqn{\gamma}{gamma} for the G-spline defining
      the distribution of the error term \eqn{\varepsilon}{epsilon}.
      
      If \sQuote{Specification} is 2, this value will not be changed
      by the MCMC and it is recommended (for easier
      interpretation of the results) to set \code{init$gamma} to zero
      (default behavior).
      
      If \sQuote{Specification} is 1 \code{init$gamma} should be
      approximately equal to the mean value of the residuals.
    }
    \item{sigma}{an~initial values of the basis
      standard deviation \eqn{\sigma}{sigma} for the G-spline defining
      the distribution of the error term \eqn{\varepsilon}{epsilon}.
      
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
      intercept term \eqn{\alpha}{alpha} for the G-spline defining
      the distribution of the error term \eqn{\varepsilon}{epsilon}.
      
      If \sQuote{Specification} is 1 this value is not changed by the
      MCMC and the initial value is always changed to zero.}
    \item{scale}{an~initial value of the scale
      parameter \eqn{\tau}{tau} for the G-spline defining
      the distribution of the error term \eqn{\varepsilon}{epsilon}.

      If \sQuote{Specification} is 1 this value is not changed by the MCMC
      and the initial value is always changed to one.}
    \item{a.b}{a~vector of length \eqn{2K+1}{2*K+1} with the initial
      values of transformed mixture weights for the G-spline defining
      the distribution of the random intercept \eqn{b}{b}.}
    \item{lambda.b}{initial values for the Markov random fields precision
      parameter for the G-spline defining
      the distribution of the random intercept \eqn{b}{b}. 
    }
    \item{gamma.b}{an~initial values for the middle
      knot \eqn{\gamma}{gamma} for the G-spline defining
      the distribution of the random intercept \eqn{b}{b}.

      Due to identifiability reasons, this value is always changed to
      zero and is for neither \sQuote{Specification} updated by the
      MCMC.
    }
    \item{sigma.b}{an~initial values of the basis
      standard deviation \eqn{\sigma}{sigma} for the G-spline defining
      the distribution of the random intercept \eqn{b}{b}.
      
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
    \item{intercept.b}{an~initial values of the
      intercept term \eqn{\alpha}{alpha} for the G-spline defining
      the distribution of the random intercept \eqn{b}{b}.

      Due to identifiability reasons, this value is always changed to
      zero and is for neither \sQuote{Specification} updated by the
      MCMC.
      }
    \item{scale.b}{an~initial value of the scale
      parameter \eqn{\tau}{tau} for the G-spline defining
      the distribution of the random intercept \eqn{b}{b}.

      If \sQuote{Specification} is 1 this value is not changed by the MCMC
      and the initial value is always changed to one.}
    \item{b}{a vector of length \eqn{N}{N} of the initial values of random effects
      \eqn{b_i,\;i=1,\dots,N}{b[i],\;i=1,..., N}
      for each cluster.}
    \item{y}{a vector of length \eqn{\sum_{i=1}^N\,n_i}{sum[i=1][N] n[i]}
      with initial values of log-event-times.}
    \item{r}{a vector of length \eqn{\sum_{i=1}^N\,n_i}{sum[i=1][N] n[i]}
      with initial
      component labels for each residual. All values must be between
      \eqn{-K}{-K} and \eqn{K.}{K.} See argument \code{init} of
      the function \code{\link{bayesHistogram}} for more details.}
    \item{r.b}{a~vector of length \eqn{N}{N}
      with initial
      component labels for each random intercept. All values must be between
      \eqn{-K}{-K} and \eqn{K.}{K.} See argument \code{init} of
      the function \code{\link{bayesHistogram}} for more details.}    
    } 
  }  %% end of item{init}  
  \item{mcmc.par}{a list specifying how some of the G-spline parameters
    related to the distribution of the error term and of the random
    intercept from \code{formula} and \code{random}
    are to be updated. See \code{\link{bayesBisurvreg}} for more
    details.

    Compared to the mcmc.par argument of the function
    \code{\link{bayesBisurvreg}} additional components related to the
    G-spline for the random intercept can be present, namely

    \tabular{l}{
      \code{type.update.a.b} \cr
      \code{k.overrelax.a.b} \cr
      \code{k.overrelax.sigma.b} \cr
      \code{k.overrelax.scale.b}
    }

    In contrast to \code{\link{bayesBisurvreg}} function arguments
    \code{mcmc.par$type.update.a} and \code{mcmc.par$type.update.a.b} can also be equal to
    \code{"block"} in which case all \eqn{a}{a} coefficients are updated
    in 1 block using the Metropolis-Hastings algorithm. 
  }  %% end of item{mcmc.par}
  \item{prior2}{a list specifying the prior distribution of the G-spline
    defining the distribution of the error term in the regression model
    given by \code{formula2} and \code{random2}. See \code{prior} argument of
    \code{\link{bayesHistogram}} function for more detail.
  }  %% end of item{prior2}
  \item{prior.b2}{prior specification for the parameters related to the
    random effects from \code{formula2} and \code{random2}. This should
    be a~list with the same structure as \code{prior.b}.

    It is ignored if the argument \code{priorinit.Nb} is given.    
  }  %% end of item{prior.b2}    
  \item{prior.beta2}{prior specification for the regression parameters
    of time-to-event in the case of doubly censored data (related to
    \code{formula2} and \code{random2}).
    This should be a~list with the same structure as \code{prior.beta}.
  }  %% end of item{prior.beta2}
  \item{init2}{an optional list with initial values for the MCMC related
    to the model given by \code{formula2} and \code{random2}.
    The list has the same structure as \code{init}.
  }  %% end of item{init2}
  \item{mcmc.par2}{a list specifying how some of the G-spline parameters
    related to \code{formula2} and \code{random2} are to be updated.
    The list has the same structure as \code{mcmc.par}.
  }  %% end of item{mcmc.par2}
  \item{priorinit.Nb}{a list specifying the prior of the random intercepts
    in the case of the AFT model with doubly-interval-censored data and
    onset, time-to-event random intercepts following bivariate normal
    distribution.

    The list should have the following components.
    \describe{
      \item{init.D}{initial value for the covariance matrix of the onset
	random intercept and time-to-event random intercept.

	It can be specified either as a vector of length 3 giving the
	lower triangle of the matrix or as a matrix 2 x 2.
      }	
      \item{df.Di}{degrees of freedom \eqn{\nu}{nu} for the Wishart prior of the
	matrix \eqn{D^{-1}}{D^(-1)}.

	Note that it must be higher than 1.
      }
      \item{scale.Di}{scale matrix \eqn{S}{S} for the Wishart prior of the
	matrix \eqn{D^{-1}}{D^(-1)}.

	It can be specified either as a vector of length 3 giving the
	lower triangle of the matrix or as a matrix 2 x 2.

	Note that a priori
	\deqn{\mbox{E}(D^{-1}.) = \nu S}{E(D^{-1}) = nu*S.}
      }      
    }  
  }    
  \item{rho}{a list specifying possible correlation between the onset
    random intercept and the time-to-event random intercept in the
    experimental version of the model. If not given correlation is fixed
    to \eqn{0}{0}.

    It is ignored if the argument \code{priorinit.Nb} is given.
    Ordinary users should not care about this argument.

    The list can have the following components.
    \describe{
      \item{type.update}{character specifying how the Fisher Z transform
	of the correlation coefficient is updated. Possible values are:
	
	\code{"fixed.zero"}:
	correlation coefficient is fixed to \eqn{0}{0} and it is not updated.

	\code{"normal.around.mode"}:
        at each iteration of MCMC, 1 Newton-Raphson step from the
	current point \eqn{Z}{Z} of the full conditional distribution is
	performed, normal approximation is formed by Taylor expansion
	and new point \eqn{Z}{Z} is sampled from that normal
	approximation.

	Note that this proposal does not work too well if the current
	point \eqn{Z}{Z} lies in the area of low posterior mass. The
	reason is that even 1 Newton-Raphson step usually leads to the
	area of high posterior probability mass and the proposal is
	``too ambisious''.
	
	\code{"langevin"}.
	 at each iteration of MCMC, new point \eqn{Z}{Z} is sampled
	 using the Langevin algorithm. A scale parameter (see below)
	 must cerefully be chosen for this algorithm to ensure that the
	 acceptance rate is about 50--60\% (Robert, Casella, 2004, p. 319).
      }		
    }  
  }  
  \item{store}{a list of logical values specifying which chains that are
    not stored by default are to be stored. The list can have the
    following components.
    \describe{
      \item{a}{if \code{TRUE} then all the transformed mixture weights
	\eqn{a_{k},}{a[k],}
	\eqn{k=-K,\dots,K,}{k=-K,..., K,}
	related to the G-spline defining the error distribution of \code{formula}
	are stored.}
      \item{a.b}{if \code{TRUE} then all the transformed mixture weights
	\eqn{a_{k},}{a[k],}
	\eqn{k=-K,\dots,K,}{k=-K,..., K,}
	related to the G-spline defining the distribution of the random
	intercept from \code{formula} and \code{random} are stored.}
      \item{a2}{if \code{TRUE} and there are doubly-censored data then
	all the transformed mixture weights
        \eqn{a_{k},}{a[k],}
	\eqn{k=-K,\dots,K,}{k=-K,..., K,}
	related to the G-spline defining the error distribution of
	\code{formula2} are stored.}
      \item{a.b2}{if \code{TRUE} then all the transformed mixture weights
	\eqn{a_{k},}{a[k],}
	\eqn{k=-K,\dots,K,}{k=-K,..., K,}
	related to the G-spline defining the distribution of the random
	intercept from \code{formula2} and \code{random2} are stored.}      
      \item{y}{if \code{TRUE} then augmented log-event times for all
	observations related to the \code{formula} are stored.}
      \item{y2}{if \code{TRUE} then augmented log-event times for all
	observations related to \code{formula2} are stored.}
      \item{r}{if \code{TRUE} then labels of mixture components for
	residuals related to \code{formula} are stored.}      
      \item{r.b}{if \code{TRUE} then labels of mixture components for
	random intercepts related to \code{formula} and \code{random}
	are stored.}      
      \item{r2}{if \code{TRUE} then labels of mixture components for
	residuals related to \code{formula2} are stored.}
      \item{r.b2}{if \code{TRUE} then labels of mixture components for
	random intercepts related to \code{formula2} and \code{random2}
	are stored.}            
      \item{b}{if \code{TRUE} then the sampled values of the random
	interceptss related to \code{formula} and \code{random} are stored.}
      \item{b2}{if \code{TRUE} then the sampled values of the random
	interceptss related to \code{formula2} and \code{random2} are stored.}
    }  %% end of describe  
  }  %% end of item{store}
  \item{dir}{a string that specifies a directory where all sampled
    values are to be stored.
  }  %% end of item{dir}  
}  %% end of arguments
\value{
  A list of class \code{bayessurvreg3} containing an information
  concerning the initial values and prior choices.
}
\section{Files created}{%%%AAA  
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
  \describe{%%%BBB
    \item{iteration.sim}{one column labeled \code{iteration} with
      indeces of MCMC iterations to which the stored sampled values
      correspond.
    }
    \item{mixmoment.sim}{this file is related to the density of the
      error term from the model given by \code{formula}.

      Columns labeled \code{k}, \code{Mean.1}, 
      \code{D.1.1}, where
      
      \bold{k} = number of mixture components that had probability
      numerically higher than zero;
      
      \bold{Mean.1} =
      \eqn{\mbox{E}(\varepsilon_{i,l})}{E(epsilon[i,l])};
            
      \bold{D.1.1} =
      \eqn{\mbox{var}(\varepsilon_{i,l})}{var(epsilon[i,l])}.
    }
    \item{mixmoment\_{}b.sim}{this file is related to the density of the
      random intercept from the model given by \code{formula} and
      \code{random}.

      The same structure as \code{mixmoment.sim}.
    }
    \item{mixmoment\_{}2.sim}{in the case of doubly-censored data. This
      file is related to the density of the error term from the model
      given by \code{formula2}.
      
      The same structure as \code{mixmoment.sim}.
    }
    \item{mixmoment\_{}b2.sim}{in the case of doubly-censored data. This
      file is related to the density of the random intercept from the model
      given by \code{formula2} and \code{random2}.
      
      The same structure as \code{mixmoment.sim}.
    }    
    \item{mweight.sim}{this file is related to the density of the
      error term from the model given by \code{formula}.

      Sampled mixture weights \eqn{w_{k}}{w[k]} of mixture components that had
      probabilities numerically higher than zero. 
    }
    \item{mweight\_{}b.sim}{this file is related to the density of the
      random intercept from the model given by \code{formula} and
      \code{random}.

      The same structure as \code{mweight.sim}.
    }
    \item{mweight\_{}2.sim}{in the case of doubly-censored data. This
      file is related to the density of the error term from the model
      given by \code{formula2}.
      
      The same structure as \code{mweight.sim}.
    }
    \item{mweight\_{}b2.sim}{in the case of doubly-censored data. This
      file is related to the density of the random intercept from the model
      given by \code{formula2} and \code{random2}.
      
      The same structure as \code{mweight.sim}.
    }    
    \item{mmean.sim}{this file is related to the density of the
      error term from the model given by \code{formula}.

      Indeces \eqn{k,}{k,}
      \eqn{k \in\{-K, \dots, K\}}{k in {-K, ..., K}}
      of mixture components that had probabilities numerically higher
      than zero. It corresponds to the weights in
      \code{mweight.sim}. 
    }
    \item{mmean\_{}b.sim}{this file is related to the density of the
      random intercept from the model given by \code{formula} and
      \code{random}.

      The same structure as \code{mmean.sim}.
    }
    \item{mmean\_{}2.sim}{in the case of doubly-censored data. This
      file is related to the density of the error term from the model
      given by \code{formula2}.
      
      The same structure as \code{mmean.sim}.
    }
    \item{mmean\_{}b2.sim}{in the case of doubly-censored data. This
      file is related to the density of the random intercept from the model
      given by \code{formula2} and \code{random2}.
      
      The same structure as \code{mmean.sim}.
    }    
    \item{gspline.sim}{this file is related to the density of the
      error term from the model given by \code{formula}.

      Characteristics of the sampled G-spline.
      This file together with \code{mixmoment.sim},
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
    \item{gspline\_{}b.sim}{this file is related to the density of the
      random intercept from the model given by \code{formula} and
      \code{random}.

      The same structure as \code{gspline.sim}.
    }
    \item{gspline\_{}2.sim}{in the case of doubly-censored data. This
      file is related to the density of the error term from the model
      given by \code{formula2}.
      
      The same structure as \code{gspline.sim}.
    }
    \item{gspline\_{}b2.sim}{in the case of doubly-censored data. This
      file is related to the density of the random intercept from the model
      given by \code{formula2} and \code{random2}.
      
      The same structure as \code{gspline.sim}.
    }    
    \item{mlogweight.sim}{this file is related to the density of the
      error term from the model given by \code{formula}.

      Fully created only if \code{store$a = TRUE}. The
      file contains the transformed weights
      \eqn{a_{k},}{a[k],}
      \eqn{k=-K,\dots,K}{k=-K,..., K}
      of all mixture components, i.e. also of components that had numerically zero
      probabilities. 
    }
    \item{mlogweight\_{}b.sim}{this file is related to the density of the
      random intercept from the model given by \code{formula} and
      \code{random}.

      Fully created only if \code{store$a.b = TRUE}.

      The same structure as \code{mlogweight.sim}.
    }
    \item{mlogweight\_{}2.sim}{in the case of doubly-censored data. This
      file is related to the density of the error term from the model
      given by \code{formula2}.

      Fully created only if \code{store$a2 = TRUE}.
      
      The same structure as \code{mlogweight.sim}.
    }
    \item{mlogweight\_{}b2.sim}{in the case of doubly-censored data. This
      file is related to the density of the random intercept from the model
      given by \code{formula2} and \code{random2}.

      Fully created only if \code{store$a.b2 = TRUE}.
      
      The same structure as \code{mlogweight.sim}.
    }    
    \item{r.sim}{this file is related to the density of the
      error term from the model given by \code{formula}.

      Fully created only if \code{store$r = TRUE}. The file
      contains the labels of the mixture components into which the
      residuals are intrinsically assigned. Instead of indeces on the
      scale \eqn{\{-K,\dots, K\}}{{-K,..., K}}
      values from 1 to \eqn{(2\,K+1)}{(2*K+1)} are stored here. Function
      \code{\link{vecr2matr}} can be used to transform it back to
      indices from \eqn{-K}{-K} to \eqn{K}{K}.
    }
    \item{r\_{}b.sim}{this file is related to the density of the
      random intercept from the model given by \code{formula} and
      \code{random}.

      Fully created only if \code{store$r.b = TRUE}.

      The same structure as \code{r.sim}.
    }
    \item{r\_{}2.sim}{in the case of doubly-censored data. This
      file is related to the density of the error term from the model
      given by \code{formula2}.

      Fully created only if \code{store$r2 = TRUE}.
      
      The same structure as \code{r.sim}.
    }
    \item{r\_{}b2.sim}{in the case of doubly-censored data. This
      file is related to the density of the random intercept from the model
      given by \code{formula2} and \code{random2}.

      Fully created only if \code{store$r.b2 = TRUE}.
      
      The same structure as \code{r.sim}.
    }    
    \item{lambda.sim}{this file is related to the density of the
      error term from the model given by \code{formula}.

      One column labeled \code{lambda}. These are the
      values of the smoothing parameter\eqn{\lambda}{lambda}
      (hyperparameters of the prior distribution of the transformed
      mixture weights \eqn{a_{k}}{a[k]}). 
    }
    \item{lambda\_{}b.sim}{this file is related to the density of the
      random intercept from the model given by \code{formula} and
      \code{random}.

      The same structure as \code{lambda.sim}.
    }
    \item{lambda\_{}2.sim}{in the case of doubly-censored data. This
      file is related to the density of the error term from the model
      given by \code{formula2}.
      
      The same structure as \code{lambda.sim}.
    }
    \item{lambda\_{}b2.sim}{in the case of doubly-censored data. This
      file is related to the density of the random intercept from the model
      given by \code{formula2} and \code{random2}.
      
      The same structure as \code{lambda.sim}.
    }    
    \item{beta.sim}{this file is related to the model given by
      \code{formula}.

      Sampled values of the regression parameters
      \eqn{\beta}{beta}.

      The columns are labeled according to the
      \code{colnames} of the design matrix.
    }
    \item{beta\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{beta.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{b.sim}{this file is related to the model given by
      \code{formula} and \code{random}.

      Fully created only if \code{store$b = TRUE}. It
      contains sampled values of random intercepts for all clusters in
      the data set. The file has \eqn{N}{N} columns.
    }
    \item{b\_{}2.sim}{fully created only if \code{store$b2 =
      TRUE} and in the case of doubly-censored data, the same
      structure as \code{b.sim}, however related to the model
      given by \code{formula2} and \code{random2}. 
    }  
    \item{Y.sim}{this file is related to the model given by \code{formula}.

      Fully created only if \code{store$y = TRUE}. It
      contains sampled (augmented) log-event times for all observations
      in the data set.
    }
    \item{Y\_{}2.sim}{fully created only if \code{store$y2 =
      TRUE} and in the case of doubly-censored data, the same
      structure as \code{Y.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{logposter.sim}{
      This file is related to the residuals of the model
      given by \code{formula}. 

      Columns labeled \code{loglik}, \code{penalty}, and \code{logprw}. 
      The columns have the following meaning.
      
      \bold{loglik}
      \eqn{=}{=} \eqn{%
	- (\sum_{i=1}^N\,n_i)\,\Bigl\{\log(\sqrt{2\pi}) + \log(\sigma) \Bigr\}-
          0.5\sum_{i=1}^N\sum_{l=1}^{n_i}
	  \Bigl\{
	  (\sigma^2\,\tau^2)^{-1}\; (y_{i,l} - x_{i,l}'\beta - b_i -
	  \alpha - \tau\mu_{r_{i,l}})^2
          \Bigr\}
      }{%
	-(sum[i=1][N] n[i]) * (log(sqrt(2*pi)) + log(sigma))
          -0.5*sum[i=1][N] sum[l=1][n[i]](
          (sigma^2*tau^2)^(-1) * (y[i,l] - x[i,l]'beta - b[i] - alpha - tau*mu[r[i,l]])^2)}
      
      where \eqn{y_{i,l}}{y[i,l]} denotes (augmented) \emph{(i,l)}th
      true log-event time.

      In other words, \code{loglik} is equal to the
      conditional log-density
      \deqn{\sum_{i=1}^N \sum_{l=1}^{n_i}\,\log\Bigl\{p\bigl(y_{i,l}\;\big|\;r_{i,l},\,\beta,\,b_i,\,\mbox{error-G-spline}\bigr)\Bigr\};}{%
      sum[i=1][N] sum[l=1][n[i]] log(p(y[i,l] | r[i,l], beta, b[i], error-G-spline));
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
      \deqn{\sum_{i=1}^N\sum_{l=1}^{n_i} \log\bigl\{p(r_{i,l}\;|\;\mbox{error-G-spline
	  weights})\bigr\}.}{%
	sum[i=1][N] sum[l=1][n[i]] log(p(r[i,l] | error-G-spline weights)).}
    }
    \item{logposter\_{}b.sim}{This file is related to the random
      intercepts from the model given by \code{formula} and
      \code{random}. 

      Columns labeled \code{loglik}, \code{penalty}, and \code{logprw}. 
      The columns have the following meaning.
      
      \bold{loglik}
      \eqn{=}{=} \eqn{%
	- N\,\Bigl\{\log(\sqrt{2\pi}) + \log(\sigma) \Bigr\}-
          0.5\sum_{i=1}^N
	  \Bigl\{
	  (\sigma^2\,\tau^2)^{-1}\; (b_i - \alpha - \tau\mu_{r_{i}})^2
          \Bigr\}
      }{%
	-N * (log(sqrt(2*pi)) + log(sigma))
          -0.5*sum[i=1][N](
          (sigma^2*tau^2)^(-1) * (b[i] - alpha - tau*mu[r[i]])^2)}
      
      where \eqn{b_{i}}{b[i]} denotes (augmented) \emph{i}th
      random intercept.

      In other words, \code{loglik} is equal to the
      conditional log-density
      \deqn{\sum_{i=1}^N \,\log\Bigl\{p\bigl(b_{i}\;\big|\;r_{i},\,\mbox{b-G-spline}\bigr)\Bigr\};}{%
      sum[i=1][N] log(p(b[i] | r[i], b-G-spline));
      }

      The columns \code{penalty} and \code{logprw} have the analogous
      meaning as in the case of logposter.sim file.
    }
    \item{logposter\_{}2.sim}{in the case of doubly-censored data, the same
      structure as \code{logposter.sim}, however related to the model
      given by \code{formula2}. 
    }
    \item{logposter\_{}b2.sim}{in the case of doubly-censored data, the same
      structure as \code{logposter\_{}b.sim}, however related to the model
      given by \code{formula2} and \code{random2}.
    }
  }
}
\references{
  Komárek, A. (2006).
  \emph{Accelerated Failure Time Models for Multivariate
    Interval-Censored Data with Flexible Distributional Assumptions}.
  PhD. Thesis, Katholieke Universiteit Leuven, Faculteit Wetenschappen.
  
  Komárek, A. and Lesaffre, E. (2008).
  Bayesian accelerated failure time model with multivariate doubly-interval-censored data
  and flexible distributional assumptions.
  \emph{Journal of the American Statistical Association}, \bold{103}, 523--533.

  Robert C. P. and Casella, G. (2004).
  \emph{Monte Carlo Statistical Methods, Second Edition.}
  New York: Springer Science+Business Media.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\examples{
## See the description of R commands for
## the cluster specific AFT model
## with the Signal Tandmobiel data,
## analysis described in Komarek and Lesaffre (2007).
##
## R commands available in the documentation
## directory of this package
## - see ex-tandmobCS.R and
##   http://www.karlin.mff.cuni.cz/~komarek/software/bayesSurv/ex-tandmobCS.pdf
##
}
\keyword{survival}
\keyword{regression}
\keyword{multivariate}
