\name{bayesHistogram}
\alias{bayesHistogram}
\title{Smoothing of a uni- or bivariate histogram using Bayesian
  G-splines}
\description{A function to estimate a density of a uni- or bivariate
  (possibly censored) sample. The density is specified as a mixture of Bayesian G-splines
  (normal densities with equidistant means and equal variances). This function performs an
  MCMC sampling from the posterior distribution of unknown quantities in
  the density specification. Other method functions are available to
  visualize resulting density estimate.

  This function served as a basis for further developed
  \code{\link{bayesBisurvreg}}, \code{\link{bayessurvreg2}} and
  \code{\link{bayessurvreg3}} functions. However, in contrast to these
  functions, \code{bayesHistogram} does not allow for doubly censoring.
  
  \bold{Bivariate case:}

  Let \eqn{Y_{i,l},\; i=1,\dots,N,\; l=1,2}{Y[i,l], i=1,..., N, l=1,2} be
  observations for the \eqn{i}{i}th cluster and the first and the second
  unit (dimension). The bivariate observations
  \eqn{Y_i=(Y_{i,1},\,Y_{i,2})',\;i=1,\dots,N}{Y[i] = (Y[i,1], Y[i,2])', i=1,..., N}
  are assumed to be i.i.d. with a~bivariate density
  \eqn{g_{y}(y_1,\,y_2)}{g[y](y[1], y[2])}. This density is expressed as
  a~mixture of Bayesian G-splines (normal densities with equidistant
  means and constant variance matrices). We distinguish two,
  theoretically equivalent, specifications.

  \describe{%%% describe
    \item{Specification 1}{
      \deqn{(Y_1,\,Y_2)' \sim
	\sum_{j_1=-K_1}^{K_1}\sum_{j_2=-K_2}^{K_2} w_{j_1,j_2} N_2(\mu_{(j_1,j_2)},\,\mbox{diag}(\sigma_1^2,\,\sigma_2^2))}{%
        (Y[1],\,Y[2])' is distributed as
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
    }
    \item{Specification 2}{
      \deqn{(Y_1,\,Y_2)' \sim (\alpha_1,\,\alpha_2)'+ \bold{S}\,(Y_1,\,Y_2)'}{%
	(Y[1],\,Y[2])' is distributed as (alpha[1], alpha[2])' + S (V[1], V[2])'}
      where \eqn{(\alpha_1,\,\alpha_2)'}{(alpha[1], alpha[2])'} is an
      \bold{unknown} intercept term and
      \eqn{\bold{S} \mbox{ is a diagonal matrix with } \tau_1 \mbox{ and }\tau_2 \mbox{ on a diagonal,}}{%
	S is a diagonal matrix with tau[1] and tau[2] on a diagonal,}
      i.e. \eqn{\tau_1,\,\tau_2}{tau[1], tau[2]} are \bold{unknown} scale
      parameters. \eqn{(V_1,\,V_2)'}{(V[1], V[2])'} is then
      standardized observational vector which is distributed according
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
  }  %% end of describe

  \bold{Univariate case:}

  It is a~direct simplification of the bivariate case.
}
\usage{
bayesHistogram(y1, y2,
   nsimul = list(niter = 10, nthin = 1, nburn = 0, nwrite = 10),
   prior, init = list(iter = 0),
   mcmc.par = list(type.update.a = "slice", k.overrelax.a = 1,
                   k.overrelax.sigma = 1, k.overrelax.scale = 1),
   store = list(a = FALSE, y = FALSE, r = FALSE),
   dir = getwd())
}
\arguments{
  \item{y1}{response for the first dimension in the form of a survival
    object created using \code{\link[survival]{Surv}}.}
  \item{y2}{response for the second dimension in the form of a survival
    object created using \code{\link[survival]{Surv}}. If the response
    is one-dimensional this item is missing.}
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
  \item{prior}{a list that identifies prior hyperparameters and prior
    choices. See the paper \eqn{\mbox{Kom{\'a}rek}}{Kom&#225rek} and
    Lesaffre (2007) and the PhD. thesis
    \eqn{\mbox{Kom\'{a}rek}}{Kom&#225rek} (2006) for more details.
    
    Some prior parameters can be guessed by the function itself. If you
    want to do so, set such parameters to \code{NULL}. Set to
    \code{NULL} also the parameters that are not needed in your model.
    \describe{
      \item{specification}{a~number giving which specification of the
	model is used. It can be one of the following numbers:
	\describe{
	  \item{1}{with this specification positions of the middle knots
	    \eqn{\gamma_1,\dots,\gamma_q,}{gamma[1],...,gamma[q],}
	    where \eqn{q}{q} is dimension of the G-spline and basis
	    standard deviations
	    \eqn{\sigma_{0,1},\dots,\sigma_{0,q}}{sigma[0,1],...,sigma[0,q]}
	    are estimated. At the same time the G-spline intercepts
	    \eqn{\alpha_1,\dots,\alpha_q}{alpha[1],...,alpha[q]} and the
	    G-spline scale parameters
	    \eqn{s_{1},\dots,s_{q}}{s[1],...,s[q]}
	    are assumed to be fixed (usually, intercepts to zero and
	    scales to 1). The user can specified the fixed quantities in
	    the \code{init} parameter of this function
	  }  
	  \item{2}{with this specification, G-spline intercepts
	    \eqn{\alpha_1,\dots,\alpha_q}{alpha[1],...,alpha[q]} and the
	    G-spline scale parameters
	    \eqn{s_{1},\dots,s_{q}}{s[1],...,s[q]} are estimated
	    at the same time positions of the middle knots
	    \eqn{\gamma_1,\dots,\gamma_q}{gamma[1],...,gamma[q]}
            and basis standard deviations
	    \eqn{\sigma_{0,1},\dots,\sigma_{0,q}}{sigma[0,1],...,sigma[0,q]}
	    are assumed to be fixed (usually, middle knots to zero ans
	    basis standard deviations to some smaller number like 0.2)
            The user can specified the fixed quantities in the
	    \code{init} parameter of this function	    
	  }	  
        }  %% end of describe
      }  %% end of item{specification}  
      \item{K}{specification of the number of knots in each dimension, i.e. \code{K} 
        is a vector of length equal to the dimension of the data \eqn{q}{q}
        and \eqn{K_j,}{K[j],} \eqn{j=1,\dots,q}{j=1,...,q} determines that 
        the subscript \eqn{k_j}{k[j]} of the knots runs over
        \eqn{-K_j,\dots,0,\dots,K_j}{-K[j],...,0,...,K[j]}. A value
        \eqn{K_j=0}{K[j]=0} is valid as well. There are only some restriction
        on the minimal value of \eqn{K_j}{K[j]} with respect to the choice
        of the neighbor system and possibly the order of the conditional
        autoregression in the prior of transformed weights (see below).
      }
      \item{izero}{subscript \eqn{k_1\dots k_q}{k[1]...k[q]} 
        of the knot whose transformed weight 
        \eqn{a_{k_1\dots k_q}}{a[k1...kq]} will constantly be equal to zero.
        This is here for identifiability. To avoid numerical problems
        it is highly recommended to set \code{izero=rep(0, q)}.
	\code{izero[j]} should be taken from the set
	\eqn{-K_j, \dots, K_j}{-K[j],...,K[j]}.
      }
      \item{neighbor.system}{identification of the neighboring system
	for the Markov random field prior of transformed mixture weights
        \eqn{a_{k_1\,k_2}}{a[k1, k2]}. This can be substring of one of
	the following strings:
        \describe{
          \item{\code{uniCAR}}{ ``univariate conditional
	    autoregression'': a~prior based on squared differences of
	    given order \eqn{m}{m} (see argument \code{order}) in each row and
	    column. 

	    For univariate smoothing:
	    \deqn{
	      p(a) \propto \exp\Bigl\{-\frac{\lambda}{2}\sum_{k=-K+m}^K\bigl(\Delta^m a_{k}\bigr)^2\Bigr\},}{%
	      p(a) propto exp(-lambda/2 * sum[k=-K+m][K] (Delta^m a[k])^2),
	    }
	    where \eqn{\Delta^m}{Delta^m} denotes the difference
	    operator of order \eqn{m}{m}, i.e.
	    \eqn{\Delta^1 a_k = a_k - a_{k-1}}{Delta^1 a[k] = a[k] -
	      a[k-1]} and
    	    \eqn{\Delta^m a_k = \Delta^{m-1}a_k - \Delta^{m-1}a_{k-1},}{%
	      Delta^m a[k] = Delta^(m-1)a[k] - Delta^(m-1)a[k-1],}
	    \eqn{m \geq 2.}{m >= 2.}
	    

            For bivariate smoothing:
	    \deqn{
	      p(a) \propto
	      \exp\Bigl\{
	      -\frac{\lambda_1}{2}\sum_{k_1=-K_1}^{K_1}\sum_{k_2=-K_2+m}^{K_2}
	      \bigl(\Delta_1^m a_{k_1,k_2}\bigr)^2
	      -\frac{\lambda_2}{2}\sum_{k_2=-K_2}^{K_2}\sum_{k_1=-K_1+m}^{K_1}
	      \bigl(\Delta_2^m a_{k_1,k_2}\bigr)^2
	      \Bigr\},}{%
	      p(a) propto
	      exp(
	      -lambda[1]/2 * sum[k[1]=-K[1]][K[1]]sum[k[2]=-K[2]+m][K[2]]
	      (Delta[1]^m a[k[1], k[2]])^2
	      -lambda[2]/2 * sum[k[2]=-K[2]][K[2]]sum[k[1]=-K[1]+m][K[1]]
	      (Delta[2]^m a[k[1], k[2]])^2),
	    }
            where \eqn{\Delta_l^m}{Delta[l]^m} denotes the difference
	    operator of order \eqn{m}{m} acting in the \eqn{l}{l}th
	    margin, e.g.
	    \deqn{\Delta_1^2 = a_{k_1,k_2} - 2a_{k_1,k_2-1} + a_{k_1,k_2-2}.}{%
               Delta[1]^2 = a[k[1], k[2]] - 2*a[k[1], k[2]-1] + a[k[1], k[2]-2].
	     }
	     The precision parameters \eqn{\lambda_1}{lambda[1]} and
	     \eqn{\lambda_2}{lambda[2]} might be forced to be equal (see
	     argument \code{equal.lambda}.)
	  }
	  \item{\code{eight.neighbors}}{this prior is based on eight
	    nearest neighbors (i.e. except on edges, each full
	    conditional depends only on eight nearest neighbors) and
	    local quadratic smoothing.
	    It applies only in the case of bivariate smoothing. The prior is then
	    defined as
	    \deqn{
              p(a) \propto \exp
	      \Bigl\{-\frac{\lambda}{2}\sum_{k_1=-K_1}^{K_1-1}\sum_{k_2=-K_2}^{K_2-1}
	      \bigl(\Delta a_{k_1,k_2} \bigr)^2\Bigr\},
	    }{%
             p(a) propto exp(-lambda/2 *
	     sum[k[1]=-K[1]][K[1]-1]sum[k[2]=-K[2]][K[2]-1] (Delta a[k[1],k[2]])^2),
	    }  
            where
	    \deqn{\Delta a_{k_1,k_2} = a_{k_1,k_2} - a_{k_1+1,k_2} -
	      a_{k_1, k_2+1} + a_{k_1+1,k_2+1}.}{%
              Delta a[k[1],k[2]] = a[k[1],k[2]] - a[k[1]+1, k[2]] -
	      a[k[1], k[2]+1] + a[k[1]+1, k[2]+1].
	    }	    
	  }
	  \item{\code{twelve.neighbors}}{!!! THIS FEATURE HAS NOT BEEN
	    IMPLEMENTED YET. !!!}
        }
      }	
      \item{order}{order of the conditional autoregression if
	\code{neighbor.system = uniCAR}. Implemented are 1, 2, 3. If
	\code{order = 0} and  \code{neighbor.system = uniCAR}
	then mixture weights are assumed to be fixed
	and equal to their initial values specified by the \code{init}
	parameter (see below). Note that the numbers \eqn{K_j,}{K[j],}
        \eqn{j=1,\dots,q}{j=1,...,q} must be all equal to or higher than
        \code{order}.
      }
      \item{equal.lambda}{\code{TRUE/FALSE} applicable in the case when
	a density of bivariate observations is estimated and
	\code{neighbor.system = uniCAR}. It specifies whether there is
	only one common Markov random field precision parameter
	\eqn{\lambda}{lambda} for all margins (dimensions) or whether
	each margin (dimension) has its own precision parameter
	\eqn{\lambda}{lambda}.
	For all other neighbor systems is \code{equal.lambda}
	automatically \code{TRUE}.
      }
      \item{prior.lambda}{specification of the prior distributions for
	the Markov random field precision parameter(s)
	\eqn{\lambda}{lambda} (when \code{equal.lambda = TRUE}) or
	\eqn{\lambda_1,\dots,\lambda_q}{lambda[1],...,lambda[q]}
	(when \code{equal.lambda = TRUE}).
	This is a vector of substring of one of the following strings
	(one substring for each margin if \code{equal.lambda = FALSE},
	otherwise just one substring):
	\describe{
	  \item{\code{fixed}}{the \eqn{\lambda}{lambda} parameter is then
	    assumed to be fixed and equal to its initial values given
	    by \code{init} (see below).}
	  \item{\code{gamma}}{a particular \eqn{\lambda}{lambda} parameter has
	    a priori gamma distribution with shape \eqn{g_j}{g[j]} and rate
	    (inverse scale) \eqn{h_j}{h[j]} where
	    \eqn{j=1}{j=1} if \code{equal.lambda=TRUE} and
	    \eqn{j=1,\dots,q}{j=1,...,q} if \code{equal.lambda=TRUE}.
	    Shape and rate parameters are specified by
	    \code{shape.lambda}, \code{rate.lambda} (see below).}
	  \item{\code{sduniform}}{a particular
	    \eqn{1/\sqrt{\lambda}}{1/sqrt(lambda)} parameter
	    (i.e.a standard deviation of the Markov random field)
	    has a priori a uniform distribution on the interval
	    \eqn{(0, S_j)}{(0, S[j])} where
	    \eqn{j=1}{j=1} if \code{equal.lambda=TRUE} and
	    \eqn{j=1,\dots,q}{j=1,...,q} if \code{equal.lambda=TRUE}.
	    Upper limit of intervals is specified by \code{rate.lambda}
	    (see below).}
        }
      }
      \item{prior.gamma}{specification of the prior distribution for a
	reference knot (intercept) \eqn{\gamma}{gamma} in each
	dimension. This is a vector of substrings of one of the following strings
	(one substring for each margin):
	\describe{
	  \item{\code{fixed}}{the \eqn{\gamma}{gamma} parameter is then
	    assumed to be fixed and equal to its initial values given
	    by \code{init} (see below).}
	  \item{\code{normal}}{the \eqn{\gamma}{gamma} parameter has a
	    priori a normal distribution with mean and variance given by
	    \code{mean.gamma} and \code{var.gamma}.}
        }	  
      }	
      \item{prior.sigma}{specification of the prior distribution for
	basis standard deviations of the G-spline in each dimension.
	This is a vector of substrings of one of the following strings
	(one substring for each margin):
	\describe{
	  \item{\code{fixed}}{the \eqn{\sigma}{sigma} parameter is then
	    assumed to be fixed and equal to its initial values given
	    by \code{init} (see below).}
	  \item{\code{gamma}}{a particular \eqn{\sigma^{-2}}{sigma^{-2}} parameter has
	    a priori gamma distribution with shape \eqn{\zeta_j}{zeta[j]} and rate
	    (inverse scale) \eqn{\eta_j}{eta[j]} where
	    \eqn{j=1,\dots,q}{j=1,...,q}.
	    Shape and rate parameters are specified by
	    \code{shape.sigma}, \code{rate.sigma} (see below).}
	  \item{\code{sduniform}}{a particular
	    \eqn{\sigma}{sigma} parameter
	    has a priori a uniform distribution on the interval
	    \eqn{(0, S_j)}{(0, S[j]) where
	    \eqn{j=1,\dots,q}{j=1,...,q}}.
	    Upper limit of intervals is specified by \code{rate.sigma}
	    (see below).}
        }
      }
      \item{prior.intercept}{specification of the prior distribution for
	the intercept terms
	\eqn{\alpha_1,\dots,\alpha_q}{alpha[1],...,alpha[q]}
	(2nd specification) in each dimension.
        This is a vector of substrings of one of the following strings
	(one substring for each margin):
	\describe{
	  \item{\code{fixed}}{the intercept parameter is then
	    assumed to be fixed and equal to its initial values given
	    by \code{init} (see below).}
	  \item{\code{normal}}{the intercept parameter has a
	    priori a normal distribution with mean and variance given by
	    \code{mean.intercept} and \code{var.intercept}.}
        }	  	
      }
      \item{prior.scale}{specification of the prior distribution for
        the scale parameter (2nd specification) of the G-spline in each dimension
	This is a vector of substrings of one of the following strings
	(one substring for each margin):
	\describe{
	  \item{\code{fixed}}{the \code{scale} parameter is then
	    assumed to be fixed and equal to its initial values given
	    by \code{init} (see below).}
	  \item{\code{gamma}}{a particular \eqn{scale^{-2}}{scale^{-2}} parameter has
	    a priori gamma distribution with shape \eqn{\zeta_j}{zeta[j]} and rate
	    (inverse scale) \eqn{\eta_j}{eta[j]} where
	    \eqn{j=1,\dots,q}{j=1,...,q}.
	    Shape and rate parameters are specified by
	    \code{shape.scale}, \code{rate.scale} (see below).}
	  \item{\code{sduniform}}{a particular
	    \eqn{scale}{scale} parameter
	    has a priori a uniform distribution on the interval
	    \eqn{(0, S_j)}{(0, S[j]) where
	    \eqn{j=1,\dots,q}{j=1,...,q}}.
	    Upper limit of intervals is specified by \code{rate.scale}
	    (see below).}
        }
      }      
      \item{c4delta}{values of \eqn{c_1,\dots,c_q}{c[1],...,c[q]} which serve
        to compute the distance \eqn{\delta_j}{delta[j]}
        between two consecutive knots in each dimension.
        The knot \eqn{\mu_{j\,k},}{mu[j,k],}
        \eqn{j=1,\dots,q,}{j=1,...,q,} 
        \eqn{k=-K_j,\dots,K_j}{k=-K[j],...,K[j]} is defined
        as \eqn{\mu_{j\,k} = \gamma_j + k\,\delta_j}{mu[j,k] = gamma[j] + k*delta[j]}
        with \eqn{\delta_j = c_j\,\sigma_j}{delta[j] = c[j]*sigma[j]}.
      }
      \item{mean.gamma}{these are means for the normal prior
	distribution of middle knots
	\eqn{\gamma_1,\dots,\gamma_q}{gamma[1],...,gamma[q]}
	in each dimension if this prior is normal. For fixed
	\eqn{\gamma}{gamma} an appropriate element of the vector
	\code{mean.gamma} may be whatever.
      }
      \item{var.gamma}{these are variances for the normal prior
	distribution of middle knots
	\eqn{\gamma_1,\dots,\gamma_q}{gamma[1],...,gamma[q]}
	in each dimension if this prior is normal. For fixed
	\eqn{\gamma}{gamma} an appropriate element of the vector
	\code{var.gamma} may be whatever.
      }
      \item{shape.lambda}{these are shape parameters for the gamma prior
	(if used) of Markov random field precision parameters
	\eqn{\lambda_1,\dots,\lambda_q}{lambda[1],...,lambda[q]} (if
	\code{equal.lambda = FALSE}) or \eqn{\lambda_1}{lambda[1]} (if
	\code{equal.lambda = TRUE}).
      }
      \item{rate.lambda}{these are rate parameters for the gamma prior
	(if \code{prior.lambda = gamma}) of Markov random field precision parameters
	\eqn{\lambda_1,\dots,\lambda_q}{lambda[1],...,lambda[q]} (if
	\code{equal.lambda = FALSE}) or \eqn{\lambda_1}{lambda[1]} (if
	\code{equal.lambda = TRUE}) or upper limits of the uniform prior
        (if \code{prior.lambda = sduniform}) of Markov random field
	standard deviation  parameters
	\eqn{\lambda_1^{-1/2},\dots,\lambda_q^{-1/2}}{lambda[1]^{-1/2},...,lambda[q]^{-1/2}} (if
	\code{equal.lambda = FALSE}) or \eqn{\lambda_1^{-1/2}}{lambda[1]^{-1/2}} (if
	\code{equal.lambda = TRUE}).
      }
      \item{shape.sigma}{these are shape parameters for the gamma prior
	(if used) of basis inverse variances
	\eqn{\sigma_1^{-2},\dots,\sigma_q^{-2}}{sigma[1]^{-2},...,sigma[q]^{-2}}.
      }
      \item{rate.sigma}{these are rate parameters for the gamma prior
	(if \code{prior.sigma = gamma}) of basis inverse variances
	\eqn{\sigma_1^{-2},\dots,\sigma_q^{-2}}{sigma[1]^{-2},...,sigma[q]^{-2}} 
	or upper limits of the uniform prior
        (if \code{prior.sigma = sduniform}) of basis standard deviations
	\eqn{\sigma_1,\dots,\sigma_q}{sigma[1],...,sigma[q]}.
      }
      \item{mean.intercept}{these are means for the normal prior
	distribution of the G-spline intercepts (2nd specification)
	\eqn{\alpha_1,\dots,\alpha_q}{alpha[1],...,alpha[q]}
	in each dimension if this prior is normal. For fixed
	\eqn{\alpha}{alpha} an appropriate element of the vector
	\code{mean.intercept} may be whatever.
      }
      \item{var.intercept}{these are variances for the normal prior
	distribution of the G-spline intercepts
	\eqn{\alpha_1,\dots,\alpha_q}{alpha[1],...,alpha[q]}
	in each dimension if this prior is normal. For fixed
	\eqn{\alpha}{alpha} an appropriate element of the vector
	\code{var.alpha} may be whatever.
      }      
      \item{shape.scale}{these are shape parameters for the gamma prior
	(if used) of the G-spline scale parameter (2nd specification)
	\eqn{scale_1^{-2},\dots,scale_q^{-2}}{scale[1]^{-2},...,scale[q]^{-2}}.
      }
      \item{rate.scale}{these are rate parameters for the gamma prior
	(if \code{prior.scale = gamma}) of the G-spline inverse variances
	\eqn{scale_1^{-2},\dots,scale_q^{-2}}{scale[1]^{-2},...,scale[q]^{-2}} 
	or upper limits of the uniform prior
        (if \code{prior.scale = sduniform}) of the G-spline scale
	\eqn{scale_1,\dots,scale_q}{scale[1],...,scale[q]}.
      }      
    }
  }
  \item{init}{a list of the initial values to start the McMC. Set to
    \code{NULL} such parameters that you want the program should itself sample
    for you or parameters that are not needed in your model.
    \describe{
      \item{iter}{the number of the iteration to which the initial values
        correspond, usually zero.}      
      \item{a}{vector/matrix of initial transformed mixture weights
	\eqn{a_{k_1},}{a[k1],} \eqn{k_1=-K_1,\dots,K_1}{k1=-K1,...,K1}
	if univariate density is estimated;
        \eqn{a_{k_1\,k_2},}{a[k1,k2],}
	\eqn{k_1=-K_1,\dots,K_1,}{k1=-K1,...,K1,}
	\eqn{k_2=-K_2,\dots,K_2,}{k2=-K2,...,K2,}
	if bivariate density is estimated. This initial value can be
	guessed by the function itself.
      }
      \item{lambda}{initial values for Markov random field precision
        parameter(s) \eqn{\lambda}{lambda} (if \code{equal.lambda = TRUE}),
        \eqn{\lambda_1,\dots,\lambda_q}{lambda[1],...,lambda[q]} 
        (if \code{equal.lambda = FALSE}.)
      }
      \item{gamma}{initial values for the middle knots in each
	dimension.

	If \code{prior$specification = 2} it is recommended
	(for easier interpretation of the results) to set
	\code{init$gamma} to zero for all dimensions.

	If \code{prior$specification = 1} \code{init$gamma} should be
	approximately equal to the mean value of the data in each margin.
      }
      \item{sigma}{initial values for basis standard deviations in each
        dimension.

	If \code{prior$specification = 2} this should be
	approximately equal to the range of standardized data (let say 4
	+ 4) divided by the number of knots in each margin and
	multiplied by something like 2/3.

	If \code{prior$specification = 1} this should be approximately equal to the range of your
	data divided by the number of knots in each margin and
	multiplied again by something like 2/3.
      }
      \item{intercept}{initial values for the intercept term in each
	dimension.

	Note that if \code{prior$specification = 1} this
	initial value is always changed to zero for all dimensions.	
      }
      \item{scale}{initial values for the G-spline scale parameter in each
        dimension.

	Note that if \code{prior$specification = 1} this
	initial value is always changed to one for all dimensions.
      }      
      \item{y}{initial values for (possibly unobserved censored) observations. This
	should be either a vector of length equal to the sample size if
	the response is univariate or a matrix with as many rows as is
	the sample size and two columns if the response is
	bivariate. Be aware that \code{init$y} must be consistent with
	data supplied. This initial can be guessed by the function
	itself. Possible missing values in \code{init$y} tells the
	function to guess the initial value.
      }
      \item{r}{initial values for labels of components to which the
	(augmented) observations belong. This initial can be guessed by
	the function itself. This
	should be either a vector of length equal to the sample size if
	the response is univariate or a matrix with as many rows as is
	the sample size and two columns if the response is
	bivariate. Values in the first column of this matrix should be
	between \code{-prior$K[1]} and \code{prior$K[1]}, values in the
	second column of this matrix between \code{-prior$K[2]} and
	\code{prior$K[2]}, e.g. when \code{init$r[i,1:2] = c(-3, 6)} it
	means that the \eqn{i}{i}th observation is initially assigned to
	the component with the mean
	\eqn{\boldsymbol{\mu}=(\mu_1, \mu_2)'}{mu = (mu[1], mu[2])'} where
	\deqn{\mu_1 = \mu_{1,\,-3} = \gamma_1 -3\,c_1\sigma_1}{%
	  mu[1] = mu[1,-3] = gamma[1] -3*c[1]*sigma[1]}
	and
	\deqn{\mu_2 = \mu_{2,\,6} = \gamma_2 +6\,c_2\sigma_2.}{%
	  mu[2] = mu[1,6] = gamma[2] +6*c[2]*sigma[2].}
      }
    }
  }
  \item{mcmc.par}{a list specifying further details of the McMC
    simulation. There are default values implemented for all components
    of this list.
    \describe{
      \item{type.update.a}{it specifies the McMC method to update
	transformed mixture weights \eqn{a}{a}. It is a~substring of one
	of the following strings:
	\describe{
	  \item{slice}{slice sampler of Neal (2003) is used (default choice);}
	  \item{ars.quantile}{adaptive rejection sampling of Gilks and
	    Wild (1992) is used with starting abscissae equal to 15\%,
	    50\% and 85\% quantiles of a~piecewise exponential
	    approximation to the full conditional from the previous
	    iteration;}
	  \item{ars.mode}{adaptive rejection sampling of Gilks and Wild
	    (1992) is used with starting abscissae equal to the mode and
	    plus/minus twice approximate standard deviation of the full
	    conditional distribution}
	}
      }
      \item{k.overrelax.a}{this specifies a frequency of overrelaxed
	updates of transformed mixture weights \eqn{a}{a} when slice
	sampler is used. Every \eqn{k}{k}th value is sampled in a usual
	way (without overrelaxation). If you do not want overrelaxation
	at all, set \code{k.overrelax.a} to 1 (default choice). Note
	that overrelaxation can be only done with the slice sampler (and
        not with adaptive rejection sampling).}
      \item{k.overrelax.sigma}{a vector of length equal to the dimension
        of the G-spline specifying a frequency of overrelaxed updates of
        basis G-spline variances. If you do not want overrelaxation at
	all, set all components of \code{k.overrelax.sigma} to 1
	(default choice).}
      \item{k.overrelax.scale}{a vector of length equal to the dimension
        of the G-spline specifying a frequency of overrelaxed updates of
        the G-spline scale parameters (2nd specification).
	If you do not want overrelaxation at
	all, set all components of \code{k.overrelax.scale} to 1
	(default choice).}      
    }       
  }
  \item{store}{a~list of logical values specifying which chains that are
    not stored by default are to be stored. The list can have the
    following components.
    \describe{
      \item{a}{if \code{TRUE} then all the transformed mixture weights
	\eqn{a_{k_1,\,k_2},}{a[k[1],k[2]],}
	\eqn{k_1=-K_1,\dots,K_1,}{k[1]=-K[1],..., K[1],}
	\eqn{k_2=-K_2,\dots,K_2,}{k[2]=-K[2],..., K[2],}
	related to the G-spline are stored.}
      \item{y}{if \code{TRUE} then augmented log-event times for all
	observations are stored.}
      \item{r}{if \code{TRUE} then labels of mixture components for
	residuals are stored.}
    }  %% end of describe  
  }  %% end of item{store}
  \item{dir}{a string that specifies a directory where all sampled
    values are to be stored.
  }   
}
\value{
  A list of class \code{bayesHistogram} containing an information
  concerning the initial values and prior choices.
}
\section{Files created}{  
  Additionally, the following files with sampled values
  are stored in a directory specified by \code{dir} argument of this
  function (some of them are created only on request, see \code{store}
  parameter of this function).

  Headers are written to all files created by default and to files asked
  by the user via the argument \code{store}.  All
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
      \code{D.1.1}, \code{D.2.1}, \code{D.2.2} in the bivariate case and
      columns labeled \code{k}, \code{Mean.1}, \code{D.1.1} in the
      univariate case, where
      
      \bold{k} = number of mixture components that had probability
      numerically higher than zero;
      
      \bold{Mean.1} =
      \eqn{\mbox{E}(Y_{i,1})}{E(Y[i,1])};
      
      \bold{Mean.2} =
      \eqn{\mbox{E}(Y_{i,2})}{E(Y[i,2])};
      
      \bold{D.1.1} =
      \eqn{\mbox{var}(Y_{i,1})}{var(Y[i,1])};
      
      \bold{D.2.1} =
      \eqn{\mbox{cov}(Y_{i,1},\,Y_{i,2})}{cov(Y[i,1], Y[i,2])};

      \bold{D.2.2} =
      \eqn{\mbox{var}(Y_{i,2})}{var(Y[i,2])}.
    }
    \item{mweight.sim}{sampled mixture weights
      \eqn{w_{k_1,\,k_2}}{w[k[1],k[2]]} of mixture components that had
      probabilities numerically higher than zero. 
    }
    \item{mmean.sim}{indeces \eqn{k_1,\;k_2,}{k[1], k[2],}
      \eqn{k_1 \in\{-K_1, \dots, K_1\},}{k[1] in {-K[1], ..., K[1]},}
      \eqn{k_2 \in\{-K_2, \dots, K_2\}}{k[2] in {-K[2], ..., K[2]}}
      of mixture components that had probabilities numerically higher
      than zero. It corresponds to the weights in
      \code{mweight.sim}. 
    }
    \item{gspline.sim}{characteristics of the sampled G-spline
      (distribution of
      \eqn{(Y_{i,1},\,Y_{i,2})'}{(Y[i,1], Y[i,2])'}).
      This file together with \code{mixmoment.sim},
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
    \item{mlogweight.sim}{fully created only if \code{store$a = TRUE}. The
      file contains the transformed weights
      \eqn{a_{k_1,\,k_2},}{a[k[1],k[2]],}
      \eqn{k_1=-K_1,\dots,K_1,}{k[1]=-K[1],..., K[1],}
      \eqn{k_2=-K_2,\dots,K_2}{k[2]=-K[2],..., K[2]} of all mixture
      components, i.e. also of components that had numerically zero
      probabilities.
    }
    \item{r.sim}{fully created only if \code{store$r = TRUE}. The file
      contains the labels of the mixture components into which the
      observations are intrinsically assigned. Instead of double indeces
      \eqn{(k_1,\,k_2)}{(k[1], k[2])}, values from 1 to \eqn{(2\,K_1+1)\times
	(2\,K_2+1)}{(2*K[1]+1)*(2*K[2]+1)} are stored here. Function
      \code{\link{vecr2matr}} can be used to transform it back to double
      indeces.
    }
    \item{lambda.sim}{either one column labeled \code{lambda} or two
      columns labeled \code{lambda1} and \code{lambda2}. These are the
      values of the smoothing parameter(s) \eqn{\lambda}{lambda}
      (hyperparameters of the prior distribution of the transformed
      mixture weights \eqn{a_{k_1,\,k_2}}{a[k[1],k[2]]}). 
    }
    \item{Y.sim}{fully created only if \code{store$y = TRUE}. It
      contains sampled (augmented) log-event times for all observations
      in the data set.
    }
    \item{logposter.sim}{columns labeled \code{loglik}, \code{penalty}
      or \code{penalty1} and
      \code{penalty2}, \code{logprw}. The columns have the following
      meaning (the formulas apply for the bivariate case).

      \bold{loglik}
      \eqn{=}{=} \eqn{%
	-N\Bigl\{\log(2\pi) + \log(\sigma_1) + \log(\sigma_2)\Bigr\}-
          0.5\sum_{i=1}^N\Bigl\{
	  (\sigma_1^2\,\tau_1^2)^{-1}\; (y_{i,1} - \alpha_1 - \tau_1\mu_{1,\,r_{i,1}})^2 +
	  (\sigma_2^2\,\tau_2^2)^{-1}\; (y_{i,2} - \alpha_2 - \tau_2\mu_{2,\,r_{i,2}})^2
          \Bigr\}
      }{%
	-N(log(2*pi) + log(sigma[1]) + log(sigma[2]))
          -0.5*sum[i=1][N](
          (sigma[1]^2*tau[1]^2)^(-1) * (y[i,1] - alpha[1] - tau[1]*mu[1,r[i,1]])^2 +
          (sigma[2]^2*tau[2]^2)^(-1) * (y[i,2] - alpha[2] - tau[2]*mu[2,r[i,2]])^2 	  
	)}
      
      where \eqn{y_{i,l}}{y[i,l]} denotes (augmented) \emph{(i,l)}th
      true log-event time. In other words, \code{loglik} is equal to the
      conditional log-density
      \eqn{\sum_{i=1}^N\,\log\Bigl\{p\bigl((y_{i,1},\,y_{i,2})\;\big|\;r_{i},\,\mbox{G-spline}\bigr)\Bigr\};}{%
      sum[i=1][N] log(p((y[i,1], y[i,2]) | r[i], G-spline));
      }

      \bold{penalty1:} If \code{prior$neighbor.system} = \code{"uniCAR"}:
        the penalty term for the first dimension not multiplied by
	\code{lambda1};
      
      \bold{penalty2:} If \code{prior$neighbor.system} = \code{"uniCAR"}:
        the penalty term for the second dimension not multiplied by
	\code{lambda2};

      \bold{penalty:} If \code{prior$neighbor.system} is different from \code{"uniCAR"}:
        the penalty term not multiplied by \code{\lambda};

      \bold{logprw} \eqn{=}{=}
      \eqn{-2\,N\,\log\bigl\{\sum_{k_1}\sum_{k_2}a_{k_1,\,k_2}\bigr\} +
	\sum_{k_1}\sum_{k_2}N_{k_1,\,k_2}\,a_{k_1,\,k_2},}{%
	-2*N*log(sum[k[1]]sum[k[2]] exp(a[k[1],k[2]])) +
	sum[k[1]]sum[k[2]] N[k[1],k[2]]*a[k[1],k[2]],}
      where \eqn{N_{k_1,\,k_2}}{N[k[1],k[2]]} is the number of observations
      assigned intrinsincally to the \eqn{(k_1,\,k_2)}{(k[1], k[2])}th
      mixture component.

      In other words, \code{logprw} is equal to the conditional
      log-density
      \eqn{\sum_{i=1}^N \log\bigl\{p(r_i\;|\;\mbox{G-spline
	  weights})\bigr\}.}{%
      sum[i=1][N] log(p(r[i] | G-spline weights)).}
    }
  }
}
\references{
  Gilks, W. R. and Wild, P. (1992).
  Adaptive rejection sampling for Gibbs sampling.
  \emph{Applied Statistics,} \bold{41}, 337 - 348.

  \eqn{\mbox{Kom\'{a}rek, A.}}{Kom&#225rek, A.} (2006).
  \emph{Accelerated Failure Time Models for Multivariate
    Interval-Censored Data with Flexible Distributional Assumptions}.
  PhD. Thesis, Katholieke Universiteit Leuven, Faculteit Wetenschappen.

  \eqn{\mbox{Kom\'{a}rek, A.}}{Kom&#225rek, A.} and Lesaffre, E. (2007).
  Bayesian accelerated failure time model with multivariate doubly-interval-censored data
  and flexible distributional assumptions.
  \emph{To appear in Journal of the American Statistical Association.}

  \eqn{\mbox{Kom\'{a}rek, A.}}{Kom&#225rek, A.} and Lesaffre, E. (2006b).
  Bayesian semi-parametric accelerated failurew time model for paired
  doubly interval-censored data.
  \emph{Statistical Modelling}, \bold{6}, 3--22.
  
  Neal, R. M. (2003).
  Slice sampling (with Discussion).
  \emph{The Annals of Statistics,} \bold{31}, 705 - 767.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{komarek@karlin.mff.cuni.cz}
}
%\examples{}
\keyword{smooth}
\keyword{survival}
\keyword{multivariate}
