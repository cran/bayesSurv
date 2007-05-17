\name{bayessurvreg1}
\alias{bayessurvreg1}
\title{
  A Bayesian survival regression with an error distribution
  expressed as a~normal mixture with unknown number of components
}
\description{
  A function to sample from the posterior distribution for a survival
  regression  model
  \deqn{\log(T_{i,l}) = \beta^T x_{i,l} + b_i^T z_{i,l} +
    \varepsilon_{i,l},\quad i=1,\dots,N,\ l=1,\dots,n_i,}{%
    log(T[i,l]) = beta'*x[i,l] + b[i]'*z[i,l] + epsilon[i,l], i=1,...,N,
    l=1,...,n[i],}
  where distribution of \eqn{\varepsilon_{i,l}}{epsilon[i,l]} is specified
  as a normal mixture with unknown number of components as in Richardson
  and Green (1997) and random effect \eqn{b_i}{b[i]} is normally distributed.

  See \eqn{\mbox{Kom\'{a}rek}}{Komarek} (2006) or
  \eqn{\mbox{Kom\'{a}rek}}{Komarek} and Lesaffre (2007)
  for more detailed description of prior assumptions.
  
  Sampled values are stored on a disk to be further worked out by e.g.
  \code{coda} or \code{boa}.
}
\usage{
bayessurvreg1(formula, random,
   data = parent.frame(), subset,
   na.action = na.fail,
   x = FALSE, y = FALSE, onlyX = FALSE,
   nsimul = list(niter = 10, nthin = 1, nburn = 0,
                 nnoadapt = 0, nwrite = 10),
   prior = list(kmax = 5, k.prior = "poisson", poisson.k = 3,
                dirichlet.w = 1,
                mean.mu = NULL, var.mu = NULL,
                shape.invsig2 = 1.5,
                shape.hyper.invsig2 = 0.8, rate.hyper.invsig2 = NULL,
                pi.split = NULL, pi.birth = NULL,
                Eb0.depend.mix = FALSE),
   prior.beta, prior.b, prop.revjump,
   init = list(iter = 0, mixture = NULL, beta = NULL,
               b = NULL, D = NULL,
               y = NULL, r = NULL, otherp = NULL, u = NULL),
   store = list(y = TRUE, r = TRUE, b = TRUE, u = TRUE,
                MHb = FALSE, regresres = FALSE),
   dir = getwd(),
   toler.chol = 1e-10, toler.qr = 1e-10, \dots)
}
\arguments{
  \item{formula}{model formula for the `fixed' part of the model, i.e. the
    part that specifies \eqn{\beta^T x_{i,l}}{beta*x[i,l]}. See
    \code{\link[survival]{survreg}} for further details. Intercept is implicitely
    included in the model by estimation of the error distribution. As a
    consequence \code{-1} in the model formula does not have any effect
    on the model.

    The left-hand side of the \code{formula} must be an~objecy created
    using \code{\link[survival]{Surv}}.
    
    If \code{random} is used then the formula must contain
    an identification of clusters in the form \code{cluster(id)}, where
    \code{id} is a name of the variable that determines clusters, e.g.
    \tabular{c}{
      \code{Surv(time, event)~gender + cluster(id)}.
    }  
  }
  \item{random}{formula for the `random' part of the model, i.e. the
    part that specifies \eqn{b_i^T z_{i,l}}{b[i]'*z[i,l]}. If omitted,
    no random part is included in the model. E.g. to specify the model with a
    random intercept, say \code{random=~1}. All effects mentioned in
    \code{random} should also be mentioned on the right-hand side of
    \code{formula}.

    When some random effects are included the random intercept is added
    by default. It can be removed using e.g. \code{random=~-1 + gender}.
  }
  \item{data}{optional data frame in which to interpret the variables
    occuring in the formulas.}
  \item{subset}{subset of the observations to be used in the fit.}
  \item{na.action}{function to be used to handle any \code{NA}s in the
    data. The user is discouraged to change a default value
    \code{na.fail}.}
  \item{x}{if \code{TRUE} then the \code{X} matrix is returned. This
    matrix contain all columns appearing in both \code{formula} and
    \code{random} parameters.}  
  \item{y}{if \code{TRUE} then the \code{y} matrix (of log-survival
    times) is returned.}
  \item{onlyX}{if TRUE, no McMC is performed. The function returns only
    a design matrix of your model (intercept excluded). It might be
    useful to set up correctly a parameter for a block update of
    \eqn{\beta}{beta} (regression parameters related to the fixed
    effects) and \eqn{\gamma}{gamma} (means of the random effects,
    random intercept excluded) parameters in the model if
    Metropolis-Hastings is to be used instead of default Gibbs.}
  \item{nsimul}{a list giving the number of iterations of the McMC and
    other parameters of the simulation.
    \describe{
      \item{niter}{total number of sampled values after discarding
	thinned ones, burn-up included.}
      \item{nthin}{thinning interval.}
      \item{nburn}{number of sampled values in a burn-up period after
	discarding thinned values. This value should be smaller than
	\code{niter}. If not, \code{nburn} is set to \code{niter - 1}. It can be set to zero.}
      \item{nnoadapt}{applicable if some blocks of parameters are
	updated using an adaptive Metropolis algorithm. This is a number
        of sampled values that are generated using an initial and fixed
	proposal covariance matrix. It should be smaller or equal to
	\code{nburn}. If not, \code{nnoadapt} is set to \code{nburn}.}
      \item{nwrite}{an interval at which sampled values are written to
	files. }
    }
  }    
  \item{prior}{a list that identifies prior hyperparameters and prior
    choices. See accompanying paper for more details.
    Some prior parameters can be guessed by the function itself. If you
    want to do so, set such parameters to \code{NULL}. Set to
    \code{NULL} also the parameters that are not needed in your model.
    \describe{
      \item{kmax}{value of \eqn{k_{max}}{k[max]}, upper limit for the
	number of mixture components. Its high values like 100 will
	usually correspond to \eqn{\infty}{infinity}.}
      \item{k.prior}{a string specifying the prior distribution of
	\eqn{k}{k}, number of mixture components. Valid are either
	``poisson'', ``uniform'', or ``fixed''. When ``fixed'' is given
	then the number of mixture components is not sampled.}
      \item{poisson.k}{prior hyperparameter \eqn{\lambda}{lambda} for
	the number of mixture components $k$ if Poisson prior for this
	parameter is used.}
      \item{dirichlet.w}{prior hyperparameter \eqn{\delta}{delta} for
	the Dirichlet distribution of mixture weights
	\eqn{w_1,\dots,w_k}{w[1],...w[k]}.}
      \item{mean.mu}{prior hyperparameter \eqn{\xi}{xi} for the mean of
	the normal prior for mixture means
	\eqn{\mu_1,\dots,\mu_k}{mu[1],...,mu[k]}.}
      \item{var.mu}{prior hyperparameter \eqn{\kappa}{kappa} for the
	variance of the normal prior for mixture means \eqn{\mu_1,\dots,\mu_k}{mu[1],...,mu[k]}.}
      \item{shape.invsig2}{prior hyperparameter \eqn{\zeta}{zeta} for
	the shape of the inverse-gamma distribution for the mixture
	variances
	\eqn{\sigma_1^2,\dots,\sigma_k^2}{sigma[1]^2,...,sigma[k]^2}.}
      \item{shape.hyper.invsig2}{prior hyperparameter (shape) \eqn{g}{g} for the
        gamma distribution of the parameter \eqn{\eta}{eta}. Remember,
	\eqn{\eta}{eta} is a scale parameter of the inverse-gamma distribution for the mixture
	variances
	\eqn{\sigma_1^2,\dots,\sigma_k^2}{sigma[1]^2,...,sigma[k]^2}.}
      \item{rate.hyper.invsig2}{prior hyperparameter (rate) \eqn{h}{h} for the
        gamma distribution of the parameter \eqn{\eta}{eta}. Remember,
	\eqn{\eta}{eta} is a scale parameter of the inverse-gamma distribution for the mixture
	variances
	\eqn{\sigma_1^2,\dots,\sigma_k^2}{sigma[1]^2,...,sigma[k]^2}.}
      \item{pi.split}{probabilities of a split move within the
	reversible jump McMC. It must be a vector of length equal to
	\code{kmax} with the first component equal to 1 and the last
	component equal to 0. If \code{NULL} 2nd to (k-1)th components
	are set to 0.5.}
      \item{pi.birth}{probabilities of a birth move within the
	reversible jump McMC. It must be a vector of length equal to
	\code{kmax} with the first component equal to 1 and the last
	component equal to 0. If \code{NULL} 2nd to (k-1)th components
	are set to 0.5.}
      \item{Eb0.depend.mix}{this will normally be \code{FALSE}. Setting
	this option to \code{TRUE} served for some experiments during
	the development of this function. In principle, when this is set
	to \code{TRUE} and the random intercept is included in the model
	then it is assumed that the mean of the random intercept is not
	zero but \eqn{\sum_{j=1}^k w_j\mu_j}{sum[j=1][k] w[j]*mu[j]},
	i.e. the mean of the random intercept depends on
	mixture. However, this did not werk too well.	
      }	
    }
  }
  \item{prior.beta}{a list defining the blocks of \eqn{\beta}{beta}
    parameters (both fixed effects and means of random effects, except
    the random intercept) that are to be updated together (in a block),
    a description of how they are updated and a specification of priors.
    The list is assumed to have the following components.
    \describe{
    \item{mean.prior}{a vector specifying a prior mean for each
      \eqn{\beta}{beta} parameter in the model.}
    \item{var.prior}{a vector specifying a prior variance for each
      \eqn{\beta}{beta} parameter. It is recommended to run the function
      \code{bayessurvreg1} first with its argument \code{onlyX} set to
      \code{TRUE} to find out how the \eqn{\beta}{beta}s are
      sorted. They must correspond to a design matrix X.}
    \item{blocks}{a list with the following components.
      \describe{
        \item{ind.block}{a list with vectors with indeces of columns of
         the design matrix X defining the effect of \eqn{\beta}{beta}s in the
         block. If not specified, all \eqn{\beta}{beta} parameters
	 corresponding to fixed effects are updated in one block and
	 remaining \eqn{\beta}{beta} parameters (means of random
	 effects) in the second block using the Gibbs move.}
       \item{cov.prop}{a list with vectors with a lower triangle of the covariance
         matrix which is used in the normal proposal (use a command
         \code{\link{lower.tri}} with \code{diag = TRUE} to get a lower
	 triangle from a matrix) when one of the Metopolis-like
	 algorithms is used for
         a given block. This matrix is used at each iteration if the given
         block is updated using a standard random-walk Metropolis-Hastings step. If the
         block is updated using an adaptive Metropolis step this matrix is
         used only at start. If not specified and Metropolis-like
	 algorith is required a diagonal matrix with prior
         variances for corresponding \eqn{\beta}{beta} on a diagonal is
         used. It is set to a vector of zeros of appropriate length when the Gibbs move is required
	 for a given block.}
      }
    }
    \item{type.upd}{a character vector specifying the type of the update
      that will be used for each block. Valid are substrings of
      either "gibbs" or "adaptive.metropolis" or "random.walk.metropolis".
      Default is "gibbs" for all blocks.}
    \item{mean.sampled}{a vector of means of up to now sampled
      values. This component is useful when the adaptive Metropolis
      algorithm is used and we do not start from the beginning
      (e.g. already several iterations of McMC have already been
      performed). Otherwise, this component does not have to be filled.}
    \item{eps.AM}{a vector with \eqn{\epsilon}{epsilon} from the
      adaptive Metropolis algorithm for each block.}
    \item{sd.AM}{a vector specifying \eqn{s_d, d = 1,\dots, D}{s[d], d =
      1, ..., D} numbers from the
      adaptive Metropolis algorithm for each dimension. This vector must
      be of length equal at least to the length of the longest
      block. Defaults values are \eqn{\frac{1}{d}2.4^2}{(1/d)*2.4^2}
      where \eqn{d}{d} denotes a length of the block.}
    \item{weight.unif}{a vector specifying the weight of the uniform
      component in the proposal for each block. If not specified, it is
      equal to 0.5 for all parameters.}
    \item{half.range.unif}{a vector of same length as the number of
      columns in the design matrix X specifying the half range of the
      uniform component of the proposal.}     
    }
  }
  \item{prior.b}{a list defining the way in which the random effects are
    to be updated and the specification of priors for random effects
    related parameters. The list is assumed to have following components.
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
      This is also the default choice.
      }
      \item{sduniform}{this can be used only when the random effect is
	univariate. Then the matrix \eqn{D}{D} is just a scalar and the
	prior of \eqn{\sqrt{D}}{sqrt(D)} (standard deviation of the
	univariate random effect) is assumed to be uniform on interval
	\eqn{(0, S)}{(0, S)}. The upper limit \eqn{S}{S} is given by \code{prior.b$scale.D}.
      }	
      }
    }
    \item{df.D}{degrees of freedom parameter \eqn{\tau}{tau} in the case
      that the prior of the matrix \eqn{D}{D} is inverse-Wishart.
    }
    \item{scale.D}{a lower triangle of the scale matrix \eqn{S}{S} in
      the case that the prior of the matrix \eqn{D}{D} is
      inverse-Wishart or the upper limit \eqn{S}{S} of the uniform distribution in
      the case that \eqn{\sqrt{D} \sim \mbox{Unif}(0, S)}{sqrt(D) ~
	Unif(0, S)}.}
    \item{type.upd}{a character vector specifying the type of the
      update. Valid are substrings of either "random.walk.metropolis" or
      "gibbs". Default is "gibbs". In contrast to \eqn{\beta}{beta}
      parameters, all random effects are updated using the same type of
      the move. If "random.walk.metropolis" is used, random effects may
      be divided into blocks in which they are updated. With "gibbs",
      there is only one block defined for all random effects.
      which are updated in one step using its full conditional distribution.}
    \item{blocks}{a list with the following components. This is set to
      NULL if \code{type.upd = "gibbs"}.      
      \describe{
         \item{ind.block}{a list with vectors with indeces of random
	   effects defining the block. Random intercept has always an
	   index 1, remaining random effects have subsequent indeces
	   according to their appearance in the design matrix X.}
         \item{cov.prop}{a list with vectors with a lower triangle of the covariance
           matrix which is used in the normal proposal (use a command
           \code{\link{lower.tri}} with \code{diag = TRUE} to get a lower triangle from a matrix) for
           a given block when
           \tabular{c}{
	     \code{type.upd = "random.walk.metropolis"}. 
	   }  
	 }      	
      }
    }	 
    \item{weight.unif}{a vector specifying the weight of the uniform
      component in the proposal for each block when
      \tabular{c}{
	\code{type.upd = "random.walk.metropolis"}. 
      }
      If not specified, it is
      equal to 0.5 for all parameters. It is set to NULL if \code{type.upd = "gibbs"}.} 
    \item{half.range.unif}{a vector of same length as the number of
      random effects specifying the half range of the uniform component
      of the proposal when \code{type.upd = "random.walk.metropolis"}.
      It is set to NULL if \code{type.upd = "gibbs"}.}       
  }
  }
  \item{prop.revjump}{a list of values defining in which way the
    reversible jumps will be performed. 
    \describe{
      \item{algorithm}{a string defining the algorithm used to generate
	canonical proposal vectors
	\eqn{u = (u_{3k+1},\dots, u_{3k_{max}})'}{u = (u[3k+1], ..., u[3kmax])}
	where \eqn{u_{3k+1}, u_{3k+2}, u_{3k+3}}{u[3k+1], u[3k+2], u[3k+3]}
	are directly used when a jump to a space of higher dimension is
	proposed. These canonical proposal vectors are further
	transformed to give desired parameters (mixture component's weight, mean and
	variance). 
	Valid values of \code{prop.revjump$algorithm} are substrings of
	"basic", "independent.av", "correlated.av". \code{"basic"} means
	that both components of vectors \eqn{u}{u} and vectors
	\eqn{u}{u} in time are generated independently from a standard
	uniform distribution. This corresponds to a basic reversible
	jumps McMC algorithm of Green (1995). Other two methods
	implement an auxiliary variable method of Brooks et
	al. (2003). The first one an independent auxiliary variable
	method where vectors \eqn{u}{u} may be correlated in time
	however their components are independent and the second one the
	correlated auxiliary method where vectors \eqn{u}{u} are
	correlated in time and also their components may be
	correlated. In both cases components of vectors \eqn{u}{u}
	follow marginally a standard uniform distribution. A moody ring
	method of Brooks et al. (2003) is used to generate \eqn{u}{u} vectors.
      }
      \item{moody.ring}{parameters for the moody ring when
	\code{algorithm} is either "independent.av" or
	"correlated.av". This is a two component vector with both
	components taking values between 0 and 0.5 defining the strength
        of a correlation in time and between the components of
	\eqn{u}{u} vectors. This vector is ignored when \code{algorithm
	  = "basic"}. The first component of this vector determines
	dependence between \eqn{u}{u} vectors in time
	(\eqn{\varepsilon}{epsilon} in Brooks et al. (2003)), the second
        component determines dependence between components of \eqn{u}{u}
	vectors (\eqn{\delta}{delta} in Brooks et al. (2003)). The
	second compoenent is ignored when \code{algorithm = "independent.av"}.
        Note that both \eqn{\varepsilon}{epsilon} and
	\eqn{\delta}{delta} do not have a meaning of correlation. They
	determine a range of additional uniform distributions. So that
	their values equal to \eqn{0}{0} mean perfect correlation and
	their values equal to \eqn{0.5}{0.5} mean
	independence. I.e. "correlated.av" with \eqn{\delta = 0.5}{delta
	= 0.5} is same as "independent.av" and "correlated.av" with
        \eqn{\delta = 0.5, \varepsilon = 0.5}{delta = 0.5, epsilon =
        0.5} is same as "basic".}
      \item{transform.split.combine}{a description of how the canonical
	variables \eqn{u}{u} are to be transformed to give new values of
        mixture component's weight, mean and variance when a split move
	is proposed. Possible values are substrings of
	"richardson.green", "brooks" and "identity". In all cases, the
	\eqn{(0, 1)}{(0, 1)} canonical variables \eqn{u}{u} are transformed 
	to \eqn{(0, 1)}{(0, 1)} variates \eqn{v}{v} that are than used
        to compute new values of mixture component's weight, mean and
	variance using a method of moments matching described in
	Richardson and Green (1997). When "identity", no further
	transformation is performed, when "richardson.green", \eqn{u}{u}
        vectors are transformed such that the components of resulting \eqn{v}{v}
	vectors follow independently beta distributions with parameters
	given further by \code{p = prop.revjump$transform.split.combine.parms}
	such that in the triplet of v's used in a particular split move,
        \eqn{v_1 \sim beta(p_1, p_2), v_2 \sim beta(p_3, p_4), v_3 \sim beta(p_5, p_6)}{ %
          v[1] ~ beta(p[1], p[2]), v[2] ~ beta(p[3], p[4]), v[3] ~ beta(p[5], p[6])}.
        When "brooks" \eqn{v_2}{v[2]} is further transformed by
	\eqn{|2v_2 - 1|}{|2*v[2] - 1|}. Default values of
	\tabular{c}{
	\code{prop.revjump$transform.split.combine$parms} 
	}
	is \code{c(2, 2, 2, 2, 1, 1)}.}
      \item{transform.split.combine.parms}{see above.}
      \item{transform.birth.death}{a description of how the canonical
	variables \eqn{u}{u} are to be transformed to give new values of
        mixture component's weight, mean and variance when a birth move
	is proposed. At this moment only one value is possible:
	"richardson.green" implementing the proposal as in Richardson
	and Green (1997).}
    }
  }      
  \item{init}{a list of the initial values to start the McMC. Set to
    \code{NULL} such parameters that you want the program should itself sample
    for you or parameters that are not needed in your model.
    \describe{
      \item{iter}{index of the iteration to which initial values
	correspond, usually zero.}
      \item{mixture}{initial mixture for the error random variable
	\eqn{\varepsilon}{epsilon}. It must a vector of length \code{1 +
        3*kmax}, where \code{mixture[1]} gives initial number of mixture
        of components \eqn{k}{k},
        \code{mixture[2:(k+1)]} gives initial mixture weights,
        \code{mixture[(2+kmax):(2+kmax+k-1)]} gives initial mixture means, 
        \code{mixture[(2+2*kmax):(2+2*kmax+k-1)]} gives initial mixture
	variances. Remaining components of this vector are ignored.}
      \item{beta}{initial values of regression parameters in the same
	order as columns of the design matrix \code{X}. Call the
	function \code{bayessurvreg1} with \code{onlyX = TRUE} to see
	how the columns are sorted. Remember, \code{beta} in this
	function contains both fixed effects \eqn{\beta}{beta} and
	means of random effect \eqn{\gamma}{gamma} in the notation of
	the accompanying paper except the mean of
	the random intercept which is always zero.}
      \item{b}{initial values of random effects \eqn{b_i}{b[i]} for each
        cluster. This must a matrix of size \eqn{q \times N}{q x N} or a
	vector of length \eqn{q*N}{q*N},
	where \eqn{q}{q} is a number of random effects and \eqn{N}{N}
	number of clusters, one column per cluster.}
      \item{D}{initial value for the covariance matrix of random effects
        \eqn{D}{D}. Only its lower triangle must be given in a
	vector, e.g. \code{c(d[1,1], d[2,1], d[3,1], d[2,2], d[3,2],
          d[3,3])} for a matrix \eqn{3 \times 3}{3 x 3}.}
      \item{y}{initial values of true log-event times. This must be a
	vector of length \eqn{\sum_{i=1}^N n_i}{sum[i=1][N] n[i]}.}
      \item{r}{initial values of component labels
	\eqn{r_{i,l}}{r_{i,l}}. This must be a vector of length
	\eqn{\sum_{i=1}^N n_i}{sum[i=1][N] n[i]}.}
      \item{otherp}{initial values for other parameters. At this moment,
        only a value of the parameter \eqn{\eta}{eta} is given here.}
      \item{u}{initial canonical proposal vector of length
	\eqn{3k_{max}}{3*kmax}. When initial number of compoents
	given by \code{init$mixture[1]} is \eqn{k}{k}, effectively only
	last \eqn{3k_{max} - 3*k}{3*kmax - 3*k} components of the
        initial \eqn{u}{u} vector are used. Further, when
	\code{prop.revjump$algorithm = "correlated.av"}, the first
	component of \code{init$u} (\code{init$u[1]}) contains an
	initial mood parameter (\eqn{C_0}{c[0]} in Brooks et al. (2003))
        for the moody ring.}
    }      
  }
  \item{store}{a list that defines which sampled values besides
    regression parameters \eqn{\beta}{beta}, means of random effects
    \eqn{\gamma}{gamma} (both stored in a file called \code{beta.sim}),
    a covariance matrix of random effects \eqn{D}{D} (stored
    in a file \code{D.sim}),
    the mixture (stored in file \code{mixmoment.sim, mweight.sim,
      mmean.sim, mvariance.sim}),
    values of other parameters - \eqn{\eta}{eta} (stored in a file \code{otherp.sim}),
    values of log-likelihoods (stored in a file \code{loglik.sim}),
    information concerning the performance of the reversible jump McMC
    and acceptance of regression parameters (stored in a file \code{MHinfo.sim}),
    iteration indeces (stored in a file \code{iteration.sim})
    are to be stored. The list \code{store} has the following
    components.
    \describe{
      \item{y}{if \code{TRUE} sampled true log-event times are  stored.}
      \item{r}{if \code{TRUE} sampled component labels are stored.}
      \item{b}{if \code{TRUE} sampled values of random effects
	\eqn{b_i}{b[i]} are stored.}
      \item{u}{if \code{TRUE} sampled values of canonical proposal
	vectors for the reversible jump McMC are stored.}
      \item{MHb}{if \code{TRUE} information concerning the performance
	of the Metropolis-Hastings algorithm for the update of random
	effects (if used instead of a dafault Gibbs) is stored.}
      \item{regresres}{if \code{TRUE} sampled values of regression
	residuals at each iteration are stored. The regression residual
	is defined as \eqn{res_{i,l} = \log(t_{i,l}) - \beta^T x_{i,l} -
	b_i^T z_{i,l}}{res[i,l] = log(t[i,l]) - beta'*x[i,l] - b[i]'*z[i,l]}.}
  }
  In the case that either \code{store$y}, or \code{store$r}, or
  \code{store$b}, or \code{store$u} are \code{FALSE}, only the last
  values of either \code{y}, or \code{r}, or \code{b}, or \code{u}
  at the time of writting of remaining quantities are stored in
  appropriate files (without headers) to be possibly used by
  \code{\link{bayessurvreg1.files2init}} function.  
  }
  \item{dir}{a string that specifies a directory where all sampled
    values are to be stored.}
  \item{toler.chol}{tolerance for the Cholesky decomposition.}
  \item{toler.qr}{tolerance for the QR decomposition.}
  \item{\dots}{who knows?}
}
\value{
  A list of class \code{bayessurvreg} containing an information
  concerning the initial values and prior choices.
}
\section{Files created}{  
  Additionally, the following files with sampled values
  are stored in a directory specified by \code{dir} parameter of this
  function (some of them are created only on request, see \code{store}
  parameter of this function).
  \describe{
    \item{iteration.sim}{one column labeled \code{iteration} with
      indeces of McMC iterations to which the stored sampled values correspond.}
    \item{loglik.sim}{two columns labeled \code{loglik} and
      \code{randomloglik}.
      \deqn{\mbox{\code{loglik}} = \sum_{i=1}^{N}\sum_{l=1}^{n_i}\Biggl[
	\biggl\{
        \log\Bigl(\frac{1}{\sqrt{2\pi\sigma_{r_{i,l}}^2}}\Bigr)
	-\frac{(y_{i,l} - \beta^T x_{i,l} - b_i^T z_{i,l} - \mu_{r_{i,l}})^2}{2\sigma_{r_{i,l}}^2}
	\biggr\}
	\Biggr],}{%
        \code{loglik} = sum[i=1][N]sum[l=1][n[i]]
        (log(1/sqrt(2*pi*sigma[r[i,l]]^2)) - 0.5*sigma[r[i,l]]^(-2)*(y[i,l] - beta'*x[i,l] - b[i]*z[i,l]-mu[r[i,l]])^2),}
      where
      \eqn{y_{i,l}}{y[i,l]} denotes (sampled) \emph{(i,l)}th true
      log-event time,
      \eqn{b_i}{b[i]} sampled value of the random effect vector for the
      \emph{i}th cluster,
      \eqn{\beta}{beta} sampled value of the regression parameter
      \eqn{\beta}{beta} and
      \eqn{k, w_j, \mu_j, \sigma_j^2, j = 1,\dots,k}{k, w[j], mu[j], sigma[j]^2, j =
	1,\dots, k} sampled mixture at each iteration.

      \deqn{\mbox{\code{randomloglik}} =
	\sum_{i=1}^{N}\log\Bigl(g(b_i)\Bigr),}{%
	\code{randomloglik} = sum[i=1][N] log(g(b[i])),}
      where \eqn{g}{g} denotes a density of
      (multivariate) normal distribution
      \eqn{N(\gamma, D),}{N(gamma, D),} where
      \eqn{\gamma}{gamma} is a sampled value of the mean of random
      effect vector and \eqn{D}{D} is a sampled value of the covariance
      matrix of the random effects at each iteration.      
    }
    \item{mixmoment.sim}{three columns labeled \code{k}, \code{Intercept} and
      \code{Scale}. These are the number of mixture components, mean and
      standard deviation of the sampled error distribution (mixture) at
      each iteration.}
    \item{mweight.sim}{each row contains mixture weights
      \eqn{w_1,\dots,w_k}{w[1],...,w[k]}
      at each iteration. From the header of this file, maximal number of mixture
      components specified in the prior can be derived.}
    \item{mmean.sim}{each row contains mixture means
      \eqn{\mu_1,\dots,\mu_k}{mu[1],...,mu[k]}
      at each iteration. From the header of this file, maximal number of mixture
      components specified in the prior can be derived.}
    \item{mvariance.sim}{each row contains mixture variances
      \eqn{\sigma^2_1,\dots,\sigma^2_k}{sigma2[1],...,sigma2[k]} at each
      iteration. From the header of this file, maximal number of mixture
      components specified in the prior can be derived.}    
    \item{beta.sim}{columns labeled according to name of the design
      matrix. These are sampled values of regression parameters
      \eqn{\beta}{beta} and means of random effects \eqn{\gamma}{gamma}
      (except the mean of the random intercept which is zero).}
    \item{b.sim}{columns labeled \code{nameb[1].id[1], \dots,
	nameb[q].id[1], \dots, nameb[1].id[N], \dots, nameb[q].id[N]},
      where \code{q} is a dimension of the random effect vector
      \eqn{b_i}{b[i]} and \code{N}{N} number of clusters. \code{nameb}
      is replaced by appropriate column name from the design matrix and
      \code{id} is replaced by identificator of the clusters. This gives
      sampled values of the random effects for each cluster.}
    \item{D.sim}{columns labeled \code{det, D.s.t, s = 1,\dots, q, t =
	s,\dots,q}, where \code{q} is dimension of the random effect
      vector \eqn{b_i}{b[i]}. Column \code{det} gives a determinant of
      the covariance matrix \eqn{D}{D} of the random effects at each
      iteration, remaining columns give a lower triangle of this matrix
      at each iteration.}
    \item{Y.sim}{columns labeled \code{Y[m]} where \code{m} goes from 1
      to \eqn{\sum_{i=1}^{N}n_i}{sum[i=1][N]n_i}. This gives sampled
      log-event times for each observation in the dataset at each
      iteration.}
    \item{r.sim}{columns labeled \code{r[m]} where \code{m} goes from 1
      to \eqn{\sum_{i=1}^{N}n_i}{sum[i=1][N]n_i}. This gives sampled
      mixture labels for each observation in the dataset at each iteration.}
    \item{otherp.sim}{Currently only one column labeled \code{eta} that
      gives sampled values of the hyperparameter \eqn{\eta}{eta}.}
    \item{MHinfo.sim}{this gives the information concerning the
      performance of reversible jump algorithm and a sampler of
      regression parameters \eqn{\beta}{beta} and means of random
      effects \eqn{\gamma}{gamma}. It has columns
      \describe{
        \item{\code{accept.spl.comb}}{relative frequency of accepted
	  split-combine moves up to that iteration.}
	\item{\code{split}}{relative frequency of proposed split moves
	  up to that iteration.}
	\item{\code{accept.birth.death}}{relative frequency of accepted
	  birth-death moves up to that iteration.}
	\item{\code{birth}}{relative frequency of proposed birth moves
	  up to that iteration.}
	\item{\code{beta.block.m}}{with \code{m} going from 1 to number
	  of defined blocks of beta parameters. This gives a relative
	  frequency of accepted proposals for each block up to that
	  iteration. When Gibbs move is used, these should be columns of
	  ones.}
      }	
    }
    \item{MHbinfo.sim}{this gives the information concerning the
      performance of a sampler for random effects (relative frequency of
      accepted values for each cluster and each block of random effects
      updated together). When Gibbs move is used only ones are seen in
      this file.}
    \item{u.sim}{Sampled values of canonical proposal variables for
      reversible jump algorithm are stored here. This file is useful
      only when trying to restart the simulation from some specific point.}
    \item{regresres.sim}{columns labeled \code{res[m]} where \code{m} goes from 1
      to \eqn{\sum_{i=1}^{N}n_i}{sum[i=1][N]n_i.} This stores so called
      regression residuals for each observation at each iteration. This
      residual is defined as
      \deqn{res_{i,l} = y_{i,l} - \beta^T x_{i,l} - b_i z_{i,l},\qquad
	i=1\dots,N,\quad l=1,\dots,n_i,}{%
	res[i,l] = y[i,l] - beta'*x[i,l] - b[i]'*z[i,l], i=1,\dots,N,
	l=1,\dots,n[i],} where \eqn{y_{i,l}}{y[i,l]} is a (sampled)
        log-event time at each iteration.
    }
  }        
}
\references{
\eqn{\mbox{Kom\'{a}rek, A.}}{Komarek, A.} (2006).
\emph{Accelerated Failure Time Models for Multivariate
  Interval-Censored Data with Flexible Distributional Assumptions}.
  PhD. Thesis, Katholieke Universiteit Leuven, Faculteit Wetenschappen.
  
\eqn{\mbox{Kom\'arek}}{Komarek}, A., Lesaffre, E. (2007).
Bayesian accelerated failure time model for correlated interval-censored data
with a normal mixture as an error distribution.
\emph{Statistica Sinica}, \bold{17}, 549--569.
  
Brooks, S. P., Giudici, P., and Roberts, G. O. (2003).
Efficient construction of reversible jump Markov chain Monte Carlo
proposal distribution (with Discussion).
\emph{Journal of the Royal Statistical Society B,} \bold{65,} 3 - 55.

Green, P. J. (1995).
Reversible jump MCMC computation and Bayesian model determination.
\emph{Biometrika,} \bold{82,} 711-732.

Richardson, S., and Green, P. J. (1997).
On Bayesian analysis of mixtures with unknown number of components (with
Discussion).
\emph{Journal of the Royal Statistical Society B,} \bold{59,} 731 - 792.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\examples{
## See the description of R commands for
## the models described in
## Komarek (2006),
## Komarek and Lesaffre (2007).
## 
## R commands available
## in the documentation
## directory of this package as
## cgd.pdf, cgd.R
## tandmobMixture.pdf, tandmobMixture.R
}
\keyword{survival}
\keyword{regression}
\keyword{multivariate}
