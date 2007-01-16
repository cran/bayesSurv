\name{predictive2}
\alias{predictive2}
\alias{predictive2.control}
\title{
  Compute predictive quantities based on a Bayesian survival regression
  model fitted using bayesBisurvreg or bayessurvreg2 or bayessurvreg3 functions.
}
\description{
  This function computes predictive densities, survivor and hazard
  curves for specified combinations of covariates.

  Firstly, either the function \code{\link{bayesBisurvreg}} or the
  function \code{\link{bayessurvreg2}} or the function \code{\link{bayessurvreg3}}
  has to be used to obtain a sample from the posterior distribution of unknown quantities.

  Function \code{predictive2.control} serves only to perform some input
  checks inside the main function \code{predictive2}.
}
\usage{
predictive2(formula, random, obs.dim, time0, data = parent.frame(),
     grid, na.action = na.fail, Gspline,
     quantile = c(0, 0.025, 0.5, 0.975, 1),
     skip = 0, by = 1, last.iter, nwrite,
     only.aver = TRUE,
     predict = list(density=FALSE, Surv=TRUE,
                    hazard=FALSE, cum.hazard=FALSE),
     dir = getwd(), extens = "", extens.random="_b", version=0)

predictive2.control(predict, only.aver, quantile, obs.dim,
     time0, Gspline, n)
}
\arguments{
\item{formula}{the same formula as that one used to sample from the
  posterior distribution of unknown quantities by the function
  \code{\link{bayesBisurvreg}} or \code{\link{bayessurvreg2}} or
  \code{\link{bayessurvreg3}}. On the left hand side whichever
  \code{Surv} object of a~proper length
  can be used (it is ignored anyway).

  REMARK: the prediction must be asked for at least two combinations of
  covariates. This is the restriction imposed by one of the internal
  functions I use.
}  %% end of item{formula}
\item{random}{the same \code{random } statement as that one used to sample from the
  posterior distribution of unknown quantities by the function
  \code{\link{bayessurvreg2}} or \code{\link{bayessurvreg3}}, not applicable if
  \code{\link{bayesBisurvreg}} was used to sample from the posterior
  distribution.
}  %% end of item{random}
\item{obs.dim}{a vector that has to be supplied if bivariate data were
  used for estimation (using the function
  \code{\link{bayesBisurvreg}}). This vector has to be of the same
  length as the number of covariate combinations for which the
  predictive quantities are to be computed. It determines to which
  dimension (1 or 2) each observation belong.
}  %% end of item{obs.dim}
\item{time0}{a~vector of length \code{Gspline$dim} giving the starting
  time for the survival model. It does not have to be supplied if equal
  to zero (usually).
  This option is used to get hazard and density functions on the
  original time scale in the case that the model was
  \eqn{\log(T-time0) = \dots}{log(T - time0) = \dots}. Note that
  \code{time0} IS NOT the starting time of doubly censored observation
  since there after subtracting the onset time, the starting time is
  (usually) equal to zero.
}  %% end of item{time0}
\item{data}{optional data frame in which to interpret the variables
  occuring in the formulas. Usually, you create a new
  \code{data.frame} similar to that one used to obtain a sample from
  the posterior distribution. In this new \code{data.frame}, put
  covariate values equal to these for which predictive quantities are
  to be obtained. If \code{cluster} statement was used, assign a
  unique cluster identification to each observation. Response variable
  and a censoring indicator may be set to arbitrary values. They are
  only used in \code{formula} but are ignored for computation.
}  %% end of item{data}
\item{grid}{a~vector giving the grid of values where predictive
  quantities are to be evaluated. The grid should normally start at some
  value slightly higher than \code{time0}.
}  %% end of item{grid}
\item{na.action}{function to be used to handle any \code{NA}s in the
  data. The user is discouraged to change a default value
  \code{na.fail}.
}  %% end of item{na.action}
\item{Gspline}{a~list specifying the G-spline used for the error
  distribution in the model. It is a~list with the following components:
  \describe{
    \item{dim}{dimension of the G-spline, in the case when the function
      \code{\link{bayesBisurvreg}} was used to fit the model this will
      usually be equal to 2, in the case when the function
      \code{\link{bayessurvreg2}} was used to fit the model this MUST be
      equal to 1.}
    \item{K}{a~vector of length \code{Gspline$dim} specifying the number
      of knots at each side of the middle knot for each dimension of the
      G-spline.}
  }
}  %% end of item{Gspline}  
\item{quantile}{a vector of quantiles that are to be computed for each
  predictive quantity.
}  %% end of item{quantile}
\item{skip}{number of rows that should be skipped at the beginning of
  each *.sim file with the stored sample.
}  %% end of item{skip}
\item{by}{additional thinning of the sample.
}  %% end of item{by}
\item{last.iter}{index of the last row from *.sim files that should be
  used. If not specified than it is set to the maximum available
  determined according to the file \code{mixmoment.sim}.
}  %% end of item{last.iter}
\item{nwrite}{frequency with which is the user informed about the
  progress of computation (every \code{nwrite}th iteration count of
  iterations change).
}  %% end of item{nwrite}
\item{only.aver}{if \code{TRUE} only posterior predictive mean is
  computed for all quantities and no quantiles.

  The word of warning: with \code{only.aver} set to \code{FALSE}, all
  quantities must be stored for all iterations of the MCMC to be able to
  compute the quantiles. This might require quite lots of memory.
}  %% end of item{only.aver}
\item{predict}{a list of logical values indicating which predictive quantities are to be computed.
  Components of the list:
  \describe{
    \item{density}{predictive density}
    \item{Surv}{predictive survivor functions}
    \item{hazard}{predictive hazard functions}
    \item{cum.hazard}{predictive cumulative hazard functions}
  }
}  %% end of item{predict}
\item{dir}{directory where to search for files (`mixmoment.sim',
  `mweight.sim', mmean.sim', gspline.sim', 'beta.sim', 'D.sim', \dots) with the McMC
  sample.
}  %% end of item{dir}
\item{extens}{an extension used to distinguish different sampled
  G-splines if more formulas were used in one MCMC simulation (e.g. with
  doubly-censored data).

  \itemize{
    \item if the data were not doubly-censored or if you wish to compute
    predictive quantities for the \emph{onset} time of the
    doubly-censored data then
    \tabular{c}{
      \code{extens = ""}
    }

    \item if the data were doubly-censored and you wish to compute
    predictive quantities for the \emph{event} time then
    \tabular{c}{
      \code{extens = "_2"}
    }
  }  
}  %% end of item{extens}
\item{extens.random}{only applicable if the function
  \code{\link{bayessurvreg3}} was used to generate the MCMC sample.

  This is an extension used to distinguish different sampled G-splines
  determining the distribution of the random intercept (under the
  presence of doubly-censored data).

  \itemize{
    \item if the data were not doubly-censored or if you wish to compute
    predictive quantities for the \emph{onset} time of the
    doubly-censored data then
    \tabular{c}{
      \code{extens.random = "_b"}
    }

    \item if the data were doubly-censored and you wish to compute
    predictive quantities for the \emph{event} time then
    \tabular{c}{
      \code{extens.random = "_b2"}
    }
  }  
}  %% end of item{extens.random}
\item{version}{ this argument indicates by which \code{bayes*survreg*} function the
  chains used by \code{bayesGspline} were created. Use the following:

  \describe{
    \item{bayesBisurvreg:}{\code{version = 0};}
    \item{bayessurvreg2:}{\code{version = 0};}
    \item{bayessurvreg3:}{with all distributions specified as G-splines:
      \code{version = 3};}
    \item{bayessurvreg3:}{with error distributions specified as
      G-splines and bivariate normal random intercepts: \code{version = 32}.}    
  }
}  %% end of item{version}
\item{n}{number of covariate combinations for which the prediction will
  be performed.
}  %% end of item{n}  
}  %% end of arguments
\value{
A list with possibly the following components (what is included depends
on the value of the arguments \code{predict} and \code{only.aver}):

  \item{grid}{a~vector with the grid values at which the survivor
    function, survivor density, hazard and cumulative hazard are computed.}
  \item{Surv}{predictive survivor functions.

    A~matrix with as many columns as length(grid) and as many
    rows as the number of covariate combinations for which the
    predictive quantities were asked. One row per covariate combination.
  } 
  \item{density}{predictive survivor densities.

    The same structure as \code{Surv} component of the list.
  }
  \item{hazard}{predictive hazard functions.

    The same structure as \code{Surv} component of the list.    
  }
  \item{cum.hazard}{predictive cumulative hazard functions.

    The same structure as \code{Surv} component of the list.    
  }
  \item{quant.Surv}{pointwise quantiles for the predictive survivor
    functions.

    This is a list with as many components as the number of covariate
    combinations. One component per covariate combination.

    Each component of this list is a~matrix with as many columns as
    length(grid) and as many rows as the length of the argument
    \code{quantile}. Each row of this matrix gives values of one
    quantile. The rows are also labeled by the probabilities (in \%) of
    the quantiles.
  }
  \item{quant.density}{pointwise quantiles for the predictive survivor
    densities.

    The same structure as \code{quant.Surv} component of the list.
  }
  \item{quant.hazard}{pointwise quantiles for the predictive hazard
    functions.

    The same structure as \code{quant.Surv} component of the list.    
  }
  \item{quant.cum.hazard}{pointwise quantiles for the predictive
    cumulative hazard functions.

    The same structure as \code{quant.Surv} component of the list.
  }
}
\references{
  \eqn{\mbox{Kom\'{a}rek, A.}}{Komarek, A.} (2006).
  \emph{Accelerated Failure Time Models for Multivariate
    Interval-Censored Data with Flexible Distributional Assumptions}.
  PhD. Thesis, Katholieke Universiteit Leuven, Faculteit Wetenschappen.

  \eqn{\mbox{Kom\'{a}rek}}{Komarek}, A. and Lesaffre, E. (2006).
  Bayesian accelerated failure time model with multivariate doubly-interval-censored data
  and flexible distributional assumptions.
  \emph{Submitted.} \cr
  See \code{Komarek_Lesaffre_2006.pdf}.

  \eqn{\mbox{Kom\'{a}rek}}{Komarek}, A. and Lesaffre, E. (2006b).
  Bayesian semi-parametric accelerated failurew time model for paired
  doubly interval-censored data.
  \emph{Statistical Modelling}, \bold{6}, 3--22.
  
  \eqn{\mbox{Kom\'{a}rek}}{Komarek}, A. Lesaffre, E., and Legrand, C. (2006).
  EORTC data article.
  \emph{To be written.}
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\examples{
## See the description of R commands for
## the models described in
## Komarek (2006),
## Komarek and Lesaffre (2006),
## Komarek and Lesaffre (2006b),
## Komarek, Lesaffre, and Legrand (2006).
##
## R commands available in the documentation
## directory of this package
## as tandmobCS.pdf, tandmobCS.R,
##    tandmobPA.pdf, tandmobPA.R.
}
\keyword{survival}
\keyword{regression}
