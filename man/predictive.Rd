\name{predictive}
\alias{predictive}
\alias{predictive.control}
\title{
  Compute predictive quantities based on a Bayesian survival regression model.
}
\description{
  This function runs additional McMC to compute predictive survivor and
  hazard curves and predictive event times for specified values of
  covariates.

  Firstly, the function \code{\link{bayessurvreg1}} has to be used to
  obtain a sample from the posterior distribution of unknown quantities.

  Directly, posterior predictive quantiles and means of asked quantities
  are computed and stored in files.
  
  Function \code{predictive.control} serves only to perform some input
  checks inside the main function \code{predictive}.
}
\usage{
predictive(
     formula,
     random,
     time0 = 0,
     data = parent.frame(),
     grid,
     type,
     subset,
     na.action = na.fail,
     quantile = c(0, 0.025, 0.5, 0.975, 1),                       
     nsimul = list(niter = 10, nwrite = 10),
     predict = list(Et=TRUE, t=FALSE, Surv=TRUE,
                    hazard=FALSE, cum.hazard=FALSE),
     store = list(Et=TRUE, t = FALSE, Surv = FALSE,
                  hazard = FALSE, cum.hazard=FALSE),
     Eb0.depend.mix = FALSE,
     dir = getwd(),
     toler.chol = 1e-10,
     toler.qr = 1e-10,
     \dots)

predictive.control(predict, store, quantile)
}
\arguments{
  \item{formula}{same formula as that one used to sample from the
    posterior distribution of unknown quantities by the function
    \code{\link{bayessurvreg1}}.}
  \item{random}{same \code{random} statement as that one used to sample from the
    posterior distribution of unknown quantities by the function
    \code{\link{bayessurvreg1}}.}
  \item{time0}{starting time for the survival model. This option is used
  to get correct hazard function in the case that the original model was
  \eqn{\log(T-time0) = \dots}{log(T - time0) = \dots}.}
  \item{data}{optional data frame in which to interpret the variables
    occuring in the formulas. Usually, you create a new
    \code{data.frame} similar to that one used to obtain a sample from
    the posterior distribution. In this new \code{data.frame}, put
    covariate values equal to these for which predictive quantities are
    to be obtained. If \code{cluster} statement was used, assign a
    unique cluster identification to each observation. Response variable
    and a censoring indicator may be set to arbitrary values. They are
    only used in \code{formula} but are ignored for computation.}
  \item{grid}{a list of length as number of observations in \code{data} or a vector
    giving grids of values where predictive survivor functions, hazards, cumulative
    hazards are to be evaluated. If it is a vector, same grid is used for all
    observations from \code{data}. Not needed if only \code{predict$t}
    or \code{predict$Et} are \code{TRUE}. If \code{time0} is different
    from zero your grid should start at \code{time0} and not at zero.}
  \item{type}{a character string giving the type of assumed error distribution. Currently, valid
    are substrings of "mixture". In the future, "spline",
    "polya.tree" might be also implemented.}
  \item{subset}{subset of the observations from the \code{data} to be
    used. This option will normally not be needed.}  
  \item{na.action}{function to be used to handle any \code{NA}s in the
    data. The user is discouraged to change a default value
    \code{na.fail}.}
  \item{quantile}{a vector of quantiles that are to be computed for each
    predictive quantity.}
  \item{nsimul}{a list giving the length of the simulation used to
    sample from posterior predictive distribution. It should be
    consistent with already simulated values obtained by
    \code{\link{bayessurvreg1}} function. The list has the following
    components.
    \describe{
      \item{niter}{number of iterations that are to be performed. It
	should not be higher than the number already simulated values.}
      \item{nwrite}{interval in which predictive quantities are written
	to files.}
    }      
  }
  \item{predict}{a list of logical values indicating which predictive quantities are to be sampled.
    Components of the list:
    \describe{
      \item{Et}{predictive expectations of survivor times}
      \item{t}{predictive survivor times}
      \item{Surv}{predictive survivor functions}
      \item{hazard}{predictive hazard functions}
      \item{cum.hazard}{predictive cumulative hazard functions}
    }
  }  
  \item{store}{a list of logical values indicating which predictive
    quantities are to be stored in files as `predET*.sim', `predT*.sim',
    `predS*.sim', `predhazard*.sim', `predcumhazard*.sim'. If you are
    interested only in posterior means or quantiles of the predictive
    quantities you do not have to store sampled values. Posterior means
    and quantiles are stored in files `quantET*.sim', `quantT*.sim',
    `quantS*.sim', `quanthazard*.sim', `quantpredhazard*.sim'.}
  \item{Eb0.depend.mix}{a logical value indicating whether the mean of
    the random intercept (if included in the model) was given in a
    hierarchical model as an overall mean of the mixture in the error
    term. With \code{FALSE} (default) you have the same model as that
    one described in an accompanying paper. An ordinary user is
    discouraged from setting this to \code{TRUE}.}
  \item{dir}{a string giving a directory where previously simulated
    values were stored and where newly obtained quantities will be
    stored. On Unix, do not use `~/' to specify your home directory.
    A full path must be given, e.g. `/home/arnost/'.}
  \item{toler.chol}{tolerance for the Cholesky decomposition.}
  \item{toler.qr}{tolerance for the QR decomposition.}  
  \item{\dots}{who knows?}
}
\value{
  An integer which should be equal to zero if everything ran fine.
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\examples{
  ## See attached files.
}
\keyword{survival}
\keyword{regression}
