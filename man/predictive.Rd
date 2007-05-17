\name{predictive}
\alias{predictive}
\alias{predictive.control}
\title{
  Compute predictive quantities based on a Bayesian survival regression
  model fitted using bayessurvreg1 function.  
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
predictive(formula, random, time0 = 0, data = parent.frame(),
     grid, type = "mixture", subset, na.action = na.fail,
     quantile = c(0, 0.025, 0.5, 0.975, 1),                       
     skip = 0, by = 1, last.iter, nwrite, only.aver = FALSE,
     predict = list(Et=TRUE, t=FALSE, Surv=TRUE, hazard=FALSE, cum.hazard=FALSE),
     store = list(Et=TRUE, t = FALSE, Surv = FALSE, hazard = FALSE, cum.hazard=FALSE),
     Eb0.depend.mix = FALSE,
     dir = getwd(), toler.chol = 1e-10, toler.qr = 1e-10)

predictive.control(predict, store, only.aver, quantile)
}
\arguments{
  \item{formula}{the same formula as that one used to sample from the
    posterior distribution of unknown quantities by the function
    \code{\link{bayessurvreg1}}.}
  \item{random}{the same \code{random} statement as that one used to sample from the
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
  \item{skip}{number of rows that should be skipped at the beginning of
    each *.sim file with the stored sample.}
  \item{by}{additional thinning of the sample.}
  \item{last.iter}{index of the last row from *.sim files that should be
    used. If not specified than it is set to the maximum available
    determined according to the file \code{mixmoment.sim}.}
  \item{nwrite}{frequency with which is the user informed about the
    progress of computation (every \code{nwrite}th iteration count of
    iterations change).}
  \item{only.aver}{if \code{TRUE} only posterior predictive mean is
    computed for all quantities and no quantiles.}
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
}
\value{
  An integer which should be equal to zero if everything ran fine.
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
