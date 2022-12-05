\name{bayesBisurvreg.help}
\alias{bayesBisurvreg.help}
\alias{bayesBisurvreg.checkStore}
\alias{bayesBisurvreg.priorInit}
\alias{bayesBisurvreg.priorBeta}
\alias{bayesBisurvreg.writeHeaders}
\title{
  Helping function for Bayesian regression with smoothed bivariate
  densities as the error term,
  based on possibly censored data
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of `bayesBisurvreg' function to make it more
  readable for the programmer.
}
\usage{
bayesBisurvreg.checkStore(store)

bayesBisurvreg.priorInit(dim, prior, init, design, mcmc.par,
                              prior2, init2, design2, mcmc.par2,
                         doubly)

bayesBisurvreg.priorBeta(prior.beta, init, design)

bayesBisurvreg.writeHeaders(dir, dim, nP, doubly, prior.init, store,
                            design, design2)
}
\arguments{
  \item{store}{a list as required by the argument \code{store} of the
    function \code{\link{bayesBisurvreg}} or \code{\link{bayesBisurvreg.checkStore}}}
  \item{dim}{dimension of the response, 1 or 2}    
  \item{prior}{a list as required by the argument \code{prior} of the
    function \code{\link{bayesBisurvreg}}}
  \item{prior2}{a list as required by the argument \code{prior2} of the
    function \code{\link{bayesBisurvreg}}}
  \item{init}{a list as required by the argument \code{init} of the
    function \code{\link{bayesBisurvreg}}}
  \item{init2}{a list as required by the argument \code{init2} of the
    function \code{\link{bayesBisurvreg}}}
  \item{mcmc.par}{a list as required by the argument \code{mcmc.par} of the
    function \code{\link{bayesBisurvreg}}}
  \item{mcmc.par2}{a list as required by the argument \code{mcmc.par2} of the
    function \code{\link{bayesBisurvreg}}}
  \item{design}{an object as returned by the function
    \code{\link{bayessurvreg.design}} related to either the onset time
    if doubly censored observations or to the event time. Remark:
    \code{design$Y} contains original times and NOT their logarithmic transformations.}
  \item{design2}{an object as returned by the function
    \code{\link{bayessurvreg.design}} related to time-to-event
    if doubly censored observations. Remark:
    \code{design2$Y} contains original times and NOT their logarithmic transformations.}    
  \item{doubly}{logical indicating whether the response is doubly
    censored or not}  
  \item{prior.beta}{a list as required by the argument \code{prior.beta} or
    \code{prior.beta2} of the function \code{\link{bayesBisurvreg}}}
  \item{dir}{path to the directory where the sampled values are to be stored}
  \item{nP}{sample size - number of observations if the univariate model
    is fitted, number of bivariate observational vectors if the
    bivariate model is fitted}
  \item{prior.init}{a list as returned by the function
    \code{\link{bayesBisurvreg.priorInit}}}
}
\value{
  Some lists.
}
\section{Value for bayesBisurvreg.priorInit}{
  A~list with the following components:
  \describe{
  \item{Gparmi}{integer arguments for the G-spline constructor in the
    C++ code related to the onset/event time}
  \item{Gparmd}{double arguments for the G-spline constructor in the
    C++ code related to the onset/event time}
  \item{y}{vector of initial values for the log(onset time)/log(event
    time), sorted as
    \eqn{y_1[1], y_1[2], \dots, y_n[1], y_n[2]}{y[1,1],y[1,2], \dots,
      y[n,1],y[n,2]} in the case of bivariate response with sample size
    equal to \eqn{n}{n}}
  \item{r}{initial component labels (vector of size \eqn{n}{n}) taking
    values from 1 to the total length of the G-spline related to the
    onset/event time}
  \item{Gparmi2}{integer arguments for the G-spline constructor in the
    C++ code related to time-to-event in the case of doubly censoring}
  \item{Gparmd2}{double arguments for the G-spline constructor in the
    C++ code related to time-to-event in the case of doubly censoring}
  \item{y2}{vector of initial values for the time-to-event in the case of doubly censoring
    sorted as
    \deqn{y_1[1], y_1[2], \dots, y_n[1], y_n[2]}{%
      y[1,1],y[1,2], \dots, y[n,1],y[n,2]} in the case of bivariate response with sample size
    equal to \eqn{n}{n}}
  \item{r2}{initial component labels (vector of size \eqn{n}{n}) taking
    values from 1 to the total length of the G-spline related to
    time-to-event in the case of doubly censoring}
  \item{iter}{index of the nullth iteration}
  \item{specification}{2 component vector (one component for onset, one
    for time-to-event), specification of the G-spline model (1 or 2), see
    \code{\link{bayesHistogram}} for more detail}
  \item{y.left}{lower limit of the log-response (or exact/right/left
    censored observation) as required by the C++ function
    \code{bayesBisurvreg}, related to the onset time in the case of
    doubly censoring and to the event time otherwise}
  \item{y.right}{upper limit of the log-response as required by the C++ function
    \code{bayesBisurvreg}, related to the onset time in the case of
    doubly censoring and to the event time otherwise}
  \item{status}{status vector as required by the C++ function
    \code{bayesBisurvreg} related to the onset time in the case of
    doubly censoring and to the event time otherwise}
  \item{t2.left}{lower limit of the response as required by the C++
    function \code{bayesBisurvreg}, related to time-to-event in the
    case of doubly censoring, equal to 0 if there is no doubly-censoring}
  \item{t2.right}{upper limit of the response as required by the C++
    function \code{bayesBisurvreg}, related to time-to-event in the
    case of doubly censoring, equal to 0 if there is no doubly-censoring}
  \item{status2}{status vector related to time-to-event in the case of
    doubly censoring, equal to 0 otherwise.}
  }

  and the following attributes:
  \tabular{l}{
    \bold{init}        \cr
    \bold{prior}       \cr
    \bold{mcmc.par}    \cr
    \bold{init2}       \cr
    \bold{prior2}      \cr
    \bold{mcmc.par2}   \cr
  }    
}
\section{Value for bayesBisurvreg.priorBeta}{
  A list with the following components:
  \describe{
    \item{parmI}{integer arguments for C++ \code{classBetaGamma} constructor}
    \item{parmD}{double arguments for C++ \code{classBetaGamma} constructor}
  }

  and the following attributes:
  \describe{
    \item{init}{a~vector with initial values of the beta parameter,
      equal to \code{numeric(0)} if there are no regressors}
    \item{prior.beta}{a~list with components \code{mean.prior} and
      \code{var.prior} containing vectors with the prior mean and prior
      variance of the \code{beta} parameters}
  }    
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\keyword{internal}
