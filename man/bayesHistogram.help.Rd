\name{bayesHistogram.help}
\alias{bayesHistogram.help}
\alias{bayesHistogram.design}
\alias{bayesHistogram.checkStore}
\alias{bayesHistogram.priorInit}
\alias{bayesHistogram.writeHeaders}
\title{
  Helping function for Bayesian smoothing of (bi)-variate densities
  based on possibly censored data
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of `bayesHistogram' function to make it more
  readable for the programmer.
}
\usage{
bayesHistogram.design(y1, y2)

bayesHistogram.checkStore(store)

bayesHistogram.priorInit(prior, init, mcmc.par, design)

bayesHistogram.writeHeaders(dir, design, prior.init, store)
}
\arguments{
  \item{y1}{response for the first dimension. This should be a~survival object
    created by \code{\link[survival]{Surv}}}
  \item{y2}{response object for the second dimension (if bivariate
    density is to be smoothed). This should be a~survival object created
    by \code{\link[survival]{Surv}}. Further, it must be consistent with
    \code{y1}, i.e. it has to have the same number of rows.}
  \item{store}{a~list with appropriate components}
  \item{prior}{a~list as required by \code{prior} argument of function \code{\link{bayesHistogram}}}
  \item{init}{a~list as required by \code{init} argument of function \code{\link{bayesHistogram}}}
  \item{mcmc.par}{a~list as required by \code{mcmc.par} argument of
    function \code{\link{bayesHistogram}}}
  \item{design}{a~list with the design information as returned by the
    function \code{bayesHistogram.design}. The list is assumed
    to have the following components:
    \describe{
      \item{y.left}{\dots}
      \item{y.right}{\dots}
      \item{status}{\dots}
      \item{dim}{\dots}
    }
  }  
  \item{dir}{string giving a~directory where to store sampled files}  
  \item{prior.init}{an~object as returned by \code{\link{bayesHistogram.priorInit}}}
}
\value{
  Some lists.
}
\section{Value for bayesHistogram.design}{
  A~list with the following components:
  \describe{
  \item{y.left}{vector or matrix with either observed, right or left
    censored observations or with the lower limits of interval censored
    observations. It is a vector if \code{dim == 1} and it is a matrix with
    2 rows and \eqn{n}{n} columns if \code{dim == 2}, where \eqn{n}{n}
    is a~sample size. In that case, the first row of the matrix are
    responses for the first dimension and the second row of the matrix
    are responses for the second dimension.}
  \item{y.right}{vector or matrix with entries equal to 1 for observed,
    right or left censored observations and entries equal to the upper
    limits of interval censored observations. The structure is the same
    as that of \code{y.left}.}
  \item{status}{a~vector or matrix with censoring indicators (1 =
    exactly observed, 0 = right censored, 2 = left censored, 3 =
    interval censored). The structure is the same as that of \code{y.left}.}
  \item{dim}{dimension of the response, i.e. 1 (univariate smoothing) or
    2 (bivariate smoothing)}
  }
}
\section{Value for bayesHistogram.priorInit}{
  A~list with the following components:
  \describe{
  \item{Gparmi}{integer parameters for the G-spline constructor in the
    C++ code}
  \item{Gparmd}{double parameters for the G-spline constructor in the
    C++ code}
  \item{iter}{index of the nullth iteration}
  \item{y}{vector of initial values for the response, sorted as
    \eqn{y_1[1], y_1[2], \dots, y_n[1], y_n[2]}{y[1,1],y[1,2], \dots,
      y[n,1],y[n,2]} in the case of bivariate response with sample size
    equal to \eqn{n}{n}}
  \item{r}{initial component labels (vector of size \eqn{n}{n}) taking
    values from 1 to the total length of the G-spline}
  \item{specification}{specification of the G-spline model (1 or 2), see
    \code{\link{bayesHistogram}} for more detail}
  }

  and the following attributes:
  \tabular{l}{
  \bold{init}      \cr
  \bold{prior}     \cr
  \bold{mcmc.par}  \cr
  }
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}
\keyword{internal}
