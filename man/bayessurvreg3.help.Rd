\name{bayessurvreg3.help}
\alias{bayessurvreg3.help}
\alias{bayessurvreg3.checkStore}
\alias{bayessurvreg3.priorInit}
\alias{bayessurvreg3.priorBeta}
\alias{bayessurvreg3.priorb}
\alias{bayessurvreg3.writeHeaders}
\title{
  Helping functions for Bayesian regression with an error distribution
  smoothed using G-splines
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of `bayessurvreg3' function to make it more
  readable for the programmer.
}
\usage{
bayessurvreg3.checkStore(store)

bayessurvreg3.priorInit(prior, init, design, mcmc.par,
   prior2, init2, design2, mcmc.par2, doubly)

bayessurvreg3.priorBeta(prior.beta, init, design)

bayessurvreg3.priorb(prior.b, init, design, mcmc.par)

bayessurvreg3.writeHeaders(dir, doubly, prior.init,
   priorb.di, priorb2.di, store, design, design2)
}
\arguments{
  \item{store}{a~list as required by the argument \code{store} of the
    function \code{\link{bayessurvreg2}}}
  \item{prior}{a~list as required by the argument \code{prior} of the
    function \code{\link{bayessurvreg3}}}
  \item{prior2}{a~list as required by the argument \code{prior2} of the
    function \code{\link{bayessurvreg3}}}
  \item{init}{a~list as required by the argument \code{init} of the
    function \code{\link{bayessurvreg3}}}
  \item{init2}{a~list as required by the argument \code{init2} of the
    function \code{\link{bayessurvreg3}}}
  \item{mcmc.par}{a~list as required by the argument \code{mcmc.par} of the
    function \code{\link{bayessurvreg3}}}
  \item{mcmc.par2}{a~list as required by the argument \code{mcmc.par2} of the
    function \code{\link{bayessurvreg3}}}
  \item{design}{an~object as returned by the function
    \code{\link{bayessurvreg.design}} related to either the onset time
    if doubly censored observations or to the event time. Remark:
    \code{design$Y} contains original times and NOT their logarithmic transformations.}
  \item{design2}{an~object as returned by the function
    \code{\link{bayessurvreg.design}} related to time-to-event
    if doubly censored observations. Remark:
    \code{design2$Y} contains original times and NOT their logarithmic transformations.}    
  \item{doubly}{logical indicating whether the response is doubly
    censored or not}  
  \item{prior.beta}{a~list as required by the argument \code{prior.beta} or
    \code{prior.beta2} of the function \code{\link{bayessurvreg3}}}
  \item{prior.b}{a~list as required by the argument \code{prior.b} or
    \code{prior.b2} of the function \code{\link{bayessurvreg3}}}  
  \item{init}{a~list as required by the argument \code{init} or \code{init2}
    of the function \code{\link{bayessurvreg3}}}
  \item{dir}{path to the directory where the sampled values are to be
    stored}
  \item{prior.init}{a~list as returned by the function
    \code{\link{bayessurvreg3.priorInit}}}
  \item{priorb.di}{a~list as returned by the function
    \code{\link{bayessurvreg3.priorb}}}
  \item{priorb2.di}{a~list as returned by the function
    \code{\link{bayessurvreg3.priorb}}}  
}
\value{
Some lists.
}
\section{Value for bayessurvreg3.priorb}{
  A list with the following components:
  \describe{
    \item{bparmI}{integer arguments for C++ \code{RandomEff}
      constructor}
    \item{bparmD}{double arguments for C++ \code{RandomEff}
      constructor}
    \item{GsplI}{integer arguments for C++ \code{Gspline}
      constructor related to the smothed density of the random intercept}
    \item{GsplD}{double arguments for C++ \code{Gspline} constructor
      related to the smoothed density of the random intercept}
    \item{specification}{1 or 2, one of the G-spline specifications
      related to the distribution of the random intercept}
    \item{r}{initial component labels (vector of size \eqn{ncluster}{ncluster}) taking
      values from 1 to the total length of the G-spline related to the
      random intercept}
  }
  and the following attributes:
  \tabular{l}{
    \bold{prior.b} \cr
    \bold{init} \cr    
    \bold{mcmc.par} \cr    
  }      
}   
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\keyword{internal}

