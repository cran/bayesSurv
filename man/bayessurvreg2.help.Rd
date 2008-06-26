\name{bayessurvreg2.help}
\alias{bayessurvreg2.help}
\alias{bayessurvreg2.checkStore}
\alias{bayessurvreg2.priorInit}
\alias{bayessurvreg2.priorBeta}
\alias{bayessurvreg2.priorb}
\alias{bayessurvreg2.writeHeaders}
\title{
  Helping functions for Bayesian regression with an error distribution
  smoothed using G-splines
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of `bayessurvreg2' function to make it more
  readable for the programmer.
}
\usage{
bayessurvreg2.checkStore(store)

bayessurvreg2.priorInit(prior, init, design, mcmc.par,
                        prior2, init2, design2, mcmc.par2, doubly)

bayessurvreg2.priorBeta(prior.beta, init, design)

bayessurvreg2.priorb(prior.b, init, design)

bayessurvreg2.writeHeaders(dir, doubly, prior.init, store, design, design2)
}
\arguments{
  \item{store}{a~list as required by the argument \code{store} of the
    function \code{\link{bayessurvreg2}}}
  \item{prior}{a~list as required by the argument \code{prior} of the
    function \code{\link{bayessurvreg2}}}
  \item{prior2}{a~list as required by the argument \code{prior2} of the
    function \code{\link{bayessurvreg2}}}
  \item{init}{a~list as required by the argument \code{init} of the
    function \code{\link{bayessurvreg2}}}
  \item{init2}{a~list as required by the argument \code{init2} of the
    function \code{\link{bayessurvreg2}}}
  \item{mcmc.par}{a~list as required by the argument \code{mcmc.par} of the
    function \code{\link{bayessurvreg2}}}
  \item{mcmc.par2}{a~list as required by the argument \code{mcmc.par2} of the
    function \code{\link{bayessurvreg2}}}
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
    \code{prior.beta2} of the function \code{\link{bayessurvreg2}}}
  \item{prior.b}{a~list as required by the argument \code{prior.b} or
    \code{prior.b2} of the function \code{\link{bayessurvreg2}}}  
  \item{init}{a~list as required by the argument \code{init} or \code{init2}
    of the function \code{\link{bayessurvreg2}}}
  \item{dir}{path to the directory where the sampled values are to be
    stored}
  \item{prior.init}{a~list as returned by the function
    \code{\link{bayessurvreg2.priorInit}}}  
}
\value{
Some lists.
}
\section{Value for bayessurvreg2.priorBeta}{
  A list with the following components:
  \describe{
    \item{parmI}{integer arguments for C++ \code{BetaGamma}
      constructor}
    \item{parmD}{double arguments for C++ \code{BetaGamma}
      constructor}
  }
  and the following attributes:
  \tabular{l}{
    \bold{init} \cr
    \bold{prior.beta} \cr    
  }      
}
\section{Value for bayessurvreg2.priorInit}{
The same object as that returned by \code{\link{bayesBisurvreg.priorInit}}.
}  
\section{Value for bayessurvreg2.priorb}{
  A list with the following components:
  \describe{
    \item{bparmI}{integer arguments for C++ \code{RandomEff}
      constructor}
    \item{bparmD}{double arguments for C++ \code{RandomEff}
      constructor}
    \item{DparmI}{integer arguments for C++ \code{CovMatrix}
      constructor}
    \item{DparmD}{double arguments for C++ \code{CovMatrix}
      constructor}    
  }
  and the following attributes:
  \tabular{l}{
    \bold{init.b} \cr
    \bold{init.D} \cr    
    \bold{prior.b} \cr    
  }      
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}
\keyword{internal}
