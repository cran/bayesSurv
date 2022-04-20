\name{bayessurvreg1.help}
\alias{bayessurvreg1.help}
\alias{bayessurvreg1.checkStore}
\alias{bayessurvreg1.priorInit}
\alias{bayessurvreg1.priorBeta}
\alias{bayessurvreg1.priorb}
\alias{bayessurvreg1.writeHeaders}
\alias{bayessurvreg1.revjump}
\title{
   Helping function for Bayesian survival regression models, version 1.
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of `bayessurvreg1' function to make it more
  readable for the programmer.
}
\usage{
bayessurvreg1.checkStore(store)

bayessurvreg1.priorInit(prior, init, Yinit, Xinit, n, nX, nrandom,
   ncluster, indb, randomInt, toler.chol)

bayessurvreg1.priorBeta(prior.beta, nX, indb, factors,
   n.factors, n.in.factors)

bayessurvreg1.priorb(prior.b, nrandom, ncluster, toler.chol)

bayessurvreg1.writeHeaders(dir, prior, store, nX, X, names.random,
   ncluster, nrandom, rnamesX, unique.cluster, nBetaBlocks, nbBlocks)

bayessurvreg1.revjump(prop.revjump)
}
\arguments{
  \item{store}{a~list as required by the \code{store} argument of the
    function \code{\link{bayessurvreg1}}}
  \item{prior, init, Yinit, Xinit, n, nX, nrandom, ncluster}{???}
  \item{indb, randomInt, toler.chol}{???}
  \item{prior.beta, factors, n.factors, n.in.factors}{???}
  \item{prior.b, nrandom, ncluster}{???}  
  \item{prop.revjump}{???}
  \item{dir, X, names.random}{???}
  \item{rnamesX, unique.cluster, nBetaBlocks, nbBlocks}{???}
}
\value{
  Some lists.
}
\section{Value for bayessurvreg1.priorb}{
  A~list with the following components:
  \describe{
  \item{double}{double vector for C++ constructor of class
    \code{bblocks}. It has the following parts:
    \tabular{l}{
      \bold{priordD}       \cr
      \bold{covparLV}      \cr
      \bold{halfRangeUnif} \cr
      \bold{weightUnif}    \cr
    }      
  }
  \item{integer}{integer vector for C++ constructor of class
    \code{bblocks}. It has the following parts:
    \tabular{l}{
        \bold{nrandom}    \cr
        \bold{ncluster}   \cr
        \bold{priorD}     \cr
        \bold{typeUpd}    \cr
        \bold{nBlocks}    \cr
        \bold{nInBlock}   \cr
        \bold{lcovparLV}  \cr
	\bold{indBlockLV} \cr
    }    
  }
  }
  and the following attributes:
  \tabular{l}{
    \bold{prior.b}        \cr
  }
}  
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\keyword{internal}
