\name{bayessurvreg1.help}
\alias{bayessurvreg1.help}
\alias{bayessurvreg1.checkStore}
\alias{bayessurvreg1.priorInit}
\alias{bayessurvreg1.revjump}
\alias{bayessurvreg1.writeHeaders}
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

bayessurvreg1.revjump(prop.revjump)

bayessurvreg1.writeHeaders(dir, prior, store, nX, X, names.random,
   ncluster, nrandom, rnamesX, unique.cluster, nBetaBlocks, nbBlocks)
}
\arguments{
  \item{store}{}
  \item{prior, init, Yinit, Xinit, n, nX, nrandom, ncluster}{}
  \item{indb, randomInt, toler.chol}{}
  \item{prop.revjump}{}
  \item{dir, X, names.random}{}
  \item{rnamesX, unique.cluster, nBetaBlocks, nbBlocks}{}
}
\value{
  Some lists.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\keyword{internal}
\keyword{connection}
