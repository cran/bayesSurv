\name{bayessurvreg.help}
\alias{bayessurvreg.help}
\alias{bayessurvreg.design}
\alias{bayessurvreg.priorBeta}
\alias{bayessurvreg.priorb}
\title{
   Helping function for Bayesian survival regression models.
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of `bayessurvreg' functions to make them more
  readable for the programmer.
}
\usage{
bayessurvreg.design(m, formula, random, data, transform, dtransform)

bayessurvreg.priorBeta(prior.beta, nX, indb, factors,
   n.factors, n.in.factors)

bayessurvreg.priorb(prior.b, nrandom, ncluster, toler.chol)
}
\arguments{
  \item{m, formula, random, data, transform, dtransform}{}
  \item{prior.beta, nX, indb, factors, n.factors, n.in.factors}{}
  \item{prior.b, nrandom, ncluster, toler.chol}{}
}
\value{
  Some lists.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\keyword{internal}
\keyword{connection}
