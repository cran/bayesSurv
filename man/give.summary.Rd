\name{give.summary}
\alias{give.summary}
\title{
  Brief summary for the chain(s) obtained using the MCMC.
}
\description{
  This function computes a sample mean, quantiles and a Bayesian
  \eqn{p}{p}-value which is defined as
  \deqn{p = 2\times\min(n_{-},\,n_{+}),}{2 * min(n[-], n[+]),}
  where
  \eqn{n_{-}}{n[-]} is the number of the sampled values which are
  negative and \eqn{n_{+}}{n[+]} is the number of sampled values which
  are positive.
}
\usage{
give.summary(sample, probs=c(0.5, 0.025, 0.975))
}
\arguments{
  \item{sample}{a data frame or a vector with sampled values}
  \item{probs}{probabilities of the quantiles that are to be computed}
}  
\value{
  A matrix or a vector with the sample mean, quantiles and a Bayesian \eqn{p}{p}-value.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\examples{
## Example with a sample stored in a vector:
sample <- rnorm(1000)
give.summary(sample)

## Example with a sample stored in a data.frame:
sample <- data.frame(x=rnorm(1000), y=rgamma(1000, shape=1, rate=1))
give.summary(sample)
}
\keyword{univar}
