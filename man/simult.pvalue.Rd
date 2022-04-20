\name{simult.pvalue}
\alias{simult.pvalue}
\alias{print.simult.pvalue}
\title{
  Compute a simultaneous p-value from a sample for a vector valued parameter.
}
\description{
  The p-value is computed as 1 - the credible level of the credible
  region which just cover the point (0, 0, ..., 0)'.

  The function returns also the simultaneous credible region (rectangle)
  with a specified credible level.
}
\usage{
simult.pvalue(sample, precision=0.001, prob=0.95)
\method{print}{simult.pvalue}(x, \dots)
}
\arguments{
  \item{sample}{a data frame or matrix with sampled values (one column =
    one parameter)}
  \item{precision}{precision with which the p-value is to be computed}
  \item{prob}{probability for which the credible region is to be
    computed}
  \item{x}{an object of class simult.pvalue}
  \item{\dots}{who knows}
}
\value{
An object of class 'simult.pvalue'.
}  
\examples{
m <- 1000
sample <- data.frame(x1=rnorm(m), x2=rnorm(m), x3=rnorm(m))
simult.pvalue(sample)

sample <- data.frame(x1=rnorm(m), x2=rnorm(m), x3=rnorm(m, mean=2))
simult.pvalue(sample)

sample <- data.frame(x1=rnorm(m), x2=rnorm(m), x3=rnorm(m, mean=5))
simult.pvalue(sample, prob=0.99, precision=0.0001)
}
\references{
  Besag, J., Green, P., Higdon, D. and Mengersen, K. (1995).
  Bayesian computation and stochastic systems (with Discussion).
  \emph{Statistical Science,} \bold{10}, 3 - 66. page 30

  Held, L. (2004).
  Simultaneous posterior probability statements from Monte Carlo output.
  \emph{Journal of Computational and Graphical Statistics,} \bold{13}, 20 - 35. 
}  
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\keyword{htest}
