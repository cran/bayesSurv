\name{credible.region}
\alias{credible.region}
\title{
  Compute a simultaneous confidence region from a sample for a vector valued parameter.
}
\description{
  See the two references for more detail.
}
\usage{
credible.region(sample, probs=c(0.90, 0.975))
}
\arguments{
\item{sample}{a data frame or matrix with sampled values (one column = one parameter)}
\item{probs}{probabilities for which the credible regions are to be computed}
}
\value{
A list (one component for each confidence region) of length equal to
\code{length(probs)}. Each component of the list is a matrix with two
rows (lower and upper limit) and as many columns as the number of
parameters giving the confidence region.
}  
\examples{
  m <- 10000
  sample <- data.frame(x1=rnorm(m), x2=rnorm(m), x3=rnorm(m))
  probs <- c(0.70, 0.90, 0.95)
  CR <- credible.region(sample, probs=probs)

  for (kk in 1:length(CR)){
    suma <- sum(sample$x1 >= CR[[kk]]["Lower", "x1"] & sample$x1 <= CR[[kk]]["Upper", "x1"] &
            sample$x2 >= CR[[kk]]["Lower", "x2"] & sample$x2 <= CR[[kk]]["Upper", "x2"] &
            sample$x3 >= CR[[kk]]["Lower", "x3"] & sample$x3 <= CR[[kk]]["Upper", "x3"])
    show <- c(suma/m, probs[kk])
    names(show) <- c("Empirical", "Desired")
    print(show)
  }
}
\references{
  Besag, J., Green, P., Higdon, D. and Mengersen, K. (1995).
  Bayesian computation and stochastic systems (with Discussion).
  \emph{Statistical Science,} \bold{10}, 3 - 66. page 30

  Held, L. (2004).
  Simultaneous inference in risk assessment; a Bayesian perspective
  \emph{In: COMPSTAT 2004, Proceedings in Computational Statistics (J. Antoch, Ed.),}
  213 - 222, page 214
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\keyword{htest}
