\name{rWishart}
\alias{rWishart}
\title{
  Sample from the Wishart distribution
}
\description{
  Sample from the Wishart distribution
  \deqn{\mbox{Wishart}(\nu, S),}{Wishart(nu, S),}
  where \eqn{\nu}{nu} are degrees of freedom of the Wishart distribution
  and \eqn{S}{S} is its scale matrix. The same parametrization as in
  Gelman (2004) is assumed, that is, if
  \eqn{W\sim\mbox{Wishart}(\nu,S)}{W~Wishart(nu,S)} then
  \deqn{\mbox{E}(W) = \nu S}{E(W) = nu*S}.

  In the univariate case, \eqn{\mbox{Wishart}(\nu,S)}{Wishart(nu,S)} is the
  same as \eqn{\mbox{Gamma}(\nu/2, 1/(2S)).}{Gamma(nu/2, 1/(2*S)).}  
  
  Generation of random numbers is performed by the algorithm described
  in Ripley (1987, pp. 99).
}
\usage{
rWishart(n, df, S)
}
\arguments{
  \item{n}{number of observations to be sampled.}    
  \item{df}{degrees of freedom of the Wishart distribution.
  }
  \item{S}{scale matrix of the Wishart distribution.
  }
}
\value{
  Matrix with sampled points (lower triangles of \eqn{W}{W}) in rows.
}
\references{
  Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B. (2004).
  \emph{Bayesian Data Analysis, Second edition}. 
  Boca Raton: Chapman and Hall/CRC.

  Ripley, B. D. (1987).
  \emph{Stochastic Simulation}.
  New York: John Wiley and Sons.  
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{komarek@karlin.mff.cuni.cz}
}  
\examples{
### The same as rgamma(n, shape=df/2, rate=1/(2*S))
n <- 1000
df <- 1
S  <- 3
w <- rWishart(n=n, df=df, S=S)
mean(w)    ## should be close to df*S
var(w)     ## should be close to 2*df*S^2

### Multivariate Wishart
n <- 1000
df <- 2
S <- matrix(c(1,3,3,13), nrow=2)
w <- rWishart(n=n, df=df, S=S)
apply(w, 2, mean)                ## should be close to df*S
df*S

df <- 2.5
S <- matrix(c(1,2,3,2,20,26,3,26,70), nrow=3)
w <- rWishart(n=n, df=df, S=S)
apply(w, 2, mean)                ## should be close to df*S
df*S
}
\keyword{distribution}
\keyword{multivariate}

