\name{rMVNorm}
\alias{rMVNorm}
\alias{C_rmvnormR2006}
\title{
  Sample from the multivariate normal distribution
}
\description{
  According to the parametrization used, sample from the multivariate
  normal distribution.

  The following parametrization can be specified
  \describe{
    \item{standard}{
      In this case we sample from either
      \eqn{\mathcal{N}(\mu, \Sigma)}{N(mu, Sigma)} or from
      \eqn{\mathcal{N}(\mu, Q^{-1}).}{N(mu, Q^(-1)).}      
    }
    \item{canonical}{
      In this case we sample from
      \eqn{\mathcal{N}(Q^{-1}b,\;Q^{-1}).}{N(Q^(-1)*b, Q^(-1)).}      
    }
  }

  Generation of random numbers is performed by Algorithms 2.3-2.5 in Rue and
  Held (2005, pp. 34-35). 
}
\usage{
  rMVNorm(n, mean=0, Sigma=1, Q, param=c("standard", "canonical"))
}
\arguments{
  \item{n}{number of observations to be sampled.}    
  \item{mean}{
    For \code{param="standard"}, it is 
    a vector \eqn{\mu}{mu} of means. If \code{length(mean)} is equal to 1, it
    is recycled and all components have the same mean.
 
    For \code{param="canonical"}, it is 
    a vector \eqn{b}{b} of canonical means. If \code{length(mean)} is equal to 1, it
    is recycled and all components have the same mean.    
  }
  \item{Sigma}{covariance matrix of the multivariate normal
    distribution. It is ignored if \code{Q} is given at the same time.
  }
  \item{Q}{precision matrix of the multivariate normal distribution.

    It does not have to be supplied provided \code{Sigma} is given and
    \code{param="standard"}.

    It must be supplied if \code{param="canonical"}.
  }
  \item{param}{a character which specifies the parametrization.
  }  
}
\value{
  Matrix with sampled points in rows.
}
\references{
  Rue, H. and Held, L. (2005).
  \emph{Gaussian Markov Random Fields: Theory and Applications}.
  Boca Raton: Chapman and Hall/CRC.  
}
\seealso{
  \code{\link{rnorm}}, \code{\link[mvtnorm]{Mvnorm}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
### Mean, covariance matrix, its inverse
### and the canonical mean
mu <- c(0, 2, 0.5)
L <- matrix(c(1, 1, 1,  0, 0.5, 0.5,  0, 0, 0.3), ncol=3)
Sigma <- L \%*\% t(L)
Q <- chol2inv(t(L))
b <- Q \%*\% mu

print(Sigma)
print(Q)
print(Sigma \%*\% Q)

### Sample using different parametrizations
set.seed(775988621)
n <- 10000

### Sample from N(mu, Sigma)
xx1 <- rMVNorm(n=n, mean=mu, Sigma=Sigma)
apply(xx1, 2, mean)
var(xx1)

### Sample from N(mu, Q^{-1})
xx2 <- rMVNorm(n=n, mean=mu, Q=Q)
apply(xx2, 2, mean)
var(xx2)

### Sample from N(Q^{-1}*b, Q^{-1})
xx3 <- rMVNorm(n=n, mean=b, Q=Q, param="canonical")
apply(xx3, 2, mean)
var(xx3)
}
\keyword{distribution}
\keyword{multivariate}

