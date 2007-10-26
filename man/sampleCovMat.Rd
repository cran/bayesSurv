\name{sampleCovMat}
\alias{sampleCovMat}
\title{
   Compute a sample covariance matrix.
}
\description{
   This function computes a sample covariance matrix. 
}
\usage{
sampleCovMat(sample)
}
\arguments{
  \item{sample}{a \code{matrix} or \code{data.frame} with sampled values
    in rows. I.e. number of rows of \code{sample} determines a sample
    size, number of columns of \code{sample} determines a dimension of
    the distribution from which it was sampled.
  }     
}
\value{This function returns a matrix.
}
\details{
  When \eqn{y_1, \dots, y_n}{y[1], ..., y[n]} is a sequence of
  \eqn{p}{p}-dimensional vectors \eqn{y_i}{y[i]} the sample covariance
  matrix \eqn{S}{S} is equal to
  \deqn{S = \frac{1}{n-1} \sum_{i=1}^n (y_i - m)(y_i - m)^T}{%
        S = (1/(n-1)) sum[i=1][n] (y[i] - m)(y[i] - m)'}
      where \deqn{m = \frac{1}{n}\sum_{i=1}^n y_i.}{%
	          m = (1/n) sum[i=1][n] y[i].}
  When \eqn{n=1}{n=1} the function returns just sum of squares.      
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{komarek@karlin.mff.cuni.cz}
}
\examples{
  ## Sample some values
  z1 <- rnorm(100, 0, 1)           ## first components of y
  z2 <- rnorm(100, 5, 2)           ## second components of y
  z3 <- rnorm(100, 10, 0.5)        ## third components of y

  ## Put them into a data.frame
  sample <- data.frame(z1, z2, z3)

  ## Compute a sample covariance matrix
  sampleCovMat(sample)
}
\keyword{multivariate}
