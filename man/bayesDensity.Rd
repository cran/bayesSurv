\name{bayesDensity}
\alias{bayesDensity}
\title{
   Summary for the density estimate based on the Bayesian model.
}
\description{
   Compute the conditional (given the number of mixture components) and
   unconditional estimate of the density function based on the values
   sampled using the reversible jumps McMC (McMC average evaluated in a
   grid of values). Give also the values of each sampled density
   evaluated at that grid (returned as the attribute of the resulting
   object). Methods for printing and plotting are also provided.
}
\usage{
bayesDensity(dir = getwd(), stgrid, grid,
    n.grid = 100, skip = 0, standard = TRUE, unstandard = TRUE)
}
\arguments{
   \item{dir}{directory where to search for files (`mixmoment.sim',
     `mweight.sim', mmean.sim', mvariance.sim') with the McMC sample.}
   \item{stgrid}{grid of values at which the sampled standardized
     densities are to be evaluated. If \code{missing}, the grid is
     automatically computed.}
   \item{grid}{grid of values at which the sampled densities are to be
     evaluated. If \code{missing},
     the grid is guessed from the first 20 sampled mixtures as the sequence starting with the minimal
     sampled mixture mean minus 3 standard deviations of the appropriate mixture
     component, ending with the maximal sampled mixture mean plus 3
     standard deviations of the appropriate mixture
     component, of the length given by \code{n.grid}.}
   \item{n.grid}{the length of the grid if \code{grid = NULL}.}
   \item{skip}{number of rows that should be skipped at the beginning of
   each file with the stored sample.}
   \item{standard}{if \code{TRUE} then also standardized (zero mean,
     unit variance) sampled densities are evaluated.}
   \item{unstandard}{of \code{TRUE} then also original (unstandardized)
     sampled densities are evaluated.}
}
\value{
  An object of class \code{bayesDensity} is returned. This object is a
  list and has potentially two components: \code{standard} and
  \code{unstandard}. Each of these two components is a \code{data.frame}
  with as many rows as number of grid points at which the density was
  evaluated and with columns called `grid', `unconditional' and `k = 1',
  \dots, `k = k.max' giving a predictive errr density, either averaged
  over all sampled \eqn{k}{k}s (unconditional) or averaged over a
  psecific number of mixture components.

  Additionally, the object of class \code{bayesDensity} has three
  attributes:
    \item{sample.size}{a vector of length \code{1 + kmax} giving the
      frequency of each \eqn{k}{k} in the sample.}
    \item{moments}{a \code{data.frame} with columns called `intercept'
      and `scale' giving the mean and variance of the sampled mixture at
      each iteration of the McMC.}
    \item{k}{a \code{data.frame} with one column called `k' giving
      number of mixture components at each iteration.}

  There exists methods to print and plot objects of the class \code{bayesDensity}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.ac.be}
}
\keyword{smooth}






