\name{traceplot2}
\alias{traceplot2}
\title{
  Trace plot of MCMC output.
}
\description{  
  Displays a plot of iterations vs. sampled values for each variable in
  the chain, with a separate plot per variable.

  This is slightly modified version of \code{\link[coda]{traceplot}} function
  of a \code{coda} package to conform to my personal preferences.
}
\usage{
traceplot2(x, chains, bty = "n", main, xlab, ...)
}
\arguments{
  \item{x}{an \code{\link[coda]{mcmc}} or \code{\link[coda]{mcmc.list}} object.}
  \item{chains}{indeces of chains from the object that are to be plotted.}
  \item{bty, main, xlab, \dots}{further arguments passed to the
    \code{\link{plot.default}} function.}
}
\value{
  No return value.
}  
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\keyword{hplot}
