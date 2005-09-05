\name{plot.bayesGspline}
\alias{plot.bayesGspline}
\title{
  Plot an object of class bayesGspline
}  
\description{
  This function plots an object created by \code{\link{bayesGspline}}.
}
\usage{
\method{plot}{bayesGspline}(x, add = FALSE, type = "l", lty=1, bty = "n",
     xlab, ylab, main, sub, \dots)
}
\arguments{
  \item{x}{an object of class \code{bayesGspline}.}
  \item{add}{if \code{TRUE} a new plot is produced, otherwise it is
    drawn to an existing plot.}
  \item{type, lty, bty, xlab, ylab, main, sub, \dots}{other
     arguments passed to the \code{plot.default} function.}
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\keyword{hplot}
