\name{plot.bayesDensity}
\alias{plot.bayesDensity}
\title{
  Plot an object of class bayesDensity
}  
\description{
  This function plots an object created by \code{\link{bayesDensity}}.
}
\usage{
\method{plot}{bayesDensity}(x, k.cond, dim.plot = TRUE, over = TRUE,
   alegend = TRUE, standard = TRUE, center = FALSE,
   type = "l", bty = "n",
   xlab = expression(epsilon), ylab = expression(f(epsilon)),
   lty, xlim, ylim, xleg, yleg, main, \dots)
}
\arguments{
  \item{x}{an object of class \code{bayesDensity}.}
  \item{k.cond}{a numerical vector giving the numbers of mixture components for which the conditional densities
     are to be plotted. 0 states for the unconditional (overall) density, averaged over the mixture with all possible
     numbers of components. If NULL, all conditional and the
     unconditional density found in \code{x} will be plotted.}
   \item{dim.plot}{an indicator whether the dimension of the plot used
     in \code{par(mfrow)} should be computed automatically. If
     \code{dim.plot = FALSE} and \code{over = FALSE} the user has to
     determine himself using \code{par(mfrow)} how to put the plots on
     the page.}
   \item{over}{an indicator whether all densities should be drawn into
     one plot using different types of lines. If \code{FALSE} a separate
     plot for each density is created.}
   \item{alegend}{an indicator whether an automatic legend should be
     added to the plot.}
   \item{standard}{logical, do we want to plot standardized density?}
   \item{center}{logical, do we want to plot centered density?, set both
     \code{standard} and \code{center} to \code{FALSE} if you wish to
     plot unstandardized density.}
   \item{xleg, yleg}{position of the legend if \code{over = TRUE}.}
   \item{type, bty, xlab, ylab, lty, xlim, ylim, main, \dots}{other
     arguments passed to the \code{plot.default} function.}
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\keyword{hplot}

