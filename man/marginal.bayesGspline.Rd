\name{marginal.bayesGspline}
\alias{marginal.bayesGspline}
\alias{C_marginal_bayesGspline}
\title{
  Summary for the marginal density estimates based on the bivariate model with Bayesian
  G-splines.
}
\description{
   Compute the estimate of the marginal density function based on the values
   sampled using the MCMC (MCMC average evaluated in a
   grid of values) in a model where density is specified as
   a bivariate Bayesian G-spline.

   This function serves to summarize the MCMC chains related to the distributional parts 
   of the considered models obtained using the functions:
   \code{\link{bayesHistogram}} and \code{\link{bayesBisurvreg}}.

   If asked, this function returns also the values of the marginal G-spline
   evaluated in a grid at each iteration of MCMC.
}
\usage{
marginal.bayesGspline(dir = getwd(), extens = "", K, grid1, grid2,
   skip = 0, by = 1, last.iter, nwrite, only.aver = TRUE)
}
\arguments{
\item{dir}{directory where to search for files (`mixmoment.sim',
  `mweight.sim', `mmean.sim', `gspline.sim') with the MCMC
  sample.}
\item{extens}{an extension used to distinguish different sampled
  G-splines if more G-splines were used in one simulation (e.g. with
  doubly-censored data). According to which
  \code{bayes*survreg*} function was used, specify the argument
  \code{extens} in the following way.
  
  \describe{
    \item{bayesHistogram:}{
      always \code{extens = ""}
    }
    \item{bayesBisurvreg:}{\eqn{\quad}{ }

      \itemize{
         \item to compute the marginals of the bivariate distribution of the \emph{error} term for
           the \emph{onset} time: \code{extens = ""};

	 \item to compute the marginals of the bivariate distribution of the \emph{error} term for
           the \emph{event} time if there was doubly-censoring:
           \code{extens = "_2"};
      }	 
    }
  }
}
\item{K}{a~vector of length 2 specifying the number of knots at each
  side of the middle knot for each dimension of the G-spline.}
\item{grid1}{grid of values from the first dimension at which the sampled 
  marginal densities are to be evaluated.}
\item{grid2}{grid of values from the second dimension
  at which the sampled marginal densities are to be
  evaluated.}
\item{skip}{number of rows that should be skipped at the beginning of
  each *.sim file with the stored sample.}
\item{by}{additional thinning of the sample.}
\item{last.iter}{index of the last row from *.sim files that should be
  used. If not specified than it is set to the maximum available
  determined according to the file \code{mixmoment.sim}.}
\item{nwrite}{frequency with which is the user informed about the
  progress of computation (every \code{nwrite}th iteration count of
  iterations change).}
\item{only.aver}{\code{TRUE/FALSE}, if \code{TRUE} only MCMC average is
  returned otherwise also values of the marginal G-spline at each iteration are
  returned (which might ask for quite lots of memory).}
}
\value{
  An object of class \code{marginal.bayesGspline} is returned. This object is a
  list with components \code{margin1} and \code{margin2} (for two margins).

  Both \code{margin1} and \code{margin2} components are data.frames with
  columns named \code{grid} and \code{average} where
  \item{grid}{is a grid of values (vector) at which the McMC average of the
    marginal G-spline was computed.}
  \item{average}{are McMC averages of the marginal G-spline (vector) evaluated in
    \code{grid}.}
  
  There exists a method to plot objects of the class \code{marginal.bayesGspline}.
}
\section{Attributes}{
  Additionally, the object of class \code{marginal.bayesGspline} has the following
  attributes:
  \describe{
    \item{\code{sample.size}}{a length of the McMC sample used to compute the McMC average.
    }
    \item{\code{sample1}}{marginal (margin = 1) G-spline evaluated in a grid of values. This attribute
      is present only if \code{only.aver = FALSE}.

      This a matrix with \code{sample.size} columns and length(grid1) rows.
    }
    \item{\code{sample2}}{marginal (margin = 2) G-spline evaluated in a grid of values. This attribute
      is present only if \code{only.aver = FALSE}.

      This a matrix with \code{sample.size} columns and length(grid2) rows.
    }    
  }    
}
\references{
  Komárek, A. (2006).
  \emph{Accelerated Failure Time Models for Multivariate
    Interval-Censored Data with Flexible Distributional Assumptions}.
    PhD. Thesis, Katholieke Universiteit Leuven, Faculteit Wetenschappen.

  Komárek, A. and Lesaffre, E. (2006).
  Bayesian semi-parametric accelerated failurew time model for paired
  doubly interval-censored data.
  \emph{Statistical Modelling}, \bold{6}, 3 - 22.
}  
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\examples{
## See the description of R commands for
## the models described in
## Komarek (2006),
## Komarek and Lesaffre (2006),
## 
## R commands available
## in the documentation
## directory of this package
## - see ex-tandmobPA.R and
##   https://www2.karlin.mff.cuni.cz/~komarek/software/bayesSurv/ex-tandmobPA.pdf
##
}
\keyword{smooth}

