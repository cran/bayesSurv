\name{bayesGspline}
\alias{bayesGspline}
\title{
  Summary for the density estimate based on the model with Bayesian
  G-splines.
}
\description{
   Compute the estimate of the density function based on the values
   sampled using the MCMC (MCMC average evaluated in a
   grid of values) in a model where density is specified as
   a Bayesian G-spline.

   This function serves to summarize the MCMC chains related to the distributional parts 
   of the considered models obtained using the functions:
   \code{\link{bayesHistogram}},
   \code{\link{bayesBisurvreg}}, \code{\link{bayessurvreg2}}, \code{\link{bayessurvreg3}}.

   If asked, this function returns also the values of the G-spline
   evaluated in a grid at each iteration of MCMC.
}
\usage{
bayesGspline(dir = getwd(), extens="", extens.adjust="_b",
   grid1, grid2, skip = 0, by = 1, last.iter, nwrite,
   only.aver = TRUE, standard = FALSE, version = 0)
}
\arguments{%%%AAA
\item{dir}{directory where to search for files (`mixmoment.sim',
  `mweight.sim', `mmean.sim', `gspline.sim') with the MCMC
  sample.}
\item{extens}{an extension used to distinguish different sampled
  G-splines if more G-splines were used in one simulation (e.g. with
  doubly-censored data or in the model where both the error term and the
  random intercept were defined as the G-splines). According to which
  \code{bayes*survreg*} function was used, specify the argument
  \code{extens} in the following way.
  
  \describe{
    \item{bayesHistogram:}{
      always \code{extens = ""}
    }
    \item{bayesBisurvreg:}{\eqn{\quad}{ }

      \itemize{
         \item to compute the bivariate distribution of the \emph{error} term for
           the \emph{onset} time: \code{extens = ""};

	 \item to compute the bivariate distribution of the \emph{error} term for
           the \emph{event} time if there was doubly-censoring:
           \code{extens = "_2"};
      }	 
    }
    \item{bayessurvreg2:}{\eqn{\quad}{ }

      \itemize{
         \item to compute the distribution of the \emph{error} term for
           the \emph{onset} time: \code{extens = ""};

         \item to compute the distribution of the \emph{error} term for
           the \emph{event} time if there was doubly-censoring:
           \code{extens = "_2"};
      }	   
    }
    \item{bayessurvreg3:}{\eqn{\quad}{ }

      \itemize{
         \item to compute the distribution of the \emph{error} term for
           the \emph{onset} time: \code{extens = ""};

         \item to compute the distribution of the \emph{error} term for
           the \emph{event} time if there was doubly-censoring:
           \code{extens = "_2"};

         \item to compute the distribution of the \emph{random intercept} for the
           \emph{onset} time: \code{extens = "_b"};

         \item to compute the distribution of the \emph{random intercept} term for
           the \emph{event} time if there was doubly-censoring:
           \code{extens = "_b2"};      
      }
    }  
  }
}
\item{extens.adjust}{this argument is applicable for the situation when
  the MCMC chains were created using the function
  \code{\link{bayessurvreg3}}, i.e. when both the distribution of the
  error term and the random intercept was specified as the G-spline.
  
  In that case the location of the error term and the random intercept
  are separately not identifiable. Only the location of the sum
  \eqn{\varepsilon + b}{epsilon + b} can be estimated. For this reason,
  the function \code{bayesGspline} always centers the distribuytion of
  the random intercept to have a zero mean and adds its original mean to
  the mean of the distribution of the error term.

  Argument \code{extens.adjust} is used to match correctly the files
  containing the G-spline of the random intercept corresponding to the
  particular error term. 

  The following values of \code{extens.adjust} should be used in the
  following situations:

  \itemize{
    \item if there are no doubly-censored data or if we are computing
    the distribution of the error term/random intercept from the model
    for the \emph{onset} time then
    \tabular{c}{
      \code{extens.adjust = "_b"}
    }

    \item if there are doubly-censored data and we are computing the
    distribution of the error term/random intercept from the model for
    the \emph{event} time then
    \tabular{c}{
      \code{extens.adjust = "_b2"}
    }
  }  
}  %% end of item{extens\_adjust} 
\item{grid1}{grid of values from the first dimension at which the sampled 
  densities are to be evaluated.}
\item{grid2}{grid of values from the second dimension (if the G-spline
  was bivariate) at which the sampled densities are to be
  evaluated. This item is \code{missing} if the G-spline is univariate.}
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
  returned otherwise also values of the G-spline at each iteration are
  returned (which might ask for quite lots of memory).}
\item{standard}{\code{TRUE/FALSE}, if \code{TRUE}, each G-spline is
  standardized to have zero mean and unit variance. Only applicable if
  \code{version} = 30 or 31, otherwise \code{standard} is always set to \code{FALSE}.}
\item{version}{
  this argument indicates by which \code{bayes*survreg*} function the
  chains used by \code{bayesGspline} were created. Use the following:

  \describe{
    \item{bayesHistogram:}{\code{version = 0};}
    \item{bayesBisurvreg:}{\code{version = 0};}
    \item{bayessurvreg2:}{\code{version = 0};}
    \item{bayessurvreg3:}{
      \code{version = 30 } or \code{ 31}.

       Use
       \code{version = 30} if you want to compute the density of the
       \emph{error} term;

       \code{version = 31} if you want to compute the density of the
       \emph{random intercept}.
    }
  }
}
}
\value{
  An object of class \code{bayesGspline} is returned. This object is a
  list with components
  \code{grid}, \code{average} for the univariate G-spline and
  components \code{grid1}, \code{grid2}, \code{average} for the bivariate G-spline.
  \item{grid}{this is a grid of values (vector) at which the McMC average of the
    G-spline was computed.}
  \item{average}{these are McMC averages of the G-spline (vector) evaluated in
    \code{grid}.}
  \item{grid1}{this is a grid of values (vector) for the first dimension
    at which the McMC average of the G-spline was computed.}
  \item{grid2}{this is a grid of values (vector) for the second dimension
    at which the McMC average of the G-spline was computed.}
  \item{average}{this is a matrix \code{length(grid1)} times
    \code{length(grid2)} with McMC averages of the G-spline evaluated in
    \tabular{lccc}{
      x1 = \tab ( grid1 \tab \dots \tab grid1 ) \cr
    }
    and
    \tabular{lrcl}{
           \tab (\tab grid2  \tab ) \cr
      x2 = \tab (\tab \dots    \tab )  \cr
	   \tab (\tab grid2  \tab ) \cr   
    }
  }

  There exists a method to plot objects of the class \code{bayesGspline}.
}
\section{Attributes}{
  Additionally, the object of class \code{bayesGspline} has the following
  attributes:
  \describe{
    \item{\code{sample.size}}{a length of the McMC sample used to compute the McMC average.
    }
    \item{\code{sample}}{G-spline evaluated in a grid of values. This attribute
      is present only if \code{only.aver = FALSE}.

      For a univariate G-spline this is a matrix with \code{sample.size} columns and
      length(grid1) rows.

      For a bivariate G-spline this is a matrix
      with \code{sample.size} columns and
      length(grid1)*length(grid2) rows.
    }
  }    
}
\references{
  \eqn{\mbox{Kom\'{a}rek, A.}}{Komarek, A.} (2006).
  \emph{Accelerated Failure Time Models for Multivariate
    Interval-Censored Data with Flexible Distributional Assumptions}.
  PhD. Thesis, Katholieke Universiteit Leuven, Faculteit Wetenschappen.

  \eqn{\mbox{Kom\'{a}rek}}{Komarek}, A. and Lesaffre, E. (2006).
  Bayesian accelerated failure time model with multivariate doubly-interval-censored data
  and flexible distributional assumptions.
  \emph{Submitted}.  \cr
  See \code{Komarek_Lesaffre_2006.pdf}.

  \eqn{\mbox{Kom\'{a}rek}}{Komarek}, A. and Lesaffre, E. (2006b).
  Bayesian semi-parametric accelerated failurew time model for paired
  doubly interval-censored data.
  \emph{Statistical Modelling}, \bold{6}, 3--22.
  
  \eqn{\mbox{Kom\'{a}rek}}{Komarek}, A. Lesaffre, E., and Legrand, C. (2006).
  EORTC data article.
  \emph{To be written.}
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\examples{
## See the description of R commands for
## the models described in
## Komarek (2006),
## Komarek and Lesaffre (2006),
## Komarek and Lesaffre (2006b),
## Komarek, Lesaffre, and Legrand (2006).
## 
## R commands available (or soon available)
## in the documentation
## directory of this package
## as tandmobCS.pdf, tandmobCS.R
##    tandmobPA.pdf, tandmobPA.R.
##
}
\keyword{smooth}
