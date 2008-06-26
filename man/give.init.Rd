\name{give.init}
\alias{give.init.Gspline}
\alias{give.init.y}
\alias{give.init.y2}
\alias{give.init.r}
\title{
  Check and possibly fill in initial values for the G-spline, augmented
  observations and allocations for Bayesian models with G-splines
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of \code{\link{bayesBisurvreg.priorInit}} and
  related functions to make them more readable for the programmer.
}
\usage{
give.init.Gspline(prior, init, mcmc.par, dim)

give.init.y(init.y, dim, y.left, y.right, status)

give.init.y2(init.y, init2.y, dim, design, design2, doubly)

give.init.r(init.r, init.y, dim, KK,
   gamma, sigma, c4delta, intcpt, scale)
}
\arguments{
  \item{prior}{a~list as required by \code{prior} argument of the function
    \code{\link{bayesHistogram}} or \code{prior} and \code{prior2}
    arguments of the function \code{\link{bayesBisurvreg}}
  }  %% end of item{prior}
  \item{init}{a~list as required by \code{init} argument of the function
    \code{\link{bayesHistogram}} or by \code{init} and \code{init2}
    arguments of the function \code{\link{bayesBisurvreg}}
  }  %% end of item{init}
  \item{mcmc.par}{a~list as required by \code{mcmc.par} argument of
    function \code{\link{bayesHistogram}} or by \code{mcmc.par} and
    \code{mcmc.par2} arguments of the function
    \code{\link{bayesBisurvreg}}
  }  %% end of item{mcmc.par}
  \item{dim}{dimension of the G-spline/response, 1 or 2.
  }  %% end of item{dim}
  \item{init.y}{initial (augmented) observations possibly given by the user.                        
      They are partially checked for consistency and these supplied by the user                   
      used in the resulting object. This should be either vector of length \eqn{n}{n}             
      where \eqn{n}{n} is a~sample size if the dimension is one or a~matrix with                  
      2 columns and \eqn{n}{n} rows if the dimension is two.
  }  %% end of item{init.y}
  \item{init2.y}{initial (augmented) times-to-event (if doubly censoring) possibly given by the user.  
      They are partially checked for consistency and these supplied by the user                   
      used in the resulting object. This should be either vector of length \eqn{n}{n}             
      where \eqn{n}{n} is a~sample size if the dimension is one or a~matrix with                  
      2 columns and \eqn{n}{n} rows if the dimension is two.
  }  %% end of item{init2.y}
  \item{design}{an~object as returned by the function
    \code{\link{bayessurvreg.design}} related to either the onset time
    if doubly censored observations or to the event time. Remark:
    \code{design$Y} contains original times and NOT their logarithmic
    transformations.
  }  %% end of item{design}
  \item{design2}{an~object as returned by the function
    \code{\link{bayessurvreg.design}} related to time-to-event
    if doubly censored observations. Remark:
    \code{design2$Y} contains original times and NOT their logarithmic
    transformations.
  }  %% end of item{design2}    
  \item{doubly}{logical indicating whether the response is doubly
    censored or not
  }  %% end of item{doubly}
  \item{y.left}{observed, left or right censored log(event time) or the lower limit                 
    of the interval censored observation. Sorted in a~transposed order compared                 
    to \code{init.y}.
  }  %% end of item{y.left}           
  \item{y.right}{upper limit of the interval censored observation, whatever if the observation      
    is not interval-censored sorted in a~transposed order compared to \code{init.y}.
  }  %% end of item{y.right}
  \item{status}{status indicator vector/matrix. 1 for observed times, 0 for right censored times,   
    2 for left censored times, 3 for interval censored times.
  }  %% end of item{status}
  \item{init.r}{initial allocations possibly given by the user. This should be a~vector of length
     \eqn{n}{n} where n is a~sample size if \code{dim} is equal to 1 and
     a~matrix with \eqn{n}{n} rows and 2 columns if \code{dim} is equal to
     2. Values should be on the scale
     \eqn{-K_j,\dots,K_j,}{-K[j],...,K[j],}
     \eqn{j=1,\dots,dim}{j=1,...,dim}
  }  %% end of item{init.r}
  \item{init.y}{correctly computed initial observations the G-spline is estimated from. In the case of regression
     these should be replaced by residuals. This must be either a~vector or matrix (in the format as returned
     by the function \code{\link{give.init.y}}).
  }  %% end of item{init.y}
  \item{KK}{vector of length \code{dim} with \eqn{K}{K} coefficients
    defining the G-spline.
  }  %% end of item{KK}
  \item{gamma}{vector of length \code{dim} with initial
    \eqn{\gamma}{gamma} parameters of the G-spline.
  }  %% end of item{gamma}
  \item{sigma}{vector of length \code{dim} with initial
    \eqn{\sigma}{sigma} parameters of the G-spline.
  }  %% end of item{sigma}
  \item{c4delta}{vector of length \code{dim} with constants to compute the distance between two knots
    defining the G-spline.
  }  %% end of item{c4delta}
  \item{intcpt}{vector of length \code{dim} with initial values of the
    intercept term of the G-spline.
  }  %% end of item{intcpt}
  \item{scale}{vector of length \code{dim} with initial values of the
    scale parameters of the G-spline.
  }  %% end of item{scale}
}
\value{
  Some lists.
}
\section{Value for give.init.Gspline}{
  A~list with the following components:
  \describe{
  \item{Gparmi}{integer parameters for the G-spline constructor in the
    C++ code}
  \item{Gparmd}{double parameters for the G-spline constructor in the
    C++ code}
  \item{specification}{specification of the G-spline model (1 or 2), see
    \code{\link{bayesHistogram}} for more detail}
  }

  and the following attributes:
  \tabular{l}{
  \bold{init}      \cr
  \bold{prior}     \cr
  \bold{mcmc.par}  \cr
  }
}
\section{Value for give.init.y}{
   A~vector or matrix with the same structure as \code{init.y}, i.e. with 2~columns and
   \eqn{n}{n} rows in the case of the bivariate data.
}
\section{Value for give.init.y2}{
  A~list with the following components:
  \describe{
    \item{init.y}{a~vector of length \eqn{n}{n} or a~matrix with 2
      columns and \eqn{n}{n} rows with initial log(onset times)}
    \item{init2.y}{a~vector of length \eqn{n}{n} or a~matrix with 2
      columns and \eqn{n}{n} rows with initial log(times-to-event). If
      the data are not doubly cdensored, this object is equal to 0.}
  \item{y1.left}{lower limit of the log-response (or exact/right/left
    censored observation) as required by the C++ function
    \code{bayesBisurvreg}, related to the onset time in the case of
    doubly censoring and to the event time otherwise}
  \item{y1.right}{upper limit of the log-response as required by the C++ function
    \code{bayesBisurvreg}, related to the onset time in the case of
    doubly censoring and to the event time otherwise}
  \item{status1}{status vector as required by the C++ function
    \code{bayesBisurvreg} related to the onset time in the case of
    doubly censoring and to the event time otherwise}
  \item{t2.left}{lower limit of the response as required by the C++
    function \code{bayesBisurvreg}, related to time-to-event in the
    case of doubly censoring, equal to 0 if there is no doubly-censoring}
  \item{t2.right}{upper limit of the response as required by the C++
    function \code{bayesBisurvreg}, related to time-to-event in the
    case of doubly censoring, equal to 0 if there is no doubly-censoring}
  \item{status2}{status vector related to time-to-event in the case of
    doubly censoring, equal to 0 otherwise.}
  }    
}  
\section{Value for give.init.r}{
  To be added somewhen...
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}
\keyword{internal}
