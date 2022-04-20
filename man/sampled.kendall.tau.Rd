\name{sampled.kendall.tau}
\alias{sampled.kendall.tau}
\alias{C_sampledKendallTau}
\title{
  Estimate of the Kendall's tau from the bivariate model
}
\description{
  This function computes an estimate of the residual (after adjustment
  for covariates) Kendall's tau for the bivariate survival model fitted
  using the functions \code{\link{bayesHistogram}} or
  \code{\link{bayesBisurvreg}}.

  For both these function their argument \code{prior$specification} must
  be equal to 2!

  When \eqn{G} is a bivariate distribution function, the population
  version of the Kendall's tau is defined as
  \deqn{\tau = 4\int G dG - 1}{tau = 4*int G dG - 1}.

  For the model estimated using one of the above mentioned functions the
  value of Kendall's tau at each iteration of MCMC is equal to
  \deqn{\tau =
    4\sum_{i=-K_1}^{K_1}\sum_{j=-K_2}^{K_2}\sum_{k=-K_1}^{K_1}\sum_{l=-K_2}^{K_2}w_{i,j} w_{k,l}
    \Phi\left(\frac{\mu_{1,i} - \mu_{1,k}}{\sqrt{2}\sigma_1}\right)
    \Phi\left(\frac{\mu_{2,j} - \mu_{2,l}}{\sqrt{2}\sigma_2}\right)
    - 1,}{%
    4*sum[i=-K[1]][K[1]]sum[j=-K[2]][K[2]]sum[k=-K[1]][K[1]]sum[l=-K[2]][K[2]] w[i,j]*w[k,l]*
    Phi((mu[1,i] - mu[1,k])/(sqrt(2)*sigma[1])) * Phi((mu[2,j] - mu[2,l])/(sqrt(2)*sigma[2]))
    - 1,}
  where \eqn{\mu_{1,-K_1},\dots,\mu_{1,K_1}}{mu[1,-K[1]],...,mu[1,K[1]]}
  are knots in the first margin, 
  \eqn{\mu_{2,-K_2},\dots,\mu_{2,K_2}}{mu[2,-K[2]],...,mu[2,K[2]]}
  are knots in the second margin,
  \eqn{\sigma_1}{sigma[1]} is the basis standard deviation in the first margin,
  \eqn{\sigma_2}{sigma[2]} is the basis standard deviation in the second margin,
  and \eqn{w_{i,j},\;i=-K_1,\dots,K_1,\;j=-K_2,\dots,K_2}{
       w[i,j], i=-K[1],...,K[1], j=-K[2],...,K[2]} are the G-spline weights.  
}
\usage{
sampled.kendall.tau(dir = getwd(), extens = "", K,
  skip = 0, by = 1, last.iter, nwrite)
}
\arguments{%%%AAA
\item{dir}{directory where to search for files (`mixmoment.sim',
  `mweight.sim', `mmean.sim', `gspline.sim') with the MCMC
  sample.}
\item{extens}{an extension used to distinguish different sampled
  G-splines if more G-splines were used in one simulation (with
  doubly-censored data) According to which
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
  }
}
\item{K}{a~vector of length 2 specifying the number of knots at each side of the middle knot for each dimension of the
      G-spline.}
\item{skip}{number of rows that should be skipped at the beginning of
  each *.sim file with the stored sample.}
\item{by}{additional thinning of the sample.}
\item{last.iter}{index of the last row from *.sim files that should be
  used. If not specified than it is set to the maximum available
  determined according to the file \code{mixmoment.sim}.}
\item{nwrite}{frequency with which is the user informed about the
  progress of computation (every \code{nwrite}th iteration count of
  iterations change).}
}
\value{
  A vector with sampled values of the Kendall's tau.
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
