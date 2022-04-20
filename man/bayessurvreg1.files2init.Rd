\name{bayessurvreg1.files2init}
\alias{bayessurvreg1.files2init}
\title{
   Read the initial values for the Bayesian survival regression model to the list.
}
\description{
   This function creates the list of initial values as required by the \code{init}
   argument of the function \code{\link{bayessurvreg1}}. The initials are taken from
   the files that are of the form of the files where the simulated values from
   the McMC run performed by the function \code{\link{bayessurvreg1}} are stored.
   The files are assumed to have the following names: "iteration.sim",
   "mixmoment.sim", "mweight.sim", "mmean.sim", "mvariance.sim",   
   "beta.sim", "b.sim", "Y.sim", "r.sim", "D.sim", "otherp.sim", "u.sim". Some of these files
   may be missing. In that case, the corresponding initial is filled by \code{NULL}.
}
\usage{
bayessurvreg1.files2init(dir = getwd(), row, kmax)
}
\arguments{
  \item{dir}{string giving the directory where it will be searched for the files with 
    initial values.}
  \item{row}{the row (possible header does not count) from the files with the
    values that will be considered to give the initial values. By default,
    it is the last row from the files.}
  \item{kmax}{maximal number of mixture components. This must be given
    only if \code{header == FALSE}.}
}
\value{
  A list with components called "iter", "mixture", "beta", "b", "D", "y",
  "r", "otherp", "u"
  in the form as required by the argument \code{init} of the function \code{\link{bayessurvreg1}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\keyword{file}
