\name{files.Gspline}
\alias{clean.Gspline}
\alias{write.headers.Gspline}
\title{
  Write headers to or clean files with sampled G-spline
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of some \code{*.writeHeaders} functions to make them more
  readable for the programmer.  
}
\usage{
clean.Gspline(dir, label, care.of.y=TRUE)

write.headers.Gspline(dir, dim, nP, label, gparmi,
    store.a, store.y, store.r, care.of.y=TRUE)
}
\arguments{
\item{dir}{a~string giving the path to the directory where to search or
  store G-spline files}
\item{label}{a~string with extension of the G-spline files}
\item{dim}{dimension of the G-spline}
\item{nP}{number of (bivariate) observations}
\item{gparmi}{a~vector with integer arguments required by the
  constructor of C++ class \code{Gspline}}
\item{store.a}{logical, store transformed weights?}
\item{store.y}{logical, store augmented observations?}
\item{store.r}{logical, store allocations?}
\item{care.of.y}{logical, do we wish to take care of the \code{Y.sim} file?}
}
\value{
  No return value.
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\keyword{internal}
