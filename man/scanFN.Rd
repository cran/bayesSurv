\name{scanFN}
\alias{scanFN}
\title{Read Data Values}
\description{
  Read numeric data into a data frame from a file. Header is assumed to
  be present in the file.
}
\usage{
scanFN(file, quiet=FALSE)
}
\arguments{
  \item{file}{the name of a file to read data values from.  If the
    specified file is \code{""}, then input is taken from the keyboard
    (or \code{stdin} if input is redirected).
    (In this case input can be terminated by a blank line or an EOF
    signal, \code{Ctrl-D} on Unix and \code{Ctrl-Z} on Windows.)

    Otherwise, the file name is interpreted \emph{relative} to the
    current working directory (given by \code{\link{getwd}()}),
    unless it specifies an \emph{absolute} path.
    Tilde-expansion is performed where supported.

    Alternatively, \code{file} can be a \code{\link{connection}},
    which will be opened if necessary, and if so closed at the end of
    the function call.  Whatever mode the connection is opened in,
    any of LF, CRLF or CR will be accepted as the EOL marker for a line
    and so will match \code{sep = "\n"}.

    \code{file} can also be a complete URL.

    To read a data file not in the current encoding (for example a
    Latin-1 file in a UTF-8 locale or conversely) use a
    \code{\link{file}} connection setting the \code{encoding} argument.
  }

  \item{quiet}{logical: if \code{FALSE} (default), scan() will print a
    line, saying how many items have been read.}

}
\details{
  See \code{\link{scan}}.
}
\value{
  \code{data.frame} with read data values.
}
\references{
  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
  \emph{The New S Language}.
  Wadsworth & Brooks/Cole.
}
\seealso{
  \code{\link{scan}}
}
\author{
  Arnošt Komárek \email{arnost.komarek@mff.cuni.cz}
}
\examples{
cat("x y z", "1 2 3", "1 4 6", "10 20 30", file="ex.data", sep="\n")
pp <- scanFN("ex.data", quiet=FALSE)
pp <- scanFN("ex.data", quiet= TRUE)
print(pp)
unlink("ex.data") # tidy up
}
\keyword{file}
\keyword{connection}
