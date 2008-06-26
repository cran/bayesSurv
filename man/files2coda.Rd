\name{files2coda}
\alias{files2coda}
\title{
   Read the sampled values from the Bayesian survival regression model
   to a coda mcmc object.
}
\description{
   This function creates a \code{coda} \code{mcmc} object from values found
   in files where sampled values from \code{bayessurvreg1} function are stored
   or from data.frames.
}
\usage{
files2coda(files, data.frames, variant = 1, dir = getwd(),
   start = 1, end, thin = 1, header = TRUE, chain)
}
\arguments{
  \item{files}{a vector of strings giving the names of files that are to
  be converted to \code{coda} objects. If missing and \code{data.frames}
  is also missing, all appropriate files
  found in a directory \code{dir} are converted to coda objects. File
  "iteration.sim" is always used (if found) to index the sampled
  values. If this file is not found the sampled values are indexed from
  1 to the sample size. If "mixture.sim" appeares here, only the column
  with number of mixture components is converted to the coda object.}
  \item{data.frames}{a vector of strings giving the names of data.frames
  that are to be converted to \code{coda} objects.
  }
  \item{variant}{a variant of \code{bayessurvreg} function used to
    generate sampled values. This argument is only used to identify
    appropriate files when \code{files} argument is missing.

    Currently only 1 is supported to cooperate with \code{\link{bayessurvreg1}}.
  }
  \item{dir}{string giving the directory where it will be searched for the files with 
    sampled values.}
  \item{start}{the first row (possible header does not count) from the files with the
    sampled values that will be converted to coda objects.}
  \item{end}{the last row from the files with the sampled values that
    will be converted to coda objects. If missing, it is the last row in
    files.}
  \item{thin}{additional thinning of sampled values (i.e. only every
    \code{thin} value from \code{files} and \code{data.frames} is considered).}
  \item{header}{TRUE or FALSE indicating whether the files with the sampled values
    contain also the header on the first line or not.}
  \item{chain}{parameter giving the number of the chain if parallel
    chains were created and sampled values stored in data.frames further
    stored in lists(). If \code{missing}, data.frames are not assumed to
    be stored in lists.}
}
\value{
  A list with \code{mcmc} objects. One object per file or data.frame.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}
\examples{
## *** illustration of usage of parameters 'data.frames' and 'chain' ***
## *********************************************************************
## Two parallel chains with four variables, stored in data.frames
## data.frames are further stored in lists
library(coda)

group1 <- list();  group2 <- list();  group3 <- list()
    ## first chain of first two variables:
  group1[[1]] <- data.frame(var1 = rnorm(100, 0, 1), var2 = rnorm(100, 5, 4))
    ## second chain of first two variables:   
  group1[[2]] <- data.frame(var1 = rnorm(100, 0, 1), var2 = rnorm(100, 5, 4))
    ## first chain of the third variable:
  group2[[1]] <- data.frame(var3 = rgamma(100, 1, 1))
    ## second chain of the third variable:
  group2[[2]] <- data.frame(var3 = rgamma(100, 1, 1))
    ## first chain of the fourth variable:
  group3[[1]] <- data.frame(var4 = rbinom(100, 1, 0.4))
    ## second chain of the fourth variable:
  group3[[2]] <- data.frame(var4 = rbinom(100, 1, 0.4))  

  ## Create mcmc objects for each chain separately
mc.chain1 <- files2coda(data.frames = c("group1", "group2", "group3"), chain = 1)
mc.chain2 <- files2coda(data.frames = c("group1", "group2", "group3"), chain = 2)

  ## Create mcmc.list to represent two parallel chains
mc <- mcmc.list(mc.chain1, mc.chain2)
rm(mc.chain1, mc.chain2)


## *** illustration of usage of parameter 'data.frames' without 'chain' ***
## ************************************************************************
## Only one chain for four variables was sampled and stored in three data.frames
    ## chain of first two variables:
  group1 <- data.frame(var1 = rnorm(100, 0, 1), var2 = rnorm(100, 5, 4))
    ## chain of the third variable:  
  group2 <- data.frame(var3 = rgamma(100, 1, 1))
    ## chain of the fourth variable:
  group3 <- data.frame(var4 = rbinom(100, 1, 0.4))                        

  ## Create an mcmc object
mc <- files2coda(data.frames = c("group1", "group2", "group3"))
}
\keyword{connection}
