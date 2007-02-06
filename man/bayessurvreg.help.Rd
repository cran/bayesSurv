\name{bayessurvreg.help}
\alias{bayessurvreg.help}
\alias{bayessurvreg.design}
\alias{bayessurvreg.checknsimul}
\title{
   Helping function for Bayesian survival regression models.
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of `bayessurvreg' functions to make them more
  readable for the programmer.
}
\usage{
bayessurvreg.design(m, formula, random, data, transform, dtransform)

bayessurvreg.checknsimul(nsimul)
}
\arguments{
  \item{m, formula, random, data, transform, dtransform}{}
  \item{nsimul}{a list}
}
\value{
  Some lists.
}
\section{Value for bayessurvreg.design}{
  A~list with the following components:
  \describe{
  \item{n}{number of observations (in the case of bivariate data, this
    is a~number of single observations, i.e. \eqn{2\times \mbox{sample
	size}}{2*sample size}) included in the dataset}
  \item{ncluster}{number of clusters included in the dataset. In the
    case of bivariate data this is equal to the number of bivariate
    observations. If there are no random effects included in the model
    and if the observations are not bivariate then \code{ncluster = n}}
  \item{nwithin}{a~vector of length equal to \code{ncluster} with
    numbers of observations within each cluster. In the case of
    bivariate observations this is a~vector filled with 2's, if there are
    no random effects and if the observations are not bivariate then
    this is a~vector filled with 1's}
  \item{nY}{number of columns in the response matrix \eqn{Y}{Y}. This is
    equal to 2 if there are no interval-censored observations and equal
    to 3 if there is at least one interval censored observation in the
    dataset}
  \item{nX}{number of columns in the design matrix \eqn{X}{X}. Note that
    the matrix \eqn{X}{X} contains covariates for both fixed and random
    effects}
  \item{nfixed}{number of fixed effects involved in the model. Note that
    possible intercept is always removed from the model}
  \item{nrandom}{number of random effects in the model, possible random
    intercept included}
  \item{randomInt}{\code{TRUE}/\code{FALSE} indicating whether the
    random intercept is included in the model}
  \item{Y}{response matrix. Its last column is always equal to the
    status indicator (1 for exactly observed event times, 0 for
    right-censored observations, 2 for left-censored observations, 3 for
    interval-censored observations).}
  \item{X}{design matrix containing covariates}
  \item{Yinit}{response matrix extracted from \code{formula} using \code{\link[stats]{model.extract}}}
  \item{Xinit}{design matrix extracted from \code{formula} using
    \code{\link[stats]{model.matrix}} function}
  \item{cluster}{a~vector of length \code{n} with identifications of
    clusters (as given by \code{cluster} in \code{formula})}
  \item{indb}{a~vector of length \code{nX} identifying fixed and random
    effects. \code{indb[j] = -1} if the \eqn{j}{j}th column of matrix
    \eqn{X}{X} is a fixed effects. it is equal to \eqn{l}{l} if the
    \eqn{j}{j}th column of matrix \eqn{X}{X} corresponds to the
    \eqn{l}{l}th random effect (in C++ indexing)}
  \item{rnames.X}{row names of \code{Xinit}}
  \item{names.random}{column names of the \eqn{X}{X} matrix corespning
    to the random effects. If there is the random intercept in the
    model, the first component of this vector is equal to "(Intercept)"}
  \item{factors}{???}
  \item{n.factors}{number of factor covariates in the model formula}
  \item{n.in.factors}{???}
  }  
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arnost Komarek} \email{arnost.komarek@med.kuleuven.be}
}
\keyword{internal}

