\name{as.data.frame.GlmBayesMfp}
\alias{as.data.frame.GlmBayesMfp}
\title{Convert a GlmBayesMfp object into a data frame}
\usage{
  \method{as.data.frame}{GlmBayesMfp} (x, row.names = NULL,
  ..., freq = TRUE)
}
\arguments{
  \item{x}{valid \code{\link{GlmBayesMfp}} object}

  \item{row.names}{optional rownames (default is to keep
  the names of the \code{\link{GlmBayesMfp}} list)}

  \item{\dots}{unused}

  \item{freq}{should empirical frequencies of the models in
  the sampling path be given? (default)}
}
\value{
  The data frame with the following columns: \describe{
  \item{posterior}{the posterior model probabilities}
  \item{logMargLik}{the log marginal likelihood of the
  models} \item{logPrior}{the log prior probabilities of
  the models} } Additionally, for each uncertain fixed form
  covariates a column with the inclusion status, and for
  each fractional polynomial a column with the powers are
  returned.
}
\description{
  Convert a GlmBayesMfp object into a data frame
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\seealso{
  \code{\link{glmBayesMfp}}
}
\keyword{methods}

