\name{getNullModelInfo}
\alias{getNullModelInfo}
\title{Internal helper function to compute all information about the null model}
\usage{
  getNullModelInfo(family)
}
\arguments{
  \item{family}{the family object including the
  loglikelihood etc.}
}
\value{
  a list including the elements \describe{
  \item{logMargLik}{log marginal likelihood of the null
  model} }
}
\description{
  This function is used by \code{\link{glmBayesMfp}}. It
  still includes a lot of unnecessary code for
  demonstration purposes (especially for the marginal
  likelihood computation).
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{internal}

