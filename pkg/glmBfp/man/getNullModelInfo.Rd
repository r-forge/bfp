\name{getNullModelInfo}
\alias{getNullModelInfo}
\title{Internal helper function to compute all information about the null model...}
\usage{getNullModelInfo(family)}
\description{Internal helper function to compute all information about the null model}
\details{This function is used by \code{\link{glmBayesMfp}}. It still includes a lot of
unnecessary code for demonstration purposes (especially for the marginal likelihood
computation).}
\value{a list including the elements
\describe{
\item{logMargLik}{log marginal likelihood of the null model}
}}
\keyword{internal}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{family}{the family object including the loglikelihood etc.}}
