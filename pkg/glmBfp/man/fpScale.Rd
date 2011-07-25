\name{fpScale}
\alias{fpScale}
\title{Shift and scale a covariate vector (if wished) to have positive and small numbers.}
\usage{fpScale(x, scaling=TRUE)}
\description{Shift and scale a covariate vector (if wished) to have positive and small numbers.}
\details{This function is (almost exactly) copied from package \dQuote{mfp} to be consistent.}
\value{list of \sQuote{shift} and \sQuote{scale} parameters appropriate for this
covariate vector.}
\keyword{internal}
\arguments{\item{x}{the covariate vector}
\item{scaling}{should shifting and scaling be performed (default)?}}
