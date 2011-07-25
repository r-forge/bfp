\name{getLogGPrior}
\alias{getLogGPrior}
\title{Helper function for glmBayesMfp: Returns the normalized log g prior density...}
\usage{getLogGPrior(gPrior)}
\description{Helper function for glmBayesMfp: Returns the normalized log g prior density}
\details{Returns the normalized log of the prior density on g as an R function.
It is checked that the result is valid, that is that exp(f) integrates to 1 over the
positive real line.}
\value{The R function.}
\keyword{internal}
\arguments{\item{gPrior}{the gPrior argument passed to \code{\link{glmBayesMfp}}}}
