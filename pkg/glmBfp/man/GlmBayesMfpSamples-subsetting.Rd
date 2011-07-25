\name{GlmBayesMfpSamples-subsetting}
\title{Subset method for GlmBayesMfpSamples objects...}
\usage{\S4method{[}{GlmBayesMfpSamples,ANY,missing,missing}(x, i)}
\description{Subset method for GlmBayesMfpSamples objects}
\details{Index the samples to select a subset of samples.}
\alias{[,GlmBayesMfpSamples,ANY,missing,missing-method}
\value{The subset of the same class.}
\note{The function call will fail if any of the saved bfpCurves or ucCoefs
does not have enough samples to be subset by \code{i} !}
\seealso{\code{\linkS4class{GlmBayesMfpSamples}}}
\keyword{methods}
\arguments{\item{x}{the \code{\linkS4class{GlmBayesMfpSamples}} object}
\item{i}{the vector defining the subset of samples}}
