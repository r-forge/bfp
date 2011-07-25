\name{computeModels}
\alias{computeModels}
\title{Compute model information for a given list of model configurations and glmBayesMfp output.}
\usage{computeModels(configurations, object, verbose=length(configurations) >
    100L, debug=FALSE)}
\description{Compute model information for a given list of model configurations and glmBayesMfp output.}
\details{If we want to compute the marginal likelihood and information necessary for
generating posterior samples for new models not encountered in the model search
done by \code{\link{glmBayesMfp}}, this function can be used: Provide it with the
models \code{configurations} to be interpreted in the context of the \code{object}
of class \code{\link{GlmBayesMfp}}. The result is again of the latter class, but contains
only the new models (similarly as the whole model space would consist of these and an
exhaustive search would have been conducted).}
\value{The \code{\link{GlmBayesMfp}} object with the new models. This can directly
be used as input for \code{\link{sampleGlm}}.}
\keyword{models}
\keyword{regression}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{configurations}{list of the model configurations}
\item{object}{the \code{\link{GlmBayesMfp}} object}
\item{verbose}{be verbose? (default: only for more than 100 configurations)}
\item{debug}{be even more verbose and echo debug-level information? (not by default)}}
