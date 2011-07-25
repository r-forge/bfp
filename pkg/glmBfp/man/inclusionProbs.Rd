\name{inclusionProbs}
\alias{inclusionProbs}
\title{Compute posterior inclusion probabilites based on GlmBayesMfp object...}
\usage{inclusionProbs(GlmBayesMfpObject,
    postProbs=posteriors(GlmBayesMfpObject, type = "normalized"))}
\description{Compute posterior inclusion probabilites based on GlmBayesMfp object}
\details{Compute (model averaged) posterior inclusion probabilites for the uncertain
variables (including FP variables) based on a GlmBayesMfp object.}
\value{the resulting inclusion probabilities vector}
\keyword{utilities}
\arguments{\item{GlmBayesMfpObject}{the GlmBayesMfp object}
\item{postProbs}{the vector of posterior model probabilities, defaults to
the normalized probabilities in \code{GlmBayesMfpObject}}}
