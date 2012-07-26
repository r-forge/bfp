\name{sampleBma}
\alias{sampleBma}
\title{Produce posterior samples from a Bayesian model average over GLMs}
\usage{
  sampleBma(object, mcmc = McmcOptions(),
    postProbs = posteriors(object), nMargLikSamples = NULL,
    verbose = TRUE, ...)
}
\arguments{
  \item{object}{valid \code{GlmBayesMfp} object containing
  the models over which to average}

  \item{mcmc}{MCMC options object with class
  \code{\linkS4class{McmcOptions}}, specifying the number
  of required BMA samples (via \code{sampleSize(mcmc)}),
  and the burn-in and thinning parameters applied to each
  model (see above).}

  \item{postProbs}{vector of posterior probabilites (will
  be normalized within the function) for the weighting of
  the models in \code{object} (defaults to the normalized
  posterior probabilities)}

  \item{nMargLikSamples}{If this is non-\code{NULL}, it
  specified the number of samples used for the marginal
  likelihood estimate for each model (see above).}

  \item{verbose}{should information on computation progress
  be given? (default)}

  \item{\dots}{optional further arguments already available
  for sampling from a single model: \code{gridList},
  \code{gridSize}, \code{newdata}, \code{weights},
  \code{marginalZApprox}, \code{debug}, \code{useOpenMP}.
  See \code{\link{sampleGlm}} for the meanings.}
}
\value{
  The result is a list with the following elements:
  \describe{ \item{modelData}{data frame containing the
  result from the \code{as.data.frame} function, and in
  addition BMA probabilities, BMA frequencies in the
  sample, acceptance ratios of the MCMC runs and optionally
  marginal likelihood estimates / standard errors.}
  \item{samples}{an object of S4 class
  \code{\linkS4class{GlmBayesMfpSamples}} containing the
  samples from the BMA.} }
}
\description{
  Based on the result list from \code{\link{glmBayesMfp}},
  sample from the Bayesian model average (BMA) over the
  models contained in this list.
}
\details{
  The sampling is based on MCMC. Therefore, instead of only
  specifying the required number of samples and the model
  probabilities, one also needs to specify the burn-in
  length and the thinning parameter, which will be applied
  to every model from which at least one sample is included
  in the average. Alternatively, you can ask for MCMC
  marginal likelihood estimates for all models in the list.
  Then at least \code{nMargLikSamples} will be produced for
  each model, whether included in the BMA sample or not.
}
\keyword{models}
\keyword{regression}

