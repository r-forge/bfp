\name{sampleGlm}
\alias{sampleGlm}
\title{Produce posterior samples from one GLM}
\usage{
  sampleGlm(object, mcmc = McmcOptions(), estimateMargLik =
  TRUE, gridList = list(), gridSize = 203L, newdata = NULL,
  fixedZ = NULL, marginalZApprox = NULL, verbose = TRUE,
  debug = FALSE, useOpenMP = TRUE)
}
\arguments{
  \item{object}{the \code{GlmBayesMfp} object, from which
  only the first model will be processed (at least for now
  \ldots)}

  \item{mcmc}{MCMC options object with class
  \code{\linkS4class{McmcOptions}}}

  \item{estimateMargLik}{shall the marginal likelihood be
  estimated in parallel? (default)}

  \item{gridList}{optional list of appropriately named grid
  vectors for FP evaluation. Default is length
  (\code{gridSize} - 2) grid per covariate additional to
  the observed values (two are at the endpoints)}

  \item{gridSize}{see above (default: 203)}

  \item{newdata}{new covariate data.frame with exactly the
  names (and preferably ranges) as before (default: no new
  covariate data)}

  \item{fixedZ}{either \code{NULL} (default) or a (single)
  fixed z value to be used, in order to sample from the
  conditional posterior given this z. If \code{object} was
  constructed by the empirical Bayes machinery, this will
  default to the estimated z with maximum conditional
  marginal likelihood.}

  \item{marginalZApprox}{method for approximating the
  marginal density of the log covariance factor z, see
  \code{\link{getMarginalZ}} for the details (default: same
  preference list as in \code{\link{getMarginalZ}})}

  \item{verbose}{should information on computation progress
  be given? (default)}

  \item{debug}{print debugging information? (not default)}

  \item{useOpenMP}{shall OpenMP be used to accelerate the
  computations? (default)}
}
\value{
  Returns a list with the following elements: \describe{
  \item{samples}{an object of S4 class
  \code{\linkS4class{GlmBayesMfpSamples}}}
  \item{coefficients}{samples of all original coefficients
  in the model (nCoefs x nSamples)}
  \item{acceptanceRatio}{proportion of accepted
  Metropolis-Hastings proposals} \item{logMargLik}{if
  \code{estimateMargLik} is \code{TRUE}, this list is
  included: it contains the elements \code{numeratorTerms}
  and \code{denominatorTerms} for the numerator and
  denominator samples of the Chib Jeliazkov marginal
  likelihood estimate,
  \code{highDensityPointLogUnPosterior} is the log
  unnormalized posterior density at the fixed parameter and
  the resulting \code{estimate} and \code{standardError}.}
  }
}
\description{
  Based on the result list from \code{\link{glmBayesMfp}},
  for the first model in the list MCMC samples are
  produced. In parallel to the sampling of coefficients and
  FP curve points, optionally the marginal likelihood of
  the model is estimated with MCMC samples. This provides a
  check of the integrated Laplace approximation used in the
  model sampling.
}
\keyword{models}
\keyword{regression}

