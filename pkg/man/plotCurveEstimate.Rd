\name{plotCurveEstimate}
\alias{plotCurveEstimate}
\alias{plotCurveEstimate.BayesMfp}
\alias{plotCurveEstimate.BmaSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{generic function for plotting a fractional polynomial curve estimate}
\description{
  Plot a fractional polynomial curve estimate for either a single model
  or a Bayesian model average over BayesMfp objects. Optionally,
  credible intervals and / or bands can be added to the plot.
}
\usage{
plotCurveEstimate(model, termName, plevel = 0.95, slevel = plevel, plot
= TRUE, legendPos = "topleft", rug = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ an object of class \code{\link{BayesMfp}} or \code{\link{BmaSamples}} }
  \item{termName}{ string denoting an fp term as written by summary method }
  \item{plevel}{credibility level for pointwise HPD, NULL means no pointwise HPD}
  \item{slevel}{credibility level for simultaneous credible band (SCB), NULL means no SCB}
  \item{plot}{if FALSE, only return values needed to produce the plot}
  \item{legendPos}{position of coefficient estimates (for BayesMfp) or
    sample size (for BmaSamples) in the plot, NULL suppresses the
    printing}
  \item{rug}{add a rug to the plot? (default: \code{FALSE})}
  \item{\dots}{ further arguments in case of a BayesMfp object (see details) and
    arguments for plotting with \code{\link{matplot}} }
}
\details{
  Further arguments for application on a BayesMfp object:
  \describe{
    \item{grid}{vector of unscaled abscissae, default is a length \code{gridSize} grid
      over observed range}
    \item{post}{list with posterior parameters of the model, which may be
      provided manually to accelerate plotting in a loop}
    \item{gridSize}{default number of grid points used when no
      \code{grid} is supplied}
    \item{numSim}{number of simulations for estimation of the SCB}
  }
}
\value{
  a list of various plotting information:
  \item{original}{grid on the original covariate scale}
  \item{grid}{grid on the transformed scale}
  \item{mode}{mode curve values, only for BayesMfp object}
  \item{median}{pointwise median curve values, only for BmaSamples object}
  \item{mean}{pointwise mean curve values, only for BmaSamples object}
  \item{plower}{lower boundaries for pointwise HPD}
  \item{pupper}{upper boundaries for pointwise HPD}
  \item{slower}{lower boundaries for SCB}
  \item{supper}{upper boundaries for SCB}
  \item{obsVals}{observed values of the covariate on the original scale}
  \item{sampleSize}{sample size underlying the curve estimate, only for BmaSamples object}
  \item{partialResids}{partial residuals}
  \item{transform}{vector of shift and scale parameter}  
}

\seealso{\code{\link{BayesMfp}}, \code{\link{BmaSamples}}}
\examples{
set.seed(19)
x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)
y <- rt (n=15, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100, method="exhaustive")
plotCurveEstimate (test[20], "x2")

testBma <- BmaSamples (test)
plotCurveEstimate (testBma, "x2")

}

\keyword{regression}

