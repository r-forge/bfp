\name{BmaSamples}
\alias{BmaSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function for Bayesian model averaging over multivariable
  fractional polynomial models}
\description{
  Samples from the Bayesian model average over the models in
  \code{object}, which must be a valid \code{\link{BayesMfp}}-object, 
  are drawn. 
}
\usage{
BmaSamples(object, sampleSize = length(object) * 10, postProbs =
posteriors(object), gridList = list(), gridSize = 203, newdata=NULL,
verbose = TRUE)  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ valid BayesMfp object including the models over which
    to average} 
  \item{sampleSize}{sample size }
  \item{postProbs}{vector of posterior probabilites that will be
    normalized within the function }
  \item{gridList}{optional list of appropriately named grid vectors for fp evaluation,
    default is a length (\code{gridSize} - 2) grid per covariate additional to the observed values
    (two are at the endpoints)}
  \item{gridSize}{see above}
  \item{newdata}{new covariate data with exactly the names (and
    preferably ranges) as before} 
  \item{verbose}{should information on sampling progress been printed?}
}

\value{
  Return an object of class \code{BmaSamples}, which is a list with
  various elements that describe the \code{BayesMfp} object over which
  was averaged, model frequencies in the samples, the samples themselves
  etc. Various methods are available, see the examples.
}

\seealso{\code{\link{BayesMfp}}}
\examples{
set.seed(19)
x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)
y <- rt (n=15, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 200, method="exhaustive")
(testBma <- BmaSamples (test))
str (testBma)
summary (testBma)
fitted (testBma)
residuals (testBma)
bmaPredict (test, newdata = list (x1 = rnorm (15), x2 = rbinom (n=15,
size=5, prob=0.2) + 1, x3 = rexp (n=15)))
}
\keyword{regression}