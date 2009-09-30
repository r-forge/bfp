\name{BayesMfp}
\alias{BayesMfp}
\alias{bfp}
\alias{uc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian model inference for multivariable fractional polynomial models}
\description{
  Bayesian model inference for multivariable fractional polynomial
  models is conducted by means of either exhaustive model space
  evaluation or posterior model sampling.
}
\usage{
BayesMfp(formula = formula(data), data = parent.frame(), family =
gaussian, priorSpecs = list(a = 4, modelPrior = "flat"), subset = NULL,
na.action = na.omit, verbose = TRUE, nModels = NULL, method = c("ask",
"exhaustive", "sampling"), chainlength = 10000) 
bfp(x, max = 2, scale = TRUE, rangeVals=NULL)
uc(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{model formula}
  \item{data}{optional data.frame for model variables}
  \item{family}{distribution and link: only gaussian("identity") supported at the moment}
  \item{priorSpecs}{prior specifications, see details}
  \item{subset}{optional subset expression}
  \item{na.action}{default is to skip rows with missing data, no other
    option supported at the moment}
  \item{verbose}{should information on computation progress be given?}
  \item{nModels}{how many best models should be saved? default: 1\% of
    the explored models, 1 would mean MAP}
  \item{method}{which method should be used to explore the  posterior model space?}
  \item{chainlength}{only has effect if sampling has been chosen as
    method}
  \item{x}{variable}
  \item{max}{maximum degree for this FP}
  \item{scale}{use pre-transformation scaling to avoid numerical
    problems?}
  \item{rangeVals}{extra numbers if the scaling should consider values
    in this range}   
  }
\details{
  The formula is of the form
  \code{y ~ bfp (x1, max = 4) + uc (x2 + x3)}, that is, the
  auxiliary functions \code{\link{bfp}} and \code{\link{uc}} must be
  used for defining the fractional polynomial and uncertain fixed form
  covariates terms, respectively. There must be an intercept, and no
  other fixed covariates are allowed.

  The prior specifications are a list:
  \describe{
    \item{a}{hyperparameter for hyper-g prior}
    \item{modelPrior}{choose if a flat model prior (default,
      \code{"flat"}) should be used or a model prior favoring explicitly
      sparse models (\code{"sparse"})}
  }

  If \code{method = "ask"}, the user is prompted with the maximum
  cardinality of the model space and can then decide whether to use
  posterior sampling or the exhaustive model space evaluation.
}
\value{
  Returns an object of class \code{BayesMfp} that inherits from list. It
  is essentially a list of models with certain attributes. Each model
  has the following components:

  \item{powers}{a list of numeric vectors, where each vector contains
    the powers of the covariate that its name denotes.}
  \item{ucTerms}{an integer vector of the indices of uncertain fixed
    form covariates that are present in the model.}

  Subsetting the object
  with \code{\link{'[.BayesMfp'}} returns again a \code{BayesMfp} object
  with the same attributes, which are

  \item{numVisited}{the number of models that have been visited}
  \item{inclusionProbs}{BMA inclusion probabilities for all uncertain
    covariates}
  \item{logNormConst}{the (estimated) log normalizing constant \eqn{f
      (D)}}
  \item{chainlength}{length of the Markov chain, only present if \code{method = "sampling"}}
  \item{call}{the original call}
  \item{formula}{the formula by which the appropriate untransformed
    design matrix can be extracted}
  \item{x}{the shifted and scaled design matrix for the data}
  \item{y}{the response vector}
  \item{yMean}{the mean of the response values}
  \item{SST}{sum of squares total}
  \item{indices}{a list with components that describe the positions of
    uncertain covariate groups, fractional polynomial terms and fixed
    variables in the design matrix}
  \item{termNames}{a list of character vectors containing the names of
    uncertain covariate groups, fractional polynomial terms and fixed
    variables}
  \item{shiftScaleMax}{matrix with 4 columns containing preliminary
    transformation parameters, maximum degrees and cardinalities of the
    powersets of the fractional polynomial terms}
  \item{priorSpecs}{the utilized prior specifications}
}

\note{\code{logNormConst} may be unusable due to necessary conversion
  from long double to double! Various methods for posterior summaries are available.
}

\seealso{
  \code{\link{BmaSamples}}, \code{\link{'[.BayesMfp'}}
}

\examples{
set.seed(19)
x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)
y <- rt (n=15, df=2)

(test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100, method="exhaustive"))
summary (test)
fitted (test)
residuals (test)
predict (test, newdata = list (x1 = rnorm (15), x2 = rbinom (n=15,
size=5, prob=0.2) + 1, x3 = rexp (15)))
}

\keyword{regression}
