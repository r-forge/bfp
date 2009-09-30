\name{findModel}
\alias{findModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{find a specific fractional polynomial model in a
  BayesMfp object}
\description{
  Returns the index of the wished model if it is present in the model
  list, and otherwise returns \code{NA}.
}
\usage{
findModel(model, BayesMfpObject)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{the specific model, a list with entries \code{powers} and
  \code{ucTerms}}
  \item{BayesMfpObject}{an object of class \code{\link{BayesMfp}}}
}
\details{
  see \code{\link{BayesMfp}} for the description of a model. 
}
\value{
  index of \code{model} in \code{BayesMfpObject} if it is present in the
  model list, otherwise \code{NA}.
}

\note{The searched model must have exactly the same construction as the
  models in \code{BayesMfpObject}. See the example for the recommended use.
}

\examples{
set.seed(92)
x1 <- rnorm (15)
x2 <- rbinom (n=15, size=20, prob=0.6)
x3 <- rexp (15)
y <- rt (15, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100, method="sampling")

myModel <- test[[1]] # copy one model 
myModel$powers["x2"] <- c (1, 2)# and modify it!
findModel (myModel, test)
}

\keyword{utilities}