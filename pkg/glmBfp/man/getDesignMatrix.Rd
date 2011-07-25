\name{getDesignMatrix}
\alias{getDesignMatrix}
\title{Construct the design matrix for a given bfp GLM model...}
\usage{getDesignMatrix(modelConfig=object[[1]]$configuration, object,
    fixedColumns=TRUE)}
\description{Construct the design matrix for a given bfp GLM model}
\details{This is an internal function to construct the (centered) design matrix for a given
bfp GLM model.}
\value{The design matrix, where the non-fixed part is columnwise centered (that is, the
colmeans are zero).}
\keyword{internal}
\keyword{utilities}
\arguments{\item{modelConfig}{the model configuration list which must have elements
\dQuote{powers} and \dQuote{powers}. Defaults to the configuration of the first element of}
\item{object}{the \code{GlmBayesMfp} object, which is needed because it contains
the covariates matrix and indices vector}
\item{fixedColumns}{return the fixed columns inside the matrix (default) or not?}}
