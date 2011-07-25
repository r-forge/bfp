\name{evalZdensity}
\alias{evalZdensity}
\title{Evaluate the (negative log) unnormalized marginal z density in a given...}
\usage{evalZdensity(config, object, zValues, conditional=FALSE, debug=FALSE,
    binaryLogisticCorrection=FALSE)}
\description{Evaluate the (negative log) unnormalized marginal z density in a given
model.}
\details{Based on the result list from \code{\link{glmBayesMfp}}, for the first model
in the list, the marginal z density can be evaluated.}
\value{the negative log marginal unnormalized density values at the
\code{zValues}. (Note the words \dQuote{negative}, \dQuote{log}, and
\dQuote{unnormalized} !!!)}
\keyword{internal}
\arguments{\item{config}{the configuration of a single \code{GlmBayesMfp} model. The
null model is not allowed. It is interpreted in the context of
\code{object}.}
\item{object}{the \code{\link{GlmBayesMfp}} object}
\item{zValues}{the z values}
\item{conditional}{return the approximate *conditional* density f(y | z)?
(not default)}
\item{debug}{print debugging information? (not default)}
\item{binaryLogisticCorrection}{should a higher-order correction of the
Laplace approximation be used, which works only (!!!) for binary logistic
regression? (not default)
todo: this user option should later be replaced
todo: by an automatic check if the higher order correction can be applied.}}
