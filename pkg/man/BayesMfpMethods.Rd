\name{BayesMfp Methods}
\alias{as.data.frame.BayesMfp}
\alias{'[.BayesMfp'}
\alias{fitted.BayesMfp}
\alias{predict.BayesMfp}
\alias{print.BayesMfp}
\alias{residuals.BayesMfp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~functions to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the functions
     as.data.frame.BayesMfp, '[.BayesMfp', fitted.BayesMfp, predict.BayesMfp, print.BayesMfp, residuals.BayesMfp do. ~~
}
\usage{
\method{as.data.frame}{BayesMfp}(x, row.names = NULL, optional = FALSE, freq = TRUE, ...)
\method{[}{BayesMfp}(x, ...)
\method{fitted}{BayesMfp}(object, design = getDesignMatrix(object), post
= getPosteriorParms(object, design = design), ...)
\method{predict}{BayesMfp}(object, newdata, ...)
\method{print}{BayesMfp}(x, ...)
\method{residuals}{BayesMfp}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{row.names}{ ~~Describe \code{row.names} here~~ }
  \item{optional}{ ~~Describe \code{optional} here~~ }
  \item{freq}{ ~~Describe \code{freq} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{design}{ ~~Describe \code{design} here~~ }
  \item{post}{ ~~Describe \code{post} here~~ }
  \item{newdata}{ ~~Describe \code{newdata} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the values returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{Daniel Saban\'es Bov\'e}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{internal}
