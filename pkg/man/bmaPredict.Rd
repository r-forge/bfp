\name{bmaPredict}
\alias{bmaPredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~functions to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the functions
     bmaPredict do. ~~
}
\usage{
bmaPredict(BayesMfpObject, postProbs = posteriors(BayesMfpObject), newdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BayesMfpObject}{ ~~Describe \code{BayesMfpObject} here~~ }
  \item{postProbs}{ ~~Describe \code{postProbs} here~~ }
  \item{newdata}{ ~~Describe \code{newdata} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The functions are currently defined as
## bmaPredict
function (BayesMfpObject, # models over which to average the predictions
                        postProbs = posteriors (BayesMfpObject),
                                        # vector of posterior probabilites that will be normalized within
                        newdata,        # new data as data.frame
                        ...             # unused
                        )
{
    if (length (postProbs) != length (BayesMfpObject))
        stop ("postProbs has wrong length")

    ## correct model matrix to new data
    oldFormula <- attr (BayesMfpObject, "formula")
    bfpTermNames <- attr (BayesMfpObject, "termNames")$bfp

    newTerms <- terms (oldFormula, data = newdata) # get model matrix
    mf <- model.frame (newTerms, data = newdata)
    tempX <- model.matrix (newTerms, data = mf)

    tr <- attr (BayesMfpObject, "shiftScaleMax")
    for (bfp in bfpTermNames){          # and scale columns like for the old data
        tempX[, bfp] <- (tempX[, bfp] + tr[bfp, "shift"]) / tr[bfp, "scale"]
    }

    ## and compute mean fit with respective original posterior coefficient mode for every model
    fitmat <- matrix (nrow = length (BayesMfpObject), ncol = nrow (tempX))
    for (i in seq_along (BayesMfpObject)){
        tempMod <- BayesMfpObject[i]
        attr (tempMod, "x") <- tempX
        design <- getDesignMatrix (tempMod)
        post <- getPosteriorParms (BayesMfpObject[i])
        fitmat[i,] <- design \%*\% post$mStar
    }

    ## average with probabilites as weights
    fitmat <- fitmat * (postProbs / sum (postProbs))

    return (colSums (fitmat))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}

