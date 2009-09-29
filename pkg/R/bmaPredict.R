#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[bmaPredict.R] by DSB Don 13/11/2008 18:18 (CET) on daniel@puc.home>
##
## Description:
## BMA prediction for new data points.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 05/09/2008   correct: centering of x matrix, check ranges
##              numerical cancellations?! (in the summation)
## 10/11/2008   don't ask for the response when creating the model matrix
## 13/11/2008   use new internal function for creating model matrix for newdata
#####################################################################################

## this is not a predict method for BmaSamples!
## it is a different predict "method" for BayesMfp

bmaPredict <- function (BayesMfpObject, # models over which to average the predictions
                        postProbs = posteriors (BayesMfpObject),
                                        # vector of posterior probabilites that will be normalized within
                        newdata,        # new data as data.frame
                        ...             # unused
                        )
{
    if (length (postProbs) != length (BayesMfpObject))
        stop ("postProbs has wrong length")

    tempX <- constructNewdataMatrix(BayesMfpObject=BayesMfpObject,
                                    newdata=newdata)

    ## and compute mean fit with respective original posterior coefficient mode for every model
    fitmat <- matrix (nrow = length (BayesMfpObject), ncol = nrow (tempX))
    for (i in seq_along (BayesMfpObject)){
        tempMod <- BayesMfpObject[i]

        post <- getPosteriorParms (BayesMfpObject[i])
        
        attr (tempMod, "x") <- tempX
        attr(tempMod, "xCentered") <- scale(tempX, center=TRUE, scale=FALSE)
        
        tempDesign <- getDesignMatrix (tempMod)

        fitmat[i,] <- tempDesign %*% post$mStar
    }

    ## average with probabilites as weights
    fitmat <- fitmat * (postProbs / sum (postProbs))

    return (colSums (fitmat))
}
