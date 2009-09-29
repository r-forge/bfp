#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[constructNewdataMatrix.R] by DSB Don 13/11/2008 18:08 (CET) on daniel@puc.home>
##
## Description:
## Internal function to construct the model matrix for new data based on the formula
## and scaling info in an existing BayesMfpObject, for passing it to prediction functions. 
##
## History:
## 13/11/2008   file creation
#####################################################################################


constructNewdataMatrix <- function(BayesMfpObject, # from BayesMfp
                                   newdata)        # new data as data.frame
{
    ## correct model matrix to new data
    covariatesOnlyFormula <- attr (BayesMfpObject, "formula")[-2]
    bfpTermNames <- attr (BayesMfpObject, "termNames")$bfp

    newTerms <- terms (covariatesOnlyFormula, data = newdata) # get model matrix
    ret <- model.matrix (newTerms, data = newdata)

    tr <- attr (BayesMfpObject, "shiftScaleMax")
    for (bfp in bfpTermNames)           # and scale columns like for the old data
    {          
        ret[, bfp] <- (ret[, bfp] + tr[bfp, "shift"]) / tr[bfp, "scale"]

        ## check range
        if(any(ret[, bfp] <= 0))
            stop(simpleError(paste("New data for covariate", bfp, "is out of range.")))
    }

    return(ret)
}



