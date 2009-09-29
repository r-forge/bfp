#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[getLogPrior.R] by DSB Die 15/09/2009 15:57 (CEST)>
##
## Description:
## Extract log prior from a model.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 29/11/2008   update for new model prior option
#####################################################################################

getLogPrior <- function (x       # a valid BayesMfp-Object of length 1 (otherwise only first element recognized)
                         )
{
    if(identical(attr(x, "priorSpecs")$modelPrior,
                 "sparse"))
    {
        x <- x[1]
        powers <- x[[1]]$powers
        maxs <- attr (x, "shiftScaleMax")[, 3]
        cards <- attr (x, "shiftScaleMax")[, 4]

        logVals <- numeric(length=length(powers))
        
        for (i in seq_along (powers))
        {
            deg <- length (powers[[i]])
            logVals[i] <- - lchoose(cards[i] - 1 + deg, deg) - log1p(maxs[i])
        }

        return (sum (logVals) - max (attr (x, "indices")$uc) * log (2))
    }
    else
    {
        return(- max (attr (x, "indices")$uc) * log (2))
    }
}
