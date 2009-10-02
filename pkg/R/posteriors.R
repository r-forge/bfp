#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[posteriors.R] by DSB Fre 02/10/2009 11:49 (CEST)>
##
## Description:
## Extract posterior model probability estimates from BayesMfp object.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 02/10/2009   checks, in particular if there are two estimates for the models when ind=2 
#####################################################################################

`posteriors` <-
function (BayesMfpObject,
                        ind = 1         # ind = 1 means normalized posteriors, ind = 2 sampling freqs
                        )
{
    stopifnot(ind %in% c(1, 2))
    if(ind == 2){
        ## does the first model have two probability estimates?
        stopifnot(identical(length(BayesMfpObject[[1]]$posterior),
                            2L))
    }

    ## now extract the required one from all models
    sapply (BayesMfpObject, function (one) one[["posterior"]] [ind])
}

