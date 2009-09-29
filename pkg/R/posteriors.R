#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[posteriors.R] by DSB Fre 04/07/2008 11:40 (CEST) on daniel@puc.home>
##
## Description:
## Extract posterior model probability estimates from BayesMfp object.
##
## History:
## 04/07/2008   copy from thesis function collection.
#####################################################################################

`posteriors` <-
function (BayesMfpObject,
                        ind = 1         # ind = 1 means normalized posteriors, ind = 2 sampling freqs
                        )
{
    sapply (BayesMfpObject, function (one) one[["posterior"]] [ind])
}

