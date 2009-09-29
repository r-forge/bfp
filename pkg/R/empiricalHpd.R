#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[empiricalHpd.R] by DSB Fre 04/07/2008 11:24 (CEST) on daniel@puc.home>
##
## Description:
## Compute a MC HPD estimate for a scalar parameter.
##
## History:
## 04/07/2008   copy from thesis function collection.
#####################################################################################

empiricalHpd <- function (theta,        # sample vector of parameter
                          level         # credibility level
                          )
{
    M <- length (theta)
    thetaorder <- sort.int (theta, method = "quick")

    alpha <- 1 - level
    maxSize <- ceiling (M * alpha)
    ind <- seq_len (maxSize)

    lower <- thetaorder[ind]
    upper <- thetaorder[M - maxSize + ind]
    size <- upper - lower
    sizeMin <- which.min(size)

    HPDLower <- thetaorder[sizeMin]
    HPDUpper <- thetaorder[M - maxSize + sizeMin]
    return (c (HPDLower, HPDUpper))
}
