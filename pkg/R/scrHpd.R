#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[scrHpd.R] by DSB Fre 04/07/2008 11:40 (CEST) on daniel@puc.home>
##
## Description:
## Calculate a series of simultaneous credible bounds from a samples matrix.
##
## History:
## 04/07/2008   copy from thesis function collection.
#####################################################################################

## all methods assume that samples is a m by n matrix where
## m is the number of samples and n the number of parameters
## hence each sample is a row in the matrix samples

## scrHpd calculates a series of simultaneous credible bounds,
## minimizing the absolute distances to the mode vector at each gridpoint

scrHpd <- function(samples,          # sample matrix
                   mode = apply (samples, 2, median), # mode vector of length = ncol (samples)
                   level = 0.95,     # credibility level
                   grid = ncol (samples) # number of elements in each sample vector
                   )
{
    n <- ncol(samples)
    m <- nrow(samples)

    if (n != length (mode))
        stop ("mode vector must have same length as samples matrix!")

    distance <- abs (sweep (samples, 2, mode)) # absolute distance from mode vector

    k <- floor (level * m)
    ## Calculates a simultaneous (k/m)*100% credible band
    ## using the ranks approach

    rankdistance <- apply (distance, 2, rank) # colwise (= elementwise) ranks of distances

    tstari <- apply (rankdistance, 1, max) # maximum ranks in each sample
    ordtstari <- sort.int (tstari, method = "quick") # and sort the maximum ranks

    tstar <- ordtstari[k]               # the required rank

    selectMat <- rankdistance <= tstar

    ret <- matrix (nrow = 2, ncol = n)
    for (i in seq_len (n)){
        ret[, i] <- range (samples[selectMat[,i], i])
    }

    return (ret)
}
