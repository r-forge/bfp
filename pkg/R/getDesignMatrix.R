#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[getDesignMatrix.R] by DSB Son 09/11/2008 19:12 (CET) on daniel@puc.home>
##
## Description:
## Extract the design matrix of the first element (== model) of a BayesMfp object.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 03/09/2008   center the (non-intercept) columns
## 04/09/2008   argument center in getFpTransforms is now used, correct centering of uc columns
## 05/09/2008   use xCentered attribute from BayesMfp object
## 09/11/2008   add parameter to control if fixed columns (intercept) should be in the
##              return matrix
#####################################################################################

getDesignMatrix <- function (x, # a valid BayesMfp-Object of length 1 (otherwise only first element
                                # recognized)
                             fixedColumns=TRUE # return the fixed columns inside the matrix or not?
                             )
{
    full <- attr (x, "x")
    fullCentered <- attr(x, "xCentered")
    
    inds <- attr (x, "indices")
    powers <- x[1][[1]]$powers
    ucSet <- x[1][[1]]$ucTerms

    nFix <- length (inds$fixed)         
    stopifnot(identical(nFix, as.integer(1)))

    ucColInds <- inds$uc %in% ucSet
    nUc <- sum (ucColInds)

    nFp <- length (unlist (powers))

    ## reserve space for return matrix
    nColumns <-
        if(fixedColumns)
            nFix + nUc + nFp
        else
            nUc + nFp
    
    ret <- matrix (nrow = nrow (full),
                   ncol = nColumns)
    retColnames <- character (ncol (ret))
    col <- 0                            # invariant: already col columns written

    if(fixedColumns)
    {
        ## fixed columns
        new <- full[, inds$fixed, drop = FALSE]
        newInds <- col + seq_along (inds$fixed)
        
        ret[, newInds] <- new
        retColnames[newInds] <- colnames (new)
        
        col <- col + nFix
    }
    
    ## fp part
    for (i in seq_along (inds$bfp)){
        pi <- powers[[i]]
        if (len <- length (pi)) {       # if there is at least one power
            new <- getFpTransforms (full[, inds$bfp[i], drop = FALSE], pi, center=TRUE)
            newInds <- col + seq_along (pi)

            ret[, newInds] <- new
            retColnames[newInds] <- colnames (new)

            col <- col + len
        }
    }

    ## uc part
    if (length (ucSet)){
        new <- fullCentered[, ucColInds, drop = FALSE]
        newInds <- col + seq_len (nUc)

        ret[, newInds] <- new
        retColnames[newInds] <- colnames (new)

        col <- col + nUc
    }

    rownames (ret) <- rownames (full)
    colnames (ret) <- retColnames

    return (ret)
}



