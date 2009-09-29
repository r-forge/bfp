#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Course:
##         from the Department of Statistics, University of Munich
## Time-stamp: <[shrinkage.R] by DSB Mit 05/11/2008 13:42 (CET) on daniel@puc.home>
##
## Description:
## test the shrinkage sampling
##
## History:
## 05/11/2008   file creation
#####################################################################################

library(bfp)

n <- 1e+5
R2 <- 0.7
nObs <- 100
p <- 20
alpha <- 3.5

tVec <- bfp:::rshrinkage(n=n,
                                    R2=R2,
                                    nObs=nObs,
                                    p=p,
                                    alpha=alpha)
histRet <- hist(tVec,
                nclass=50,
                prob=TRUE)

## compare with unn density:
grid <- histRet$mids

vals <- (1 - grid)^((p + alpha - 2) / 2 - 1) * (1 - R2 * grid)^(-(nObs - 1) / 2)

## scale
vals <- vals / max(vals) * max(histRet$density)

lines(grid,
      vals,
      col=2)
