## Start

library (Hmisc)
library(bfp)
data(ozone)

describe (ozone)
names (ozone)

symnum (cor (subset (na.omit (ozone), select = - weekday))) # correlations between two variables
ozoneFull <- na.omit (subset (ozone, select = - tempElMonte))

workday <- ozoneFull$weekday
levels (workday)
levels (workday) <- c (rep ("workday", 5), rep ("weekend", 2))
ozoneFull$workday <- workday

## day of the year

library (doBy)
ozoneFull <- orderBy (~ month + day, ozoneFull)
ozoneFull$dayOfYear <- seq_len (nf <- nrow (ozoneFull))

set.seed (321)                          # determine test and training set
whichTest <- sample (nf, 30, replace = FALSE)
ozoneTest <- ozoneFull [whichTest, ]
ozoneTraining <- ozoneFull[-whichTest, ]

ozoneTest <- orderBy (~month + day, ozoneTest)
testDates <- with (ozoneTest, paste (month, day, sep = "/"))
testDates <- matrix (testDates, nrow = 3, ncol = 10, byrow = TRUE)


debug(BayesMfp)
## run sampler
set.seed (432)
system.time (
             ozoneModels <- BayesMfp (hourAverageMax ~
                                      bfp (dayOfYear) +
                                      bfp (pressure500Height) +
                                      bfp (windSpeed) +
                                      bfp (humidity) +
                                      bfp (tempSandburg) +
                                      bfp (inversionBaseHeight) +
                                      bfp (pressureGradientDaggett) +
                                      bfp (inversionBaseTemp) +
                                      bfp (visibility),
                                      data = ozoneTraining,
                                      nModels = 3000,
                                      method = "sampling",
                                      chainlength = 2e+5
                                      )
             )

sumOzoneModels <- summary (ozoneModels)
str (sumOzoneModels)

save (ozoneModels, sumOzoneModels, file = "ozoneModels.RData")
load ("ozoneModels.RData")

sumOzoneModels$dataframe[1:50,]

## test: add workday variable
set.seed (567)
system.time (
             ozoneModelsWork <- BayesMfp (hourAverageMax ~
                                      bfp (dayOfYear) +
                                      bfp (pressure500Height) +
                                      bfp (windSpeed) +
                                      bfp (humidity) +
                                      bfp (tempSandburg) +
                                      bfp (inversionBaseHeight) +
                                      bfp (pressureGradientDaggett) +
                                      bfp (inversionBaseTemp) +
                                      bfp (visibility) +
                                          uc (workday),
                                      data = ozoneTraining,
                                      nModels = 3000,
                                      method = "sampling",
                                      chainlength = 2e+5
                                      )
             )
save (ozoneModelsWork, file = "ozoneModelsWork.RData")


## MAP model
mapModel <- ozoneModels[1]
summary (mapModel)

myLty <- c (1, 2, 2, 3, 3)
myLwd <- 2

post <- getPosteriorParms (mapModel)
plotCurveEstimate (mapModel, "dayOfYear", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[0]^2~paste("(", x[0], ",",~alpha, ",", ~p[0],")")),
                   xlab = expression (z[0]),
                   post = post,
                   legendPos = "topright"
                   )
plotCurveEstimate (mapModel, "humidity",  lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[6]^1~paste("(", x[6], ",",~alpha, ",", ~p[6],")")),
                   xlab = expression (z[6]),
                   post = post
                   )
plotCurveEstimate (mapModel, "tempSandburg", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[7]^1~paste("(", x[7], ",",~alpha, ",", ~p[7],")")),
                   xlab = expression (z[7]),
                   post = post
                   )
plotCurveEstimate (mapModel, "inversionBaseHeight", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[8]^1~paste("(", x[8], ",",~alpha, ",", ~p[8],")")),
                   xlab = expression (z[8]),
                   post = post
                   )
plotCurveEstimate (mapModel, "pressureGradientDaggett", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[9]^1~paste("(", x[9], ",",~alpha, ",", ~p[9],")")),
                   xlab = expression (z[9]),
                   post = post
                   )
plotCurveEstimate (mapModel, "visibility", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[10]^2~paste("(", x[10], ",",~alpha, ",", ~p[10],")")),
                   xlab = expression (z[10]),
                   post = post,
                   legendPos = "topright"
                   )
mapPredictions <- predict (mapModel, newdata = ozoneTest)
mapRMSE <- sqrt (mean ((mapPredictions - ozoneTest$hourAverageMax)^2))

## BMA

set.seed (456)
system.time(
            BmaModelSamples <- BmaSamples (ozoneModels, gridSize = 0)
            )
save (BmaModelSamples, file = "BmaOzoneModels.RData")
load ("BmaOzoneModels.RData")

sumBma <- summary (BmaModelSamples)
print (sumBma, table = FALSE)

hist (BmaModelSamples$fixed, nclass = 100)
hist (BmaModelSamples$sigma2, nclass = 100)
hist(BmaModelSamples$shrinkage, nclass=100)

plotCurveEstimate (BmaModelSamples, "dayOfYear")
plotCurveEstimate (BmaModelSamples, "humidity")
plotCurveEstimate (BmaModelSamples, "visibility")

plotCurveEstimate (BmaModelSamples, "tempSandburg")
x11 ()
tempSamples <- BmaModelSamples$bfp[["tempSandburg"]]
matplot (attr (tempSamples, "scaledGrid"), t (as.vector (BmaModelSamples$fixed) + tempSamples), type = "l",
         ylim = c (-1000, 1000))


bmaPredictions <- bmaPredict (ozoneModels, newdata = ozoneTest)
bmaRMSE <- sqrt (mean ((bmaPredictions - ozoneTest$hourAverageMax)^2))
bmaRMSE

## mfp?
library (mfp)
ozoneMfp <- mfp (hourAverageMax ~
                         fp (dayOfYear) +
                         fp (pressure500Height) +
                         fp (windSpeed) +
                         fp (humidity) +
                         fp (tempSandburg) +
                         fp (inversionBaseHeight) +
                         fp (pressureGradientDaggett) +
                         fp (inversionBaseTemp) +
                         fp (visibility),
                         data = ozoneTraining
                         )
ozoneMfp
str (ozoneMfp)
ozoneMfpLm <- lm(ozoneMfp$formula, data = ozoneTraining)
mfpPredictions <- predict (ozoneMfpLm, newdata = ozoneTest)
mfpRMSE <- sqrt (mean ((mfpPredictions - ozoneTest$hourAverageMax)^2))

sumMfp <- summary (ozoneMfp)
str (sumMfp)

ozoneMfpModel <- ozoneModels[1]
ozoneMfpModel[[1]][c ("powers", "ucTerms")] <-
    list (
          powers = list (
          dayOfYear = c (1, 1),
          humidity = 1,
          inversionBaseHeight = 1,
          inversionBaseTemp = 1,
          pressure500Height = 1,
          pressureGradientDaggett = 3,
          tempSandburg = c (0.5, 2),
          visibility = c (-0.5, -0.5),
          windSpeed = 1
          ),
          ucTerms = integer (0)
          )
findModel (ozoneMfpModel, ozoneModels)

mean (as.data.frame (ozoneModels)$logMargLik < getLogMargLik (ozoneMfpModel))
exp (getLogMargLik (ozoneMfpModel) + getLogPrior (ozoneMfpModel))

## third best Casella & Moreno model?

casella <- lm (hourAverageMax ~ humidity + I (windSpeed^2) + I (tempSandburg^2) + I (pressureGradientDaggett^2) +
               factor (month):visibility + pressure500Height:tempSandburg + windSpeed:visibility +
               humidity:inversionBaseHeight,
               data = ozoneTraining)
summary (casella)
casellaPredictions <- predict (casella, newdata = ozoneTest)
casellaRMSE <- sqrt (mean ((casellaPredictions - ozoneTest$hourAverageMax)^2))


## graph comparing prediction performance
predMat <- cbind(casellaPredictions, mfpPredictions, bmaPredictions, mapPredictions)

save (predMat, casellaRMSE, mfpRMSE, bmaRMSE, mapRMSE, file = "ozonePrediction.RData")
load ("ozonePrediction.RData")

myPch <- c (3, 4, 19, 19)
myCol <- c (rep("black", 2), "grey", "black")

table (ozoneTest$hourAverageMax)
x <- jitter (ozoneTest$hourAverageMax)

matplot (x,
         predMat,
         pch = myPch,
         col = myCol,
         type = "p",
         xlim = (lims <- c (0, max (ozoneTest$hourAverageMax)) * 1.1),
         ylim = lims,
         xlab = "actual 1-hour-average ozone level [ppm]",
         ylab = "predicted value"
         )
rug (x)
abline (0, 1)
legend ("bottomright", pch = myPch, col = myCol, legend = c ("Casella", "mfp", "BMA", "MAP"), bty = "n")




