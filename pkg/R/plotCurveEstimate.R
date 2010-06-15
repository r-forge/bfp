`plotCurveEstimate` <-
function (model, termName, plevel = 0.95, slevel = plevel, plot = TRUE, legendPos = "topleft",
          rug=FALSE, ..., main = NULL)
    UseMethod ("plotCurveEstimate")     # define generic function: BayesMfp and BmaSamples methods following

