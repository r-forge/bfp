#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[glmBayesMfp.R] by DSB Fre 29/07/2011 15:12 (CEST)>
##
## Description:
## Main user interface for Bayesian inference for fractional polynomials in generalized linear
## models. 
##
## History:
## 26/10/2009   file creation: copy and modify the old BayesMfp.R
## 27/10/2009   first only allow the binomial family, and write the code in an
##              extensible way, so that it later can be extended easily when binomial works
##              (11/12/2009 note: actually we now have generality for all GLM families.)
## 06/11/2009   include group sizes
## 09/11/2009   add Gauss Hermite quantiles and weights
## 18/11/2009   progress towards first testable version of the marginal likelihood
##              approximation,
##              use roxygen after cleaning the package from old hyper-g stuff.
## 01/12/2009   do not coerce model response to type double automatically,
##              because it could be a factor in the binomial case.
## 11/12/2009   rewrite passing of info to C++ and storage of attributes info,
##              to make it clearer and facilitate reuse in other functions, e.g. sampleGlm.
## 06/01/2010   extend getFamily to return self-crafted "simulate" function in the
##              return list. This is used e.g. by "sampleGlm".
## 12/02/2010   for clarity, split off helper functions into their separate files,
##              so we need to include them in the preamble
## 17/02/2010   split off null model information computation, and include the returned
##              list in the things passed to C++ (instead of only the log marginal
##              likelihood for the null model)
## 15/03/2010   also pass weights to C++ via the family list,
##              new g-prior class is expected in the prior list
## 12/04/2010   do not coerce totalNumber to integer but to double, as it is done in C++
##              in the same manner and keeps the number space large enough.
## 14/04/2010   add "useBfgs" and "largeVariance" options.
## 17/05/2010   useBfgs=FALSE ("optimize") is now the default because it is more
##              robust than Bfgs. Add "useOpenMP" option, which makes it easier
##              to switch than setting the environment variable OMP_NUM_THREADS.
## 21/05/2010   be more careful when passing "fpnames" to C++
## 25/05/2010   add "nObs" to "data" attribute of return list.
## 08/07/2010   add an empirical Bayes option which ranks the models in terms of
##              approximate *conditional* marginal likelihoods
## 29/07/2010   add the new option to get a better Laplace approximation in the
##              case of binary logistic regression.
## 08/07/2011   add the new modelPrior option "dependent"
## 29/07/2011   now "higherOrderCorrection"
#####################################################################################

##' @include helpers.R
##' @include formula.R
##' @include fpScale.R
##' @include GPrior-classes.R
##' @include getFamily.R
##' @include getNullModelInfo.R
{}


##' Bayesian model inference for fractional polynomial GLMs
##'
##' Bayesian model inference for fractional polynomial models from the generalized linear model
##' family is conducted by means of either exhaustive model space evaluation or posterior model
##' sampling. The approach is based on analytical marginal likelihood approximations, using
##' integrated Laplace approximation.
##'
##' The formula is of the form
##' \code{y ~ bfp (x1, max = 4) + uc (x2 + x3)}, that is, the
##' auxiliary functions \code{\link{bfp}} and \code{\link{uc}} must be
##' used for defining the fractional polynomial and uncertain fixed form
##' covariates terms, respectively. There must be an intercept, and no
##' other fixed covariates are allowed. All \code{max} arguments of the
##' \code{\link{bfp}} terms must be identical.
##' 
##' The prior specifications are a list:
##' \describe{
##'   \item{gPrior}{A g-prior class object. Defaults to a hyper-g prior. See
##'   \code{\linkS4class{GPrior}} for more information.} 
##'   \item{modelPrior}{choose if a flat model prior (\code{"flat"}), a
##'   model prior favoring sparse models explicitly (default, \code{"sparse"}),
##'   or a dependent model prior (\code{"dependent"}) should be used.}
##' }
##' 
##' If \code{method = "ask"}, the user is prompted with the maximum
##' cardinality of the model space and can then decide whether to use
##' posterior sampling or the exhaustive model space evaluation.
##'
##' Note that if you specify only one FP term, the exhaustive model search
##' must be done, due to the structure of the model sampling algorithm.
##' However, in reality this will not be a problem as the model space will
##' typically be very small.
##' @param formula model formula
##' @param data optional data.frame for model variables (defaults to the parent
##' frame)
##' @param weights optionally a vector of positive weights (if not provided, a
##' vector of one's)
##' @param family distribution and link (as in the glm function)
##' @param phi value of the dispersion parameter (defaults to 1)
##' @param empiricalBayes rank the models in terms of \emph{conditional}
##' marginal likelihood, using an empirical Bayes estimate of g? (not default)
##' Due to coding structure, the prior on g must be given in \code{priorSpecs}
##' although it does not have an effect when \code{empiricalBayes==TRUE}.
##' @param priorSpecs prior specifications, see details
##' @param method which method should be used to explore the  posterior model
##' space? (default: ask the user)
##' @param subset optional subset expression
##' @param na.action default is to skip rows with missing data, and no other
##' option supported at the moment 
##' @param verbose should information on computation progress be given?
##' (default)
##' @param debug print debugging information? (not default)
##' @param nModels how many best models should be saved? (default: 1\% of the
##' total number of (cached) models). Must not be larger than \code{nCache} if
##' \code{method == "sampling"}.
##' @param nCache maximum number of best models to be cached at the same time
##' during the model sampling, only has effect if method = sampling 
##' @param chainlength length of the model sampling chain (only has an effect if
##' sampling has been chosen as method) 
##' @param nGaussHermite number of quantiles used in Gauss Hermite quadrature
##' for marginal likelihood approximation (and later in the MCMC sampler for the
##' approximation of the marginal covariance factor density). If
##' \code{empiricalBayes}, this option has no effect.
##' @param useBfgs Shall the BFGS algorithm be used in the internal maximization
##' (not default)? Else, the default Brent optimize routine is used, which seems
##' to be more robust. If \code{empiricalBayes}, this option has no effect and
##' always the Brent optimize routine is used.
##' @param largeVariance When should the BFGS variance estimate be considered
##' \dQuote{large}, so that a reestimation of it is computed? (Only has an
##' effect if \code{useBfgs == TRUE}, default: 100)
##' @param useOpenMP shall OpenMP be used to accelerate the computations?
##' (default)
##' @param higherOrderCorrection should a higher-order correction of the
##' Laplace approximation be used, which works only for canonical GLMs? (not
##' default) 
##'
##' @aliases glmBayesMfp GlmBayesMfp
##' @return An object of S3 class \code{GlmBayesMfp}.
##' 
##' @keywords models regression
##' @export
glmBayesMfp <-
    function (formula = formula(data), 
              data = parent.frame(),   
              weights,            
              family = gaussian,       
              phi=1,
              empiricalBayes=FALSE,
              priorSpecs =            
              list(gPrior=HypergPrior(), 
                   modelPrior="sparse"),
              method = c ("ask", "exhaustive", "sampling"),
              subset,           
              na.action = na.omit,     
              verbose = TRUE,
              debug=FALSE,
              nModels,
              nCache=1e9,
              chainlength = 1e4,  
              nGaussHermite=20,
              useBfgs=FALSE,
              largeVariance=100,
              useOpenMP=TRUE,
              higherOrderCorrection=FALSE)
{
    ## checks
    stopifnot(is.bool(verbose),
              is.bool(debug),
              is.bool(useBfgs),
              is.bool(empiricalBayes),
              is(priorSpecs$gPrior, "GPrior"),
              is.bool(useOpenMP),
              is.bool(higherOrderCorrection))
    
    ## save call for return object
    call <- match.call()
    method <- match.arg (method)

    ## check and evaluate Gauss Hermite stuff
    nGaussHermite <- as.integer(nGaussHermite)
    gaussHermite <- statmod::gauss.quad(n=nGaussHermite, kind="hermite")
    
    ## evaluate family, this list then also includes the dispersion
    family <- getFamily(family, phi)    

    ## get model prior choice
    priorSpecs$modelPrior <- match.arg(priorSpecs$modelPrior,
                                       choices=c("flat", "sparse", "dependent"))

    
    ## evaluate call for model frame building
    m <- match.call(expand = FALSE)

    ## select normal parts of the call
    temp <- c("", "formula", "data", "weights", "subset", "na.action") # "" is the function name
    m <- m[match(temp, names(m), nomatch = 0)]
   
    ## sort formula, so that bfp comes before uc
    Terms <- if (missing(data))
        terms(formula)
    else
        terms(formula, data = data)

    ## check if intercept is present
    if (! attr(Terms, "intercept"))
        stop(simpleError("there must be an intercept term in the model formula"))

    ## now sort the formula
    sortedFormula <- paste(deparse (Terms[[2]]),
                           "~ 1 +", 
                           paste(sort(attr(Terms, "term.labels")),
                                 collapse = "+"))
    sortedFormula <- as.formula (sortedFormula)

    ## filter special parts in formula: uncertain covariates (uc) and (Bayesian) fractional polynomials (bfp)
    special <- c("uc", "bfp")
    Terms <- if (missing(data))
        terms(sortedFormula, special)
    else
        terms(sortedFormula, special, data = data)

    ucTermInd <- attr (Terms, "specials")$uc  # special indices in original formula (beginning with 1 = response!)
    nUcGroups <- length (ucTermInd)
    bfpTermInd <- attr (Terms, "specials")$bfp
    nFps <- length (bfpTermInd)
    
    ## check if bfp's are present
    if (nFps == 0)
        warning(simpleWarning("no fractional polynomial terms in formula"))
    
    
    ## get vector with covariate entries
    vars <- attr (Terms, "variables")    # language object
    varlist <- eval (vars, env = data)              # list
    covariates <- paste(as.list (vars)[-c(1,2)]) # vector with covariate entries (no list or response or Intercept)

    ## remove bfp() from entries and save the inner arguments
    bfpInner <- varlist[bfpTermInd]     # saved for later use
    covariates[bfpTermInd - 1] <- unlist(bfpInner) # remove bfp( ) from formula; -1 because of reponse column

    ## if ucs are present:
    if (nUcGroups){
        ## remove uc() from entries and save the inner arguments
        ucInner <- unlist(varlist[ucTermInd])
        covariates[ucTermInd - 1] <- ucInner
        
        ## determine association of terms with uc groups
        ucTermLengths <- sapply (ucInner, function (oneUc)
                                 length (attr (terms (as.formula (paste ("~", oneUc))), "term.labels"))
                                 )
        ucTermLengthsCum <- c(0, cumsum (ucTermLengths - 1)) # how much longer than 1, accumulated
        ucTermList <- lapply (seq (along = ucTermInd), function (i) # list for association uc group and assign index
                              as.integer(
                                         ucTermInd[i] - 1 + # Starting assign index
                                         ucTermLengthsCum[i]+ # add lengths from before
                                         0:(ucTermLengths[i]-1) # range for this uc term
                                         )
                              )
    } else {
        ucInner <- ucTermList <- NULL
    }
    ## consistency check:
    stopifnot(identical(length(ucTermList), nUcGroups))

    ## build new formula from the cleaned covariate entries
    newFormula <-                       # is saved for predict method at the end
        update (sortedFormula,
                paste (".~ 1 +",  
                       paste (covariates, collapse = "+")) # only update RHS
                )
    newTerms <- if (missing(data))
        terms(newFormula)
    else
        terms(newFormula, data = data)

    ## build model frame
    m$formula <- newTerms
    m$scale <- m$family <- m$verbose <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    ## build design matrix
    X <- model.matrix (newTerms, m)
    Xcentered <- scale(X, center=TRUE, scale=FALSE)

    ## get and check weights
    weights <- as.vector(model.weights(m))
    if(is.null(weights))
        weights <- rep(1, nrow(X))
    
    if (!is.null(weights) && !is.numeric(weights)) 
        stop(simpleError("'weights' must be a numeric vector"))
    if (!is.null(weights) && any(weights < 0)) 
        stop(simpleError("negative weights not allowed"))
    
    ## get response
    Y <- model.response(m)

    ## initialize (here e.g. the binomial matrix Y case is handled as in 'glm')
    init <- family$init(y=Y, weights=weights)

    Y <- init$y
    family$weights <- as.double(init$weights)
    family$dispersions <- as.double(family$phi / init$weights) # and zero weights ?!
    family$linPredStart <- as.double(init$linPredStart)
    
  
    ## which terms gave rise to which columns?
    ## (0 = intercept, 1 = first term)
    termNumbers <- attr (X, "assign") 

    ## vector of length col (X) giving uc group indices or 0 (no uc)
    ## for associating uc groups with model matrix columns: ucIndices
    ucIndices <- fpMaxs <- integer (length (termNumbers))

    ## list for mapping group -> columns in model matrix: ucColList
    if(nUcGroups){
        for (i in seq (along=ucTermList)){
            ucIndices[termNumbers %in% ucTermList[[i]]] <-  i
        }

        ucColList <- lapply (seq (along=ucTermList), function (ucGroup) which (ucIndices == ucGroup))
    } else {
        ucColList <- NULL
    }


    ## vectors of length col (X) giving maximum fp degrees or 0 (no bfp)
    ## and if scaling is wanted (1) or not (0)
    ## for associating bfps with model matrix columns
    ## In addition, scale Columns or exit if non-positive values occur
    bfpInds <- bfpTermInd - 1 + attr (Terms, "intercept")           # now indexing matrix column
    for (i in seq (along=bfpInner)){
        colInd <- bfpInds[i]
        fpObj <- bfpInner[[i]]

        fpMaxs[colInd] <- attr (fpObj, "max")

        ## get scaling info
        scaleResult <- fpScale (c(attr(fpObj, "rangeVals"), # extra values not in the data
                                  X[,colInd]),              # covariate data
                                scaling = attr(fpObj, "scale")) # scaling wished?
        attr (bfpInner[[i]], "prescalingValues") <- scaleResult
        
        ## do the scaling
        X[,colInd] <- X[,colInd] + scaleResult$shift
        X[,colInd] <- X[,colInd] / scaleResult$scale

        ## check positivity
        if (min(X[,colInd]) <= 0) 
            stop (simpleError(paste("prescaling necessary for negative values in variable", fpObj)))
    }

    ## check that all maximum FP degrees are equal, so that the SWITCH move
    ## in the model sampling algorithm will always be possible.
    ## This assumption could potentially be removed later on, or otherwise the "max" option in bfp()
    ## could be removed.
    if (length(unique(fpMaxs[fpMaxs != 0])) > 1L)
        stop(simpleError("all maximum FP degrees must be identical"))

    ## check that only the intercept (one column) is a fixed term
    if (sum(! (fpMaxs | ucIndices)) > 1)
        stop (simpleError("only the intercept can be a fixed term"))

    
    ## attach a loglik-function, which is then called from the C++ code.
    ## It gives then the loglikelihood of the mu vector of means.
    family$loglik <- function(mu)
    {
        return(- 0.5 * sum(family$dev.resids(y=Y, mu=mu, wt=weights)) / phi)
    }
    ## note that this does not include normalizing constants of
    ## the sampling density, e.g. - 0.5 * log(2 * pi * phi) is *not*
    ## included in the Gaussian case!

    ## get all information on the null model
    nullModelInfo <- getNullModelInfo(family=family) 
    
    ## compute and print cardinality of the model space to guide decision
    fpSetCards <- ifelse (fpMaxs[fpMaxs != 0] <= 3, 8, 5 + fpMaxs[fpMaxs != 0])
    getNumberPossibleFps <- function (  # computes number of possible univariate fps (including omission)
                                      maxDegree # maximum fp degree
                                      ){
        s <- ifelse (maxDegree <= 3, 8, 5 + maxDegree) # card. of power set
        singleDegreeNumbers <- sapply (0:maxDegree, function (m)
                                       choose (s - 1 + m, m))
        return (sum (singleDegreeNumbers))
    }
    singleNumbers <- sapply (fpMaxs, getNumberPossibleFps)
    totalNumber <- prod (singleNumbers) * 2^(nUcGroups) # maximum number of possible models

    
    ## process the nModels argument
    if(missing(nModels))
    {
        ## then we would like to have the default number of models:
        nModels <- max(1L, floor(totalNumber / 100))
    }
    ## check nModels is at least 1
    if(nModels < 1)
    {
        stop(simpleError("nModels must at least be 1"))
    }
        
        
    ## decide if we are going to do model sampling or an exhaustive search
    if (identical(method, "ask")){
        cat ("The cardinality of the model space is at most ", totalNumber, ".\n", sep = "")
        decision <- substr (readline(paste("Do you want to do a deterministic search for the best model (y)",
                                           "or sample from the model space (n) or abort (else) ?\n")),
                            1, 1)

        ## ensure that correct decision string has been entered, otherwise abort
        if(! decision %in% c("y", "n"))
        {
            cat("Aborting.\n")
            return()
        }

    } else {
        decision <- switch(method, exhaustive = "y", sampling = "n")
    }
    ## and the match.arg before ensures that no further problems can occur.

    ## translate the decision to logical variable 
    doSampling <- identical(decision, "n")
    
    ## ensure that we only do model sampling if there is more than 1 FP term in the model
    if(doSampling && identical(nFps, 1L))
    {
        warning(simpleWarning(paste("We need to do an exhaustive computation of all models,",
                                    "because there is only 1 FP term in the model!")))
        doSampling <- FALSE
    }

    ## if sampling, we possibly ask for the chainlength
    if (doSampling)
    {
        ## get chainlength?
        if (identical(method, "ask"))
        {
            chainlength <- as.numeric (readline ("How long do you want the Markov chain to run?\n"))
        }

        ## compute the default number of models to be saved
        if(is.null(nModels))
        {
            nModels <- as.integer(max(chainlength / 100, 1L))
        }
        else
            stopifnot(nModels >= 1L)

        ## check the chosen cache size
        nCache <- as.integer(nCache)
        stopifnot(nCache >= nModels)        
        
        if (verbose)
        {
            cat("Starting sampler...\n")            
        }
    }
    else                                
    {
        if (verbose)
        {
            cat("Starting with computation of every model...\n")
        }
    }

    ## start the progress bar (is continued in the C++ code)
    if(verbose)
    {
        cat ("0%", rep ("_", 100 - 6), "100%\n", sep = "")
    }

    ## pack the data together
    data <- list(x=X,                   # design matrix
                 xCentered=Xcentered,   # centered design matrix
                 y=as.double(Y),        # response vector
                 nObs=nrow(X))          # number of observations

    ## pack the FP info things together
    fpInfos <- list(fpmaxs=as.integer(fpMaxs[fpMaxs != 0]), # vector of maximum fp degrees)
                    fppos=as.integer(bfpInds),              # vector of fp columns
                    fpcards=as.integer(fpSetCards), # cardinality of corresponding power sets
                    fpnames=if (nFps == 0) character(0) else unlist(bfpInner)) # names of fp terms. Note that
                                        # it is necessary to have this if-else
                                        # construction (ifelse does not work,
                                        # because it would take the first
                                        # element of the length-0 vector
                                        # character(0), which is NA!!)
    
    
    
    ## pack the UC info things together
    ucInfos <- list(ucIndices=as.integer(ucIndices), # vector giving uncertainty custer indices
                                        # (column -> which group) 
                    ucColList=ucColList) # list for group -> which columns mapping

    ## pack model search configuration:
    searchConfig <- list(totalNumber=as.double(totalNumber), # cardinality of model space                         
                         nModels=as.integer(nModels),         # number of best
                                        # models returned
                         empiricalBayes=empiricalBayes, # use EB for g and
                                        # conditional marginal likelihoods?
                         doSampling=doSampling,               # shall model sampling be done? If
                                        # false, then exhaustive search.
                         chainlength=as.double(chainlength),  # how many times should a jump be
                                        # proposed?
                         nCache=nCache, # how many models to cache at the same time
                         largeVariance=as.double(largeVariance), # what is a "large" variance output
                                        # of BFGS?
                         useBfgs=useBfgs) # should we use the BFGS algorithm (or
                                        # Brent's optimize)?                        

    ## pack prior and likelihood information:
    distribution <- list(nullModelInfo=nullModelInfo, # all information about the null model.
                         gPrior=priorSpecs$gPrior, # prior on the covariance
                                        # factor g (S4 class object)
                         modelPrior=priorSpecs$modelPrior, # model prior string                         
                         family=family)    # GLM family and link                                 

    ## pack other options
    options <- list(verbose=verbose,           # should progress be displayed?
                    debug=debug,               # echo debug-style messages?
                    gaussHermite=gaussHermite,   # nodes and weights for Gauss
                                        # Hermite quadratures
                    useOpenMP=useOpenMP, # should we use openMP for speed up?
                    higherOrderCorrection=higherOrderCorrection) # should
                                        # the higher-order Laplace correction be used?    
    
    ## then go C++
    Ret <-
        .External (cpp_glmBayesMfp,
                   data,
                   fpInfos,
                   ucInfos,
                   searchConfig,
                   distribution,
                   options)

    ## C++ attaches the following attributes:

    ## numVisited
    ## inclusionProbs
    ## logNormConst

    ## name the inclusion probabilities
    names (attr (Ret, "inclusionProbs")) <- c(unlist (bfpInner), ucInner)

    ## name the models with the model index
    names (Ret) <- 1:length(Ret)

    ## attach additional information:

    ## information passed to C++, which is important for e.g. the function "sampleGlm"
    attr (Ret, "data") <- data
    attr(Ret, "fpInfos") <- fpInfos
    attr(Ret, "ucInfos") <- ucInfos
    attr(Ret, "searchConfig") <- searchConfig
    attr(Ret, "distribution") <- distribution
    attr(Ret, "options") <- options

    ## original call and formula
    attr (Ret, "call") <- call
    attr (Ret, "formula") <- newFormula

    ## prior specs argument
    attr (Ret, "priorSpecs") <- priorSpecs


    ## list with index info
    fixedInds <- setdiff (1:ncol (X), c (bfpInds, which (ucIndices > 0)))
    attr (Ret, "indices") <- list (uc = ucIndices,
                                   ucList = ucColList,
                                   bfp = bfpInds,
                                   fixed = fixedInds)

    
    ## names of the terms
    fixedNamesInds <- setdiff(2:length (varlist), unlist (attr (Terms, "specials"))) + 1
    interceptName <- ifelse (attr (Terms, "intercept"), "(Intercept)", NULL)
    attr (Ret, "termNames") <- list (fixed=
                                     c(interceptName,
                                       sapply(as.list(vars[fixedNamesInds]),
                                               deparse)),
                                     bfp=unlist(bfpInner),
                                     uc=ucInner)

    ## matrix with shift/scale info, maximum degree and cardinality of powerset
    shiftScaleMaxMat <- matrix (nrow = nFps, ncol = 4)
    colnames (shiftScaleMaxMat) <- c ("shift", "scale", "maxDegree",
                                      "cardPowerset")

    if(nFps > 0L)
    {
        shiftScaleMaxMat[, 1:2] <- matrix(unlist(lapply(bfpInner,
                                                        attr,
                                                        "prescalingValues")),
                                          ncol = 2,
                                          byrow = TRUE)
        shiftScaleMaxMat[, 3] <- fpMaxs[fpMaxs != 0]
        shiftScaleMaxMat[, 4] <- fpSetCards

        rownames (shiftScaleMaxMat) <- unlist (bfpInner)
    }

    attr (Ret, "shiftScaleMax") <- shiftScaleMaxMat

    ## set class and return
    class (Ret) <- c("GlmBayesMfp", "list")
    return(Ret)
}

