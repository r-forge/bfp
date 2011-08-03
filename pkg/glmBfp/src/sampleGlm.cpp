/*
 * sampleGlm.cpp
 *
 *  Created on: 10.12.2009
 *      Author: daniel
 */

#include <rcppExport.h>
#include <combinatorics.h>
#include <dataStructure.h>
#include <types.h>
#include <iwls.h>
#include <bfgs.h>
#include <optimize.h>
#include <fpUcHandling.h>
#include <linalgInterface.h>
#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ***************************************************************************************************//

struct MarginalZ
{
    MarginalZ(const RFunction& logDens,
              const RFunction& gen) :
                  logDens(logDens),
                  gen(gen)
                  {
                  }


    const RFunction logDens;
    const RFunction gen;
};

// ***************************************************************************************************//

struct Options
{
    Options(bool estimateMargLik,
            bool verbose,
            bool debug,
            bool isNullModel,
            bool useFixedZ,
            PosInt iterations,
            PosInt burnin,
            PosInt step) :
        estimateMargLik(estimateMargLik), verbose(verbose),
                debug(debug), isNullModel(isNullModel), useFixedZ(useFixedZ),
                nSamples(ceil((iterations - burnin) * 1.0 / step)),
                iterations(iterations), burnin(burnin),
                step(step)
    {
    }

    const bool estimateMargLik;
    const bool verbose;
    const bool debug;
    const bool isNullModel;
    const bool useFixedZ;

    const PosInt nSamples;
    const PosInt iterations;
    const PosInt burnin;
    const PosInt step;
};

// ***************************************************************************************************//


class Mcmc
{

public:
    // ctr
    Mcmc(const MarginalZ& marginalz, PosInt nObs, PosInt nCoefs) :
        sample(nCoefs),
        proposalInfo(nObs, nCoefs),
        marginalz(marginalz)
        {
        }

    // the current parameter sample
    Parameter sample;

    // the unnormalized log posterior of this sample
    double logUnPosterior;

    // info about the normal proposal distribution given the sampled z
    IwlsResults proposalInfo;

    // compute the log of the normalized proposal density when the z log density is provided
    // normalize correctly in order not to get problems (perhaps) with the Chib-Jeliazkov estimate
    // computation.
    double
    computeLogProposalDens() const
    {
        // Be careful: qFactor is in fact lower-triangular, so a simple multiplication would fail!
        // use instead directly a BLAS routine for this multiplication.
        AVector tmp = sample.coefs - proposalInfo.coefs;
        trmv(false, true, proposalInfo.qFactor, tmp);

        return 0.5 * (proposalInfo.logPrecisionDeterminant - arma::dot(tmp, tmp)) -
               M_LN_SQRT_2PI * proposalInfo.qFactor.n_rows +
               marginalz.logDens(sample.z);
    }

    // non-default assignment operator
    Mcmc&
    operator=(const Mcmc& rhs)
    {
        if(this == &rhs)
        {
            return *this;
        }
        else
        {
            sample = rhs.sample;
            logUnPosterior = rhs.logUnPosterior;
            proposalInfo = rhs.proposalInfo;

            return *this;
        }
    }


private:
    // the marginal z info: same for all Mcmc objects,
    // therefore it is not assigned by the assignment operator
    const MarginalZ marginalz;
    // important: copy the object, because otherwise (reference/pointer)
    // we are not sure that the functions are still available if we do not use "new"
};

// ***************************************************************************************************//


struct Samples
{
    // constructor: basically allocates beta matrix
    Samples(PosInt nCoefs, PosInt nSamples) :
        coefsSamples(nCoefs, nSamples),
        nSaved(0)
        {
        }

    // save a sample consisting of coefs and z
    void
    storeParameters(const Parameter& sample)
    {
        coefsSamples.col(nSaved++) = sample.coefs;
        zSamples.push_back(sample.z);
    }

    // save terms for marginal likelihood estimate
    void
    storeMargLikTerms(double num, double denom)
    {
        numerator.push_back(num);
        denominator.push_back(denom);
    }

    // output everything to an R list
    List
    convert2list() const;

private:
    // nCoefs x nSamples:
    AMatrix coefsSamples;

    // counts the number of saved parameters, so that we know where to store the
    // next coefficients vector
    PosInt nSaved;

    // is gradually extended:
    DoubleVector zSamples;

    // possibly stays empty if not required by the user:
    // the numerator and denominator terms for the marginal likelihood estimate
    DoubleVector numerator;
    DoubleVector denominator;
};

List
Samples::convert2list() const
{
    return List::create(_["coefficients"] = coefsSamples,
                        _["z"] = zSamples,
                        _["margLikNumerator"] = numerator,
                        _["margLikDenominator"] = denominator);
}


// ***************************************************************************************************//


// get a vector with normal variates from N(mean, sd^2)
AVector
drawNormalVariates(PosInt n, double mean, double sd)
{
    AVector ret(n);

    // use R's random number generator
    GetRNGstate();

    for (PosInt i = 0; i < n; ++i)
    {
        ret(i) = Rf_rnorm(mean, sd);
    }

    // no RNs required anymore
    PutRNGstate();

    return ret;
}

// draw a single random normal vector from N(mean, (precisionCholeskyFactor * t(precisionCholeskyFactor))^(-1))
AVector
drawNormalVector(const AVector& mean,
                 const AMatrix& precisionCholeskyFactor)
{
    // get vector from N(0, I)
    AVector w = drawNormalVariates(mean.n_rows, // as many normal variates as required by the dimension.
                                   0,
                                   1);

    // then solve L' * ret = w, and overwrite w with the result:
    trs(false,
        true,
        precisionCholeskyFactor,
        w);

    // return the shifted vector
    return (w + mean);
}


// draw a single uniform random variable:
// be careful with the seed because the z generator function also uses it (via R)
double
unif()
{
    GetRNGstate();

    double ret = unif_rand();

    PutRNGstate();

    return ret;
}


// ***************************************************************************************************//


// R call is:
//
//    samples <- .External(cpp_sampleGlm,
//                         model,
//                         attrs$data,
//                         attrs$fpInfos,
//                         attrs$ucInfos,
//                         attrs$distribution,
//                         newdata,
//                         options,
//                         marginalz)

SEXP
cpp_sampleGlm(SEXP r_interface)
{
    // ----------------------------------------------------------------------------------
    // extract arguments
    // ----------------------------------------------------------------------------------

    r_interface = CDR(r_interface);
    List rcpp_model(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_data(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_fpInfos(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_ucInfos(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_distribution(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_options(CAR(r_interface));

    r_interface = CDR(r_interface);
    List rcpp_marginalz(CAR(r_interface));

    // ----------------------------------------------------------------------------------
    // unpack the R objects
    // ----------------------------------------------------------------------------------

    // data:
    const NumericMatrix n_x = rcpp_data["x"];
    const AMatrix x(n_x.begin(), n_x.nrow(),
                   n_x.ncol(), false);

    const NumericMatrix n_xCentered = rcpp_data["xCentered"];
    const AMatrix xCentered(n_xCentered.begin(), n_xCentered.nrow(),
                           n_xCentered.ncol(), false);

    const NumericVector n_y = rcpp_data["y"];
    const AVector y(n_y.begin(), n_y.size(),
                   false);

    // FP configuration:

    // vector of maximum fp degrees
    const PosIntVector fpmaxs = as<PosIntVector>(rcpp_fpInfos["fpmaxs"]);
    // corresponding vector of fp column indices
    const PosIntVector fppos = rcpp_fpInfos["fppos"];
    // corresponding vector of power set cardinalities
    const PosIntVector fpcards = rcpp_fpInfos["fpcards"];
    // names of fp terms
    const StrVector fpnames = rcpp_fpInfos["fpnames"];


    // UC configuration:

    const PosIntVector ucIndices = rcpp_ucInfos["ucIndices"];
    List rcpp_ucColList = rcpp_ucInfos["ucColList"];

    std::vector<PosIntVector> ucColList;
    for (R_len_t i = 0; i != rcpp_ucColList.length(); ++i)
    {
        ucColList.push_back(as<PosIntVector>(rcpp_ucColList[i]));
    }


    // distributions info:

    List rcpp_nullModelInfo = rcpp_distribution["nullModelInfo"];
    S4 rcpp_gPrior = rcpp_distribution["gPrior"];
    List rcpp_family = rcpp_distribution["family"];


    // options:

    const bool estimateMargLik = as<bool>(rcpp_options["estimateMargLik"]);
    const bool verbose = as<bool>(rcpp_options["verbose"]);
    const bool debug = as<bool>(rcpp_options["debug"]);
    const bool isNullModel = as<bool>(rcpp_options["isNullModel"]);
    const bool useFixedZ = as<bool>(rcpp_options["useFixedZ"]);
#ifdef _OPENMP
    const bool useOpenMP = as<bool>(rcpp_options["useOpenMP"]);
#endif

    S4 rcpp_mcmc = rcpp_options["mcmc"];
    const PosInt iterations = rcpp_mcmc.slot("iterations");
    const PosInt burnin = rcpp_mcmc.slot("burnin");
    const PosInt step = rcpp_mcmc.slot("step");


    // z density stuff:

    const RFunction logMarginalZdens(as<SEXP>(rcpp_marginalz["logDens"]));
    const RFunction marginalZgen(as<SEXP>(rcpp_marginalz["gen"]));


    // ----------------------------------------------------------------------------------
    // further process arguments
    // ----------------------------------------------------------------------------------

    // data:

     // only the intercept is always included, that is fixed, in the model
     IntSet fixedCols;
     fixedCols.insert(1);

     // totalnumber is set to 0 because we do not care about it.
     const DataValues data(x, xCentered, y, 0, fixedCols);

     // FP configuration:
     const FpInfo fpInfo(fpcards, fppos, fpmaxs, fpnames, x);

     // UC configuration:

     // determine sizes of the UC groups, and the total size == maximum size reached together by all
     // UC groups.
     PosIntVector ucSizes;
     PosInt maxUcDim = 0;
     for (std::vector<PosIntVector>::const_iterator cols = ucColList.begin(); cols != ucColList.end(); ++cols)
     {
         PosInt thisSize = cols->size();

         maxUcDim += thisSize;
         ucSizes.push_back(thisSize);
     }
     const UcInfo ucInfo(ucSizes, maxUcDim, ucIndices, ucColList);

     // model configuration:
     GlmModelConfig config(rcpp_family, rcpp_nullModelInfo, rcpp_gPrior,
                           data.response, debug);


     // model config/info:
     const Model thisModel(ModelPar(rcpp_model["configuration"],
                                   fpInfo),
                          GlmModelInfo(as<List>(rcpp_model["information"])));


     // the options
     const Options options(estimateMargLik,
                           verbose,
                           debug,
                           isNullModel,
                           useFixedZ,
                           iterations,
                           burnin,
                           step);

     // marginal z stuff
     const MarginalZ marginalZ(logMarginalZdens,
                               marginalZgen);


     // use only one thread if we do not want to use openMP.
#ifdef _OPENMP
     if(! useOpenMP)
     {
         omp_set_num_threads(1);
     } else {
         omp_set_num_threads(omp_get_num_procs());
     }
#endif


     // ----------------------------------------------------------------------------------
     // prepare the sampling
     // ----------------------------------------------------------------------------------

     // construct IWLS object, which can be used for all IWLS stuff,
     // and also contains the design matrix etc
     Iwls iwlsObject(thisModel.par,
                     data,
                     fpInfo,
                     ucInfo,
                     config,
                     config.linPredStart,
                     options.useFixedZ,
                     EPS,
                     options.debug);

     // check that we have the same answer about the null model as R
     assert(iwlsObject.isNullModel == options.isNullModel);

     // allocate sample container
     Samples samples(iwlsObject.nCoefs, options.nSamples);

     // count how many proposals we have accepted:
     PosInt nAccepted(0);

     // at what z do we start?
     double startZ = useFixedZ ? rcpp_options["fixedZ"] : thisModel.info.zMode;

     // get the mode for beta given the mode of the approximated marginal posterior as z
     PosInt iwlsIterations = iwlsObject.startWithNewLinPred(30,
                                                            // this is the corresponding g
                                                            exp(startZ),
                                                            // and the start value for the linear predictor is taken from the Glm model config
                                                            config.linPredStart);

     // echo debug-level message?
     if(options.debug)
     {
         Rprintf("\ncpp_sampleGlm: Initial IWLS for high density point finished after %d iterations",
                 iwlsIterations);
     }

     // start container with current things
     Mcmc now(marginalZ, data.nObs, iwlsObject.nCoefs);

     // this is the current proposal info:
     now.proposalInfo = iwlsObject.getResults();

     // and this is the current parameters sample:
     now.sample = Parameter(now.proposalInfo.coefs,
                            startZ);

     // compute the (unnormalized) log posterior of the proposal
     now.logUnPosterior = iwlsObject.computeLogUnPosteriorDens(now.sample);

     // so the parameter object "now" is then also the high density point
     // required for the marginal likelihood estimate:
     const Mcmc highDensityPoint(now);

     // we accept this starting value, so initialize "old" with the same ones
     Mcmc old(now);

     // ----------------------------------------------------------------------------------
     // start sampling
     // ----------------------------------------------------------------------------------

     // echo debug-level message?
     if(options.debug)
     {
         Rprintf("\ncpp_sampleGlm: Starting MCMC loop");
     }


     // i_iter starts at 1 !!
     for(PosInt i_iter = 1; i_iter <= options.iterations; ++i_iter)
     {
         // echo debug-level message?
         if(options.debug)
         {
             Rprintf("\ncpp_sampleGlm: Starting iteration no. %d", i_iter);
         }

         // ----------------------------------------------------------------------------------
         // store the proposal
         // ----------------------------------------------------------------------------------

         // sample one new log covariance factor z (other arguments than 1 are not useful
         // with the current setup of the RFunction wrapper class)
         now.sample.z = marginalZ.gen(1);

         // then do 1 IWLS step, starting from the last linear predictor and the new z
         // (here the return value is not very interesting, as it must be 1)
         iwlsObject.startWithNewCoefs(1,
                                      exp(now.sample.z),
                                      now.sample.coefs);

         // get the results
         now.proposalInfo = iwlsObject.getResults();

         // draw the proposal coefs:
         now.sample.coefs = drawNormalVector(now.proposalInfo.coefs,
                                             now.proposalInfo.qFactor);

         // compute the (unnormalized) log posterior of the proposal
         now.logUnPosterior = iwlsObject.computeLogUnPosteriorDens(now.sample);

         // ----------------------------------------------------------------------------------
         // get the reverse jump normal density
         // ----------------------------------------------------------------------------------

         // copy the old Mcmc object
         Mcmc reverse(old);

         // do again 1 IWLS step, starting from the sampled linear predictor and the old z
         iwlsObject.startWithNewCoefs(1,
                                      exp(reverse.sample.z),
                                      now.sample.coefs);

         // get the results for the reverse jump Gaussian:
         // only the proposal has changed in contrast to the old container,
         // the sample stays the same!
         reverse.proposalInfo = iwlsObject.getResults();


         // ----------------------------------------------------------------------------------
         // compute the proposal density ratio
         // ----------------------------------------------------------------------------------

         // first the log of the numerator, i.e. log(f(old | new)):
         double logProposalRatioNumerator = reverse.computeLogProposalDens();

         // second the log of the denominator, i.e. log(f(new | old)):
         double logProposalRatioDenominator = now.computeLogProposalDens();

         // so the log proposal density ratio is
         double logProposalRatio = logProposalRatioNumerator - logProposalRatioDenominator;

         // ----------------------------------------------------------------------------------
         // compute the posterior density ratio
         // ----------------------------------------------------------------------------------

         double logPosteriorRatio = now.logUnPosterior - old.logUnPosterior;

         // ----------------------------------------------------------------------------------
         // accept or reject proposal
         // ----------------------------------------------------------------------------------

         double acceptanceProb = exp(logPosteriorRatio + logProposalRatio);

         if(unif() < acceptanceProb)
         {
             old = now;

             ++nAccepted;
         }
         else
         {
             now = old;
         }

         // ----------------------------------------------------------------------------------
         // store the sample?
         // ----------------------------------------------------------------------------------

         // if the burnin was passed and we are at a multiple of step beyond that, then store
         // the sample.
         if((i_iter > options.burnin) &&
            (((i_iter - options.burnin) % options.step) == 0))
         {
             // echo debug-level message
             if(options.debug)
             {
                 Rprintf("\ncpp_sampleGlm: Storing samples of iteration no. %d", i_iter);
             }

             // store the current parameter sample
             samples.storeParameters(now.sample);

             // ----------------------------------------------------------------------------------
             // compute marginal likelihood terms
             // ----------------------------------------------------------------------------------

             // compute marginal likelihood terms and save them?
             if(options.estimateMargLik)
             {
                 // echo debug-level message?
                 if(options.debug)
                 {
                     Rprintf("\ncpp_sampleGlm: Compute marginal likelihood estimation terms");
                 }

                 // ----------------------------------------------------------------------------------
                 // compute next term for the denominator
                 // ----------------------------------------------------------------------------------

                 // draw from the high density point proposal distribution
                 Mcmc denominator(highDensityPoint);
                 denominator.sample.z = marginalZ.gen(1);

                 iwlsObject.startWithNewLinPred(1,
                                                exp(denominator.sample.z),
                                                highDensityPoint.proposalInfo.linPred);

                 denominator.proposalInfo = iwlsObject.getResults();

                 denominator.sample.coefs = drawNormalVector(denominator.proposalInfo.coefs,
                                                             denominator.proposalInfo.qFactor);

                 // get posterior density of the sample
                 denominator.logUnPosterior = iwlsObject.computeLogUnPosteriorDens(denominator.sample);

                 // get the proposal density at the sample
                 double denominator_logProposalDensity = denominator.computeLogProposalDens();

                 // then the reverse stuff:
                 // first we copy again the high density point
                 Mcmc revDenom(highDensityPoint);

                 // but choose the new sampled coefficients as starting point
                 iwlsObject.startWithNewCoefs(1,
                                              exp(revDenom.sample.z),
                                              denominator.sample.coefs);
                 revDenom.proposalInfo = iwlsObject.getResults();

                 // so the reverse proposal density is
                 double revDenom_logProposalDensity = revDenom.computeLogProposalDens();


                 // so altogether the next term for the denominator is the following acceptance probability
                 double denominatorTerm = denominator.logUnPosterior - highDensityPoint.logUnPosterior +
                                          revDenom_logProposalDensity - denominator_logProposalDensity;
                 denominatorTerm = exp(fmin(0.0, denominatorTerm));

                 // ----------------------------------------------------------------------------------
                 // compute next term for the numerator
                 // ----------------------------------------------------------------------------------

                 // compute the proposal density of the current sample starting from the high density point
                 Mcmc numerator(now);

                 iwlsObject.startWithNewLinPred(1,
                                                exp(numerator.sample.z),
                                                highDensityPoint.proposalInfo.linPred);
                 numerator.proposalInfo = iwlsObject.getResults();

                 double numerator_logProposalDensity = numerator.computeLogProposalDens();

                 // then compute the reverse proposal density of the high density point when we start from the current
                 // sample
                 Mcmc revNum(highDensityPoint);

                 iwlsObject.startWithNewCoefs(1,
                                              exp(revNum.sample.z),
                                              now.sample.coefs);
                 revNum.proposalInfo = iwlsObject.getResults();

                 double revNum_logProposalDensity = revNum.computeLogProposalDens();

                 // so altogether the next term for the numerator is the following guy:
                 double numeratorTerm = exp(fmin(revNum_logProposalDensity,
                                                 highDensityPoint.logUnPosterior - now.logUnPosterior +
                                                 numerator_logProposalDensity));

                 // ----------------------------------------------------------------------------------
                 // finally store both terms
                 // ----------------------------------------------------------------------------------

                 samples.storeMargLikTerms(numeratorTerm, denominatorTerm);

             }
         }

         // ----------------------------------------------------------------------------------
         // echo progress?
         // ----------------------------------------------------------------------------------

         // echo debug-level message?
         if(options.debug)
         {
             Rprintf("\ncpp_sampleGlm: Finished iteration no. %d", i_iter);
         }

         if((i_iter % std::max(static_cast<int>(options.iterations / 100), 1) == 0) &&
             options.verbose)
         {
             // display computation progress at each percent
             Rprintf("-");

         } // end echo progress

     } // end MCMC loop


     // echo debug-level message?
     if(options.debug)
     {
         Rprintf("\ncpp_sampleGlm: Finished MCMC loop");
     }


     // ----------------------------------------------------------------------------------
     // build up return list for R and return that.
     // ----------------------------------------------------------------------------------

     return List::create(_["samples"] = samples.convert2list(),
                         _["nAccepted"] = nAccepted,
                         _["highDensityPointLogUnPosterior"] = highDensityPoint.logUnPosterior);

} // end cpp_sampleGlm

// ***************************************************************************************************//

// End of sampleGlm.cpp
