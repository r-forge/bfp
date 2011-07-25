#include <zdensity.h>
#include <stdexcept>

#include <rcppExport.h>
#include <linalgInterface.h>

extern "C" {

double F77_NAME(laplace)(const int* nObs, const int* nCoefs, const double* linPred, const double* L);

}


// call the function object
double
NegLogUnnormZDens::operator()(double z)
{
    // map back to the original covariance factor scale
    double const g = exp(z);

    // optional status message
    if(verbose)
    {
        Rprintf("\nNegLogUnnormZDens: Computing function call NegLogUnnormZDens(%f) ...", z);
    }

    // if g is numerically zero, we cannot proceed
    if(g == 0.0)
    {
        std::ostringstream stream;
        stream << "g numerically zero";
        throw std::domain_error(stream.str().c_str());
    }

    // protect against errors in IWLS and non-finite result
    try
    {
        // compute Gaussian approximation:
        PosInt requiredIter = iwlsObject.startWithLastLinPred(nIter, // iterate at most nIter times (because we want convergence)
                                                              g);    // for this covariance factor, for the current linear predictor start.

        // if we did not converge in the maximum number of iterations, then start again from the original linear predictor,
        // and allow even more iterations.
        if(requiredIter > nIter)
        {
            PosInt higherIter = 2 * nIter;
            requiredIter = iwlsObject.startWithNewLinPred(higherIter,
                                                          g,
                                                          linPredStart);

            // if we still do not have converged, print a warning message if we are verbose.
            if(requiredIter > higherIter && verbose)
            {
                Rprintf("\nnNegLogUnnormZDens: IWLS did not converge in %d iterations (even after restart)!", higherIter);
            }

            // do not warn the user, because this need not be a serious problem.
        }

        if(verbose)
        {
            Rprintf("\nNegLogUnnormZDens: finished after %d IWLS iterations", requiredIter);
        }

        // get iwls results
        const IwlsResults iwlsResults = iwlsObject.getResults();

        // then the return value is:
        double ret = 0.5 * iwlsResults.logPrecisionDeterminant -
                     iwlsObject.nCoefs * M_LN_SQRT_2PI -    // This is dependent on the model dimension! (we had to add it here because
                                                            // computeLogUnPosteriorDens also contains now the analogous normalization constant.
                     iwlsObject.computeLogUnPosteriorDens(Parameter(iwlsResults.coefs, z));


        // ------------- begin test

        if(binaryLogisticCorrection)
        {

            // the Raudenbush & Yang & Yosef correction to the Laplace approximation:
            // (now only specific for binary logistic regression !!)

            // compute the matrix L with columns l_i = (iwlsResults.qFactor)^-1 * x_i :

            AMatrix L = trans(iwlsObject.design);

            trs(false,
                false,
                iwlsResults.qFactor,
                L);

            int nObs = iwlsObject.nObs;
            int nCoefs = iwlsObject.nCoefs;

            const double correctionFactor =
                    F77_CALL(laplace)(& nObs,
                                      & nCoefs,
                                      iwlsResults.linPred.memptr(),
                                      L.memptr());

            if(verbose)
            {
                Rprintf("\nNegLogUnnormZDens: Higher-order correction factor is %f", 1.0 + correctionFactor);
            }

            // since we return here the negative log of the conditional marg lik:
            ret -= log1p(correctionFactor);


        } // end if(binaryLogisticCorrection)

        // ------------- end test

        // check finiteness of return value
        if(! R_finite(ret))
        {
            std::ostringstream stream;
            stream << "NegLogUnnormZDens() got non-finite result " << ret << " for z=" << z;
            throw std::domain_error(stream.str().c_str());
        }

        // and only if it is finite we return
        return ret;

    }
    // catch errors in IWLS and non-finite result
    catch (std::domain_error error)
    {
        if(verbose)
        {
            Rprintf("For z=%f, the density value could not be computed because\n%s",
                    z, error.what());
        }
        Rf_warning("for z=%f, the density value could not be computed for the following model:\n%s\nCheck for near-collinearity of covariates.",
                   z, mod.print(fpInfo).c_str());

        // return NaN. This can be handled by the Brent optimize routine! It apparently replaces it (implicitely) by
        // the highest positive number.
        return R_NaN;
    }
}
