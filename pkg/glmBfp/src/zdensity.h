/*
 * zdensity.h
 *
 *  Created on: 16.11.2009
 *      Author: daniel
 */

#ifndef ZDENSITY_H_
#define ZDENSITY_H_

#include <iwls.h>
#include <types.h>

class NegLogUnnormZDens {
public:

    // call the function object:
    // z is the argument,
    double
    operator()(double z);

    // constructor
    NegLogUnnormZDens(const ModelPar &mod,
                      const DataValues& data,
                      const FpInfo& fpInfo,
                      const UcInfo& ucInfo,
                      const GlmModelConfig& config,
                      // return the approximate *conditional* density f(y | z, mod) by operator()?
                      // otherwise return the approximate unnormalized *joint* density f(y, z | mod).
                      bool conditional,
                      bool verbose,
                      bool higherOrderCorrection,
                      PosInt nIter=40) :
                          mod(mod),
                          fpInfo(fpInfo),
                          config(config),
                          linPredStart(config.linPredStart),
    iwlsObject(mod, data, fpInfo, ucInfo, config,
               linPredStart,
               // take the same original start value for each model, but then update
               // it inside the iwls object when new calls to the functor are made.
               conditional, // fixed z if conditional density of y given z is wished
               EPS, // take EPS as the convergence epsilon
               verbose), // and debug if verbose is wished.
               verbose(verbose),
               higherOrderCorrection(higherOrderCorrection),
               nIter(nIter)
               {
               }

private:
    // save the model reference and fp info, so we can write a nice warning message if the
    // IWLS fails for some z
    const ModelPar& mod;
    const FpInfo& fpInfo;
    const GlmModelConfig& config;

    // also save the original start linear predictor
    const AVector linPredStart;

    // the IWLS object
    Iwls iwlsObject;

    // be verbose?
    const bool verbose;

    const bool higherOrderCorrection;

    // number of IWLS iterations
    PosInt nIter;
};


#endif /* ZDENSITY_H_ */
