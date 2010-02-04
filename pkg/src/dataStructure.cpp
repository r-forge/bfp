#include "dataStructure.h"
#include "RnewMat.h"
#include "sum.h"
#include <set>
#include <types.h>

#include <algorithm>

using std::lexicographical_compare;
using std::pair;

// model info functions
modelInfo& modelInfo::operator=(const modelInfo& m)
{ // assignment operator 
	if (this != &m){
		logMargLik = m.logMargLik;
		logPrior = m.logPrior;
		postExpectedg = m.postExpectedg;
		postExpectedShrinkage = m.postExpectedShrinkage;		
		hits = m.hits;
		R2 = m.R2;
	}
	return *this;	
}
	
SEXP
modelInfo::convert2list(double addLogMargLikConst,
                        double subLogPriorConst,
                        long double logNormConst,
                        const book& bookkeep) const
{
    // allocate return list
    SEXP ret;
    Rf_protect(ret = Rf_allocVector(VECSXP, 6));

    SEXP names;
    Rf_protect(names = Rf_allocVector(STRSXP, 6));

    Rf_setAttrib(ret, R_NamesSymbol, names);
    Rf_unprotect(1);

    // assemble components
    SET_VECTOR_ELT(ret, 0, Rf_ScalarReal(logMargLik + addLogMargLikConst));
    SET_STRING_ELT(names, 0, Rf_mkChar("logM"));

    SET_VECTOR_ELT(ret, 1, Rf_ScalarReal(logPrior - subLogPriorConst));
    SET_STRING_ELT(names, 1, Rf_mkChar("logP"));

    SEXP posterior;
    Rf_protect(posterior = Rf_allocVector(REALSXP, 2));

    SET_VECTOR_ELT(ret, 2, posterior);
    SET_STRING_ELT(names, 2, Rf_mkChar("posterior"));

    REAL(posterior)[0] = exp(logPost - logNormConst);
    REAL(posterior)[1] = hits * 1.0 / bookkeep.chainlength;

    Rf_unprotect(1);

    SET_VECTOR_ELT(ret, 3, Rf_ScalarReal(postExpectedg));
    SET_STRING_ELT(names, 3, Rf_mkChar("postExpectedg"));

    SET_VECTOR_ELT(ret, 4, Rf_ScalarReal(postExpectedShrinkage));
    SET_STRING_ELT(names, 4, Rf_mkChar("postExpectedShrinkage"));

    SET_VECTOR_ELT(ret, 5, Rf_ScalarReal(R2));
    SET_STRING_ELT(names, 5, Rf_mkChar("R2"));

    // unprotect and return
    Rf_unprotect(1);
    return ret;
}

// model parameter functions
	
bool modelPar::operator<(const modelPar& m) const
{ 
	// return size() < m.size(); way too easy... lexicographical comparison, starting with uc indices:
	if (ucPars < m.ucPars)
		return true;
	else if (ucPars > m.ucPars)
		return false;
	else // uc indices are equal
		return lexicographical_compare(fpPars.begin(), fpPars.end(), m.fpPars.begin(), m.fpPars.end());	
}

modelPar& modelPar::operator=(const modelPar& m)
{
	if (this != &m){
		fpPars = m.fpPars;
		ucPars = m.ucPars;
		fpSize = m.fpSize;
		ucSize = m.ucSize;
		nFps = m.nFps;
	}
	return *this;	
}

int modelPar::size() const
{
	return fpSize + ucSize;
}

SEXP
modelPar::convert2list(const fpInfo& currFp) const
{
    // allocate return list
    SEXP ret;
    Rf_protect(ret = Rf_allocVector(VECSXP, 2));

    SEXP names;
    Rf_protect(names = Rf_allocVector(STRSXP, 2));

    Rf_setAttrib(ret, R_NamesSymbol, names);
    Rf_unprotect(1);

    // powers
    SEXP powers;
    Rf_protect(powers = Rf_allocVector(VECSXP, nFps));
    Rf_setAttrib(powers, R_NamesSymbol, currFp.fpnames);

    for (PosInt i = 0; i != nFps; i++)
    {
        SEXP thesePowers;
        Rf_protect(thesePowers = Rf_allocVector(REALSXP, fpPars[i].size()));

        currFp.inds2powers(fpPars[i], REAL(thesePowers));

        SET_VECTOR_ELT(powers, i, thesePowers);
        Rf_unprotect(1);
    }

    SET_VECTOR_ELT(ret, 0, powers);
    Rf_unprotect(1);
    SET_STRING_ELT(names, 0, Rf_mkChar("powers"));

    // ucTerms
    SEXP ucTerms;
    Rf_protect(ucTerms = Rf_allocVector(INTSXP, ucSize));

    copy(ucPars.begin(), ucPars.end(), INTEGER(ucTerms));

    SET_VECTOR_ELT(ret, 1, ucTerms);
    Rf_unprotect(1);
    SET_STRING_ELT(names, 1, Rf_mkChar("ucTerms"));

    Rf_unprotect(1);
    return ret;
}


// internal helper function to combine the contents of
// two named R lists into one named R list
static SEXP
combineLists(SEXP firstList, SEXP secondList)
{
    // allocate result list with names
    const R_len_t nFirst = Rf_length(firstList);
    const R_len_t nSecond = Rf_length(secondList);

    SEXP ret;
    Rf_protect(ret = Rf_allocVector(VECSXP, nFirst + nSecond));

    SEXP names;
    Rf_protect(names = Rf_allocVector(STRSXP, nFirst + nSecond));

    Rf_setAttrib(ret, R_NamesSymbol, names);
    Rf_unprotect(1);

    // now fill in contents of the first list
    SEXP firstNames = Rf_getAttrib(firstList, R_NamesSymbol);
    for(R_len_t i = 0; i < nFirst; ++i)
    {
        SET_VECTOR_ELT(ret, i, VECTOR_ELT(firstList, i));
        SET_STRING_ELT(names, i, STRING_ELT(firstNames, i));
    }

    // now fill in contents of the second list
    SEXP secondNames = Rf_getAttrib(secondList, R_NamesSymbol);
    for(R_len_t i = 0; i < nSecond; ++i)
    {
        SET_VECTOR_ELT(ret, i + nFirst, VECTOR_ELT(secondList, i));
        SET_STRING_ELT(names, i + nFirst, STRING_ELT(secondNames, i));
    }

    // unprotect and return
    Rf_unprotect(1);
    return ret;
}


// for combination of modelPar and modelInfo: model //

model& model::operator=(const model& m) // assignment operator
{ 
	if (this != &m){
		info = m.info;
		par = m.par;		
	}
	return *this;
}

bool model::operator<(const model& m) const  // less		
{
	double thisLogPost = info.logMargLik + info.logPrior;
	double mLogPost = m.info.logMargLik + m.info.logPrior;
	if (thisLogPost < mLogPost)
		return true;
	else if (thisLogPost > mLogPost)
		return false;
	else  // posteriors are equal, then the parameter makes the decision
		return m.par < par;
}

SEXP model::convert2list(const fpInfo& currFp,
                         double addLogMargLikConst,
                         double subLogPriorConst,
                         long double logNormConst,
                         const book& bookkeep) const // convert model into list for export to R
{
    return combineLists(par.convert2list(currFp),
                        info.convert2list(addLogMargLikConst,
                                          subLogPriorConst,
                                          logNormConst,
                                          bookkeep));
}

// dataValues //

dataValues::dataValues(const Matrix &x, const Matrix &xcentered, const ColumnVector &y, const double &totalNum) : 
	design(x), centeredDesign(xcentered), response(y), totalNumber(static_cast<map<modelPar, modelInfo>::size_type>(totalNum)) 
{
	// number of observations
	nObs = design.Nrows();

	// nObs long vector of ones
	onesVector = ColumnVector(nObs);
	onesVector = 1.0;
	
	// and the SST
	ColumnVector centeredResponse = response - (response.sum() / nObs) * onesVector;
	sumOfSquaresTotal = centeredResponse.sum_square();   	
}	


// fpInfo //

void fpInfo::inds2powers(const multiset<int> &m, double* p) const // convert inds m into power array p
{
	unsigned int i = 0;
	for (multiset<int>::const_iterator j = m.begin(); j != m.end(); i++, j++){
		p[i] = powerset[*j];
	}
}

// safeSum //

void safeSum::add(const long double &val) 
{
	vals.push_back(val);
}

long double safeSum::sum()
{
	long double ret	= modified_deflation(vals);
	return ret;
}

long double safeSum::logSumExp()
{
    // the maximum of the log contributions is:
    long double maxLogContrib = *std::max_element(vals.begin(), vals.end());

    // now compute the constant which is added to all log contributions,
    // in order to avoid infinite contributions and at the same time use
    // the whole number space (i.e. possibly avoid zero contributions)
    long double constant = log(LDBL_MAX) - 100.0L - maxLogContrib;
    // 100 is for safety.

    // so now the contributions, offset by the constant
    LongDoubleVector expVals;
    for(LongDoubleVector::const_iterator
            l = vals.begin();
            l != vals.end();
            ++l)
    {
        expVals.push_back(exp(*l + constant));
    }

    // the result is the log of the sum, corrected with the constant:
    long double ret = log(modified_deflation(expVals)) - constant;
    return ret;
}

long double safeSum::simpleSum()
{
        long double ret = 0.0;
        for(LongDoubleVector::const_iterator
                v = vals.begin();
                v != vals.end();
                ++v)
        {
            ret += *v;
        }
        return ret;
}


// indexSafeSum //

void indexSafeSum::add(const std::vector<long double>::size_type& ind)
{
	indices.insert(ind);	
}

long double indexSafeSum::sum(const safeSum& s) const
{
	vector<long double> tempVec;
	for(set<indexType>::const_iterator i = indices.begin(); i != indices.end(); i++){
		tempVec.push_back(s.vals.at(*i));	
	}
	return modified_deflation(tempVec);
}


// ModelCache //

// insert model parameter and corresponding info into cache,
// with caring about the maximum number of elements in the map.
bool
ModelCache::insert(const modelPar& par, const modelInfo& info)
{
    // first check size of cache
    if(isFull())
    {
        // if we are full, then check if this log posterior is better than
        // the worst cached model, which is pointed to by
        MapType::iterator worstModelIter = *(modelIterSet.begin());

        // the comparison
        if((worstModelIter->second.logPost) < info.logPost)
        {
            // new model is better than worst model cached.
            // so we delete the worst model from the cache.

            // first from the map
            modelMap.erase(worstModelIter);
            // and then from the set
            modelIterSet.erase(modelIterSet.begin());
        }
        else
        {
            // the new model is not better than the worst model cached,
            // so we do not cache it.
            return false;
        }
    }

    // so now we know that we want to insert the model into the cache,
    // either because the cache was not full or because the new model was better
    // than the worst model cached.

    // -> try inserting into the map:
    pair<MapType::iterator, bool> ret = modelMap.insert(MapType::value_type(par, info));

    // if we were successful:
    if(ret.second)
    {
        // then also insert the iterator pointing to the map element into the set.
        modelIterSet.insert(ret.first);

        // return success
        return true;
    }
    else
    {
        return false;
        Rf_error("Should not happen: model already contained in model cache!");
    }
}

// search for the log marginal likelihood of a model config in the map,
// and return NA if not found
double
ModelCache::getLogMargLik(const modelPar& par) const
{
    // search for the config in the map
    MapType::const_iterator ret = modelMap.find(par);

    // if found, return the log marg lik
    if(ret != modelMap.end())
        return ret->second.logMargLik;
    else
        return R_NaReal;
}

// increment the sampling frequency for a model configuration
// (of course, if this config is not cached nothing is done)
void
ModelCache::incrementFrequency(const modelPar& par)
{
    // search for the config in the map
    MapType::iterator ret = modelMap.find(par);

    // if found, increment the hits
    if(ret != modelMap.end())
        ret->second.hits++;
}

// compute the log normalising constant from all cached models
long double
ModelCache::getLogNormConstant() const
{
    // use safe summation
    safeSum vec;

    // traverse the cache
    for(MapType::const_iterator
            m = modelMap.begin();
            m != modelMap.end();
            ++m)
    {
        // and add all unnormalized log posteriors
        vec.add(m->second.logPost);
    }

    // return the log of the sum of the exp'ed saved elements
    return vec.logSumExp();
}

// compute the inclusion probabilities from all cached models,
// taking the log normalising constant, the number of FPs and the number of UC groups
DoubleVector
ModelCache::getInclusionProbs(long double logNormConstant, PosInt nFps, PosInt nUcs) const
{
    // abbreviation
    typedef std::vector<safeSum> SafeSumVector;
    // allocate vector of safeSum objects for all FPs
    SafeSumVector fps(nFps);

    // and all UC groups
    SafeSumVector ucs(nUcs);

    // now process each model in the cache
    for(MapType::const_iterator
            m = modelMap.begin();
            m != modelMap.end();
            ++m)
    {
        // abbrevs
        const modelPar& thisPar = m->first;
        const modelInfo& thisInfo = m->second;

        // first process the FPs
        {
        SafeSumVector::iterator s = fps.begin();
        for (powervecType::const_iterator
                p = thisPar.fpPars.begin();
                p != thisPar.fpPars.end();
                ++p, ++s)
        {
            // is this FP in the model m?
            if (! p->empty())
            {
                // then add the normalized model probability onto his FP stack
                s->add(exp(thisInfo.logPost - logNormConstant));
            }
        }
        }

        // then process the UC groups
        {
        SafeSumVector::iterator s = ucs.begin();
        for (PosInt i = 1; i <= nUcs; ++i, ++s)
        {
            // is this UC group in the model m?
            if (thisPar.ucPars.find(i) != thisPar.ucPars.end())
            {
                // then add the normalized model probability onto his UC stack
                s->add(exp(thisInfo.logPost - logNormConstant));
            }
        }
        }
    } // end processing all models in the cache

    // so now we can sum up safesum-wise to the return double vector
    DoubleVector ret;

    for(SafeSumVector::iterator
            s = fps.begin();
            s != fps.end();
            ++s)
    {
        ret.push_back(s->sum());
    }

    for(SafeSumVector::iterator
            s = ucs.begin();
            s != ucs.end();
            ++s)
    {
        ret.push_back(s->sum());
    }

    return ret;
}

// convert the best nModels from the cache into an R list
SEXP
ModelCache::getListOfBestModels(const fpInfo& currFp,
                                double addLogMargLikConst,
                                double subLogPriorConst,
                                long double logNormConst,
                                const book& bookkeep) const
{
    // allocate the return R-list
    SEXP ret;

    // cast is necessary for gcc-4.2 on Mac on R-forge.
    Rf_protect(ret = Rf_allocVector(VECSXP, min(bookkeep.nModels, static_cast<PosInt>(modelIterSet.size()))));

    // process the ordered list of best models from the end (because the set is ordered increasingly)
    PosInt i = 0;
    for(SetType::const_reverse_iterator
            s = modelIterSet.rbegin();
            (i < bookkeep.nModels) && (s != modelIterSet.rend());  // so the return list has min(nModels, modelIterSet.size()) elements.
            ++s, ++i)
    {
        // and for this model, combine the config and info lists to one list and
        // put that in the i-th slot of the return list.
        SET_VECTOR_ELT(ret, i, combineLists((**s).first.convert2list(currFp),
                                            (**s).second.convert2list(addLogMargLikConst,
                                                                      subLogPriorConst,
                                                                      logNormConst,
                                                                      bookkeep)));
    }

    // unprotect and return
    Rf_unprotect(1);
    return ret;
}
