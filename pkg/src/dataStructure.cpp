#include "dataStructure.h"
#include "RnewMat.h"
#include "sum.h"
#include <set>

#include <algorithm>

using std::lexicographical_compare;

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

SEXP model::convert2list(	const fpInfo& currFp, 
							const double& addLogMargLikConst, 
							const double& subLogPriorConst, 
							const long double& normConst) const // convert model into list for export to R  
{
	// names of components
	const char* components[] = {"powers", "ucTerms", "logM", "logP", "posterior", "postExpectedg", "postExpectedShrinkage", "R2"};
	unsigned int nProtect = 0;
	unsigned int nComponents = sizeof(components) / sizeof(*components);
	
	SEXP ret, powers, ucTerms, logMargLik, logPrior, posterior, postExpectedg, names, R2, postExpectedShrinkage;
	
	 // allocate return list
	Rf_protect(ret = Rf_allocVector(VECSXP, nComponents));
	nProtect++;
	
	// assemble components
	Rf_protect(ucTerms = Rf_allocVector(INTSXP, par.ucSize)); // ucTerms
	nProtect++;
	copy(par.ucPars.begin(), par.ucPars.end(), INTEGER(ucTerms));

	double correctLogMargLik = info.logMargLik + addLogMargLikConst;
	Rf_protect(logMargLik = Rf_ScalarReal(correctLogMargLik)); // logM	
	nProtect++;
	
	double correctLogPrior = info.logPrior - subLogPriorConst;
	Rf_protect(logPrior = Rf_ScalarReal(correctLogPrior)); // logPrior
	nProtect++;
	
	long double propToPost = expl (info.logMargLik + info.logPrior);
	double post = propToPost / normConst;
	Rf_protect(posterior = Rf_ScalarReal(post)); // posterior	
	nProtect++;
	
	Rf_protect(powers = Rf_allocVector(VECSXP, currFp.nFps)); // powers
	nProtect++;
	for (unsigned int i = 0; i != currFp.nFps; i++){ // for
		SEXP thesePowers;
		Rf_protect(thesePowers = Rf_allocVector(REALSXP, par.fpPars[i].size()));
		nProtect++;
		currFp.inds2powers(par.fpPars[i], REAL(thesePowers));
		SET_VECTOR_ELT(powers, i, thesePowers);	
	}	
	Rf_setAttrib(powers, R_NamesSymbol, currFp.fpnames);
	
	Rf_protect(postExpectedg = Rf_ScalarReal(info.postExpectedg)); // postExpectedg	
	nProtect++;

	Rf_protect(postExpectedShrinkage = Rf_ScalarReal(info.postExpectedShrinkage)); // postExpectedShrinkage	
	nProtect++;
	
	Rf_protect(R2 = Rf_ScalarReal(info.R2)); // R^2	
	nProtect++;	
	
	// set names
	Rf_protect(names = Rf_allocVector(STRSXP, nComponents));
	nProtect++;
	for (unsigned int i = 0; i != nComponents; i++){
		SET_STRING_ELT(names, i, Rf_mkChar(components[i]));	
	}
	Rf_setAttrib(ret, R_NamesSymbol, names);
	
	// set values
	SET_VECTOR_ELT(ret, 0, powers);	
	SET_VECTOR_ELT(ret, 1, ucTerms);	
	SET_VECTOR_ELT(ret, 2, logMargLik);	
	SET_VECTOR_ELT(ret, 3, logPrior);	
	SET_VECTOR_ELT(ret, 4, posterior);	
	SET_VECTOR_ELT(ret, 5, postExpectedg);
	SET_VECTOR_ELT(ret, 6, postExpectedShrinkage);
	SET_VECTOR_ELT(ret, 7, R2);
	
	Rf_unprotect(nProtect);	
	return ret;
}

SEXP model::convert2listMcmc(	const fpInfo& currFp, 
								const double& addLogMargLikConst, 
								const double& subLogPriorConst, 
								const long double& normConst, 
								const book& bookkeep) const
{
	SEXP ret, newPosterior;
	unsigned int nProtect = 0;
	
	Rf_protect(ret = convert2list(currFp, addLogMargLikConst, subLogPriorConst, normConst)); // at first the same as for non-mcmc
	nProtect++;
	
	Rf_protect(newPosterior = Rf_allocVector(REALSXP, 2)); // assemble new posterior vector
	nProtect++;
	
	REAL(newPosterior)[0] = REAL(VECTOR_ELT(ret, 4))[0]; // copy old value
	REAL(newPosterior)[1] = info.hits * 1.0 / bookkeep.chainlength; // and attach hit estimate
	SET_VECTOR_ELT(ret, 4, newPosterior);
	
	Rf_unprotect(nProtect);
	return ret;
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




