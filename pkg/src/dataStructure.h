#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

#include <set>
#include <map>
#include <vector>
#include "RnewMat.h"
#include <iterator>

typedef std::multiset<int> powers;

typedef std::vector < powers >  powervecType;

// data structure ###
struct modelInfo{ // contents: must be assignable 
	double logMargLik;
	double logPrior;
	double postExpectedg; // posterior expected factor g given this model
	double postExpectedShrinkage; // posterior expected shrinkage factor g/(1+g) given this model
	double R2; // coefficient of determination for this model
	
	unsigned long int hits; // only for MCMC
	
	modelInfo() : // default
		logMargLik(0), logPrior(0), postExpectedg(0), postExpectedShrinkage(0), R2(0), hits(0) {} 
	modelInfo(const double &v, const double &w, const double &x, const double &y, const double &z) : // initialize without hits
		logMargLik(v), logPrior(w), postExpectedg(x), postExpectedShrinkage(y), R2(z) {}	
	modelInfo(const double &v, const double &w, const double &x, const double &x2, const double &y, const unsigned long int &z) : 
		logMargLik(v), logPrior(w), postExpectedg(x), postExpectedShrinkage(x2), R2(y), hits(z) {}	// initialize with hits
	modelInfo(const modelInfo& m) : // copy constructor
		logMargLik(m.logMargLik), logPrior(m.logPrior), postExpectedg(m.postExpectedg), postExpectedShrinkage(m.postExpectedShrinkage), R2(m.R2), hits(m.hits) {} ; 
	
	modelInfo& operator=(const modelInfo& m); // assignment operator
};
	
struct modelPar{ // key: must have a strict weak ordering
	powervecType fpPars; // vector of multisets	
	unsigned int nFps; // length of vector
	unsigned int fpSize; // number of fp powers
	std::set<int> ucPars; // set of group indices, starting from 1 (!)
	int ucSize; // number of uc Groups included
	
	modelPar() : nFps(0), fpSize(0), ucSize(0) {} // default ctor
	//modelPar(const modelPar& m) : nFps(m.nFps), fpSize(m.fpSize), ucPars(m.ucPars), ucSize(m.ucSize) // copy ctor
	//	{} // default ctor is synthesized
	modelPar(const unsigned int &n, const unsigned int &fp, const int &uc) : nFps(n), fpSize(fp), ucSize(uc) {}
	
	bool operator<(const modelPar& m) const;
	modelPar& operator=(const modelPar& m); // assignment operator
	
	int size() const;
};

struct modelmcmc{ // all information needed in mcmc function
	modelPar modPar;	
	std::set<unsigned int> freeCovs; // indices of free covs (starting from first fp with index 1 up to uc index = nFps + 1)
	std::set<unsigned int> presentCovs; // analogue
	std::set<int> freeUcs; // indices within uc groups, denoting the birthable ones
	unsigned int dim; // number of columns in this model's design matrix
	double birthprob, deathprob, moveprob; // move type probabilites, switchprob is 1-bprob-dprob-mprob.
	double logMargLik;
	map<modelPar, modelInfo>::iterator mapPos;
};
 
struct hyperPriorPars{ 
	double a; // hyperparameter for hyper-g prior on g
	bool useSparsePrior; // use a sparse model prior?
	
	hyperPriorPars(const double &a, const bool &useSparsePrior) :
        a(a),
        useSparsePrior(useSparsePrior)
    {
    }
};

struct dataValues{
	Matrix design;
	Matrix centeredDesign;
	
	ColumnVector response;
	double sumOfSquaresTotal;
	
	int nObs;
	
	ColumnVector onesVector;
	
	map<modelPar, modelInfo>::size_type totalNumber; // cardinality of model space
	
	dataValues(const Matrix &x, const Matrix &xcentered, const ColumnVector &y, const double &totalNum);
};

struct fpInfo{ // collects all information on fractional polynomials needed to be passed down
	unsigned int nFps;
	double* powerset;
	int* fpcards;
	int* fppos;
	int* fpmaxs;
	SEXP fpnames;
	vector<ColumnVector>* tcols;
	void inds2powers(const multiset<int> &m, double* p) const;
	unsigned int maxFpDim;
};

struct safeSum{
	std::vector<long double> vals;
	void add(const long double &val);	
	long double sum();
};

struct indexSafeSum{
	typedef std::vector<long double>::size_type indexType;
	std::set<indexType > indices;
	void add(const indexType&);
	long double sum(const safeSum&) const;
};

struct book{
	unsigned long long int modelCounter;
	safeSum modelPropToPosteriors;
	indexSafeSum *covGroupWisePosteriors; // for computation of covariate inclusion probs: array (bfp, uc)
	bool verbose;
	unsigned long long int chainlength;	
	unsigned long long int nanCounter;
	unsigned int nModels;
	
	book() : modelCounter(0), nanCounter(0) {};
};


struct model{
	modelPar par;
	modelInfo info;
	
	model(const modelPar& p, const modelInfo& i) : par(p), info(i) {} // initialize
	model(const model& m) : par(m.par), info(m.info) {}; // copy ctor
	
	model& operator=(const model& m); // assignment operator
	bool operator<(const model& m) const; // less		
	
	SEXP convert2list(	const fpInfo& currFp, 
						const double& addLogMargLikConst, 
						const double& subLogPriorConst, 
						const long double& normConst) const; // model to list
	SEXP convert2listMcmc(	const fpInfo& currFp, 
							const double& addLogMargLikConst, 
							const double& subLogPriorConst, 
							const long double& normConst, 
							const book&) const; // mcmc model to list, includes hits
};

// delete a number from a set
template <class T>
typename std::set<T> removeElement(std::set<T> input, T element)
{
	typename std::set<T>::iterator iter = input.begin();
    while( iter != input.end() )
    {
      if (*iter == element)
        // A copy of iter is passed into erase(), ++ is executed after erase().
        // Thus iter remains valid
        input.erase( iter++ );
      else
        ++iter;
    }

    return input;
}

// construct a sequence 1:maximum
template <class T>
typename std::set<T> constructSequence(T maximum)
{
	std::set<T> ret;

	for(T i = 1; i <= maximum; ++i)
	{
		ret.insert(ret.end(), i);
	}

	return ret;
}


#endif /*DATASTRUCTURE_H_*/
