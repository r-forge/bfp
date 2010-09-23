#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

#include <set>
#include <map>
#include <vector>
#include "RnewMat.h"
#include <iterator>
#include <types.h>


struct safeSum
{
    LongDoubleVector vals;

    void
    add(const long double &val);

    // compute the sum of the elements using accurate algorithm
    long double
    sum();

    // compute the log of the sum of the exp of the elements using accurate algorithm,
    // and avoiding infinite contributions.
    long double
    logSumExp();

    // compute the sum of the elements using very simple algorithm
    long double
    simpleSum();
};


struct indexSafeSum{
        typedef std::vector<long double>::size_type indexType;
        std::set<indexType > indices;
        void add(const indexType&);
        long double sum(const safeSum&) const;
};



struct book{

    PosLargeInt modelCounter;
    safeSum modelPropToPosteriors;
    std::vector<indexSafeSum> covGroupWisePosteriors; // for computation of covariate inclusion probs: array (bfp, uc)
        bool verbose;
        PosLargeInt chainlength;
        PosLargeInt nanCounter;
        PosInt nModels;
        book() : modelCounter(0), nanCounter(0) {};
};


struct fpInfo{ // collects all information on fractional polynomials needed to be passed down
    PosInt nFps;
    DoubleVector powerset;
    int* fpcards;
    int* fppos;
    int* fpmaxs;
    SEXP fpnames;

    ColumnVectorArray tcols;

    void inds2powers(const Powers &m, double* p) const;

    PosInt maxFpDim;
};


// data structure ###
struct modelInfo{ // contents: must be assignable 
	double logMargLik;
	double logPrior;
	double logPost;
	double postExpectedg; // posterior expected factor g given this model
	double postExpectedShrinkage; // posterior expected shrinkage factor g/(1+g) given this model
	double R2; // coefficient of determination for this model
	
	unsigned long int hits; // only for MCMC, else NA
	
	modelInfo() : // default
		logMargLik(0), logPrior(0), logPost(0), postExpectedg(0), postExpectedShrinkage(0), R2(0), hits(0) {}
	modelInfo(const double &v, const double &w, const double &x, const double &y, const double &z) : // initialize without hits
		logMargLik(v), logPrior(w), logPost(v + w), postExpectedg(x), postExpectedShrinkage(y), R2(z), hits(R_NaReal) {}
	modelInfo(const double &v, const double &w, const double &x, const double &x2, const double &y, const unsigned long int &z) : 
		logMargLik(v), logPrior(w), logPost(v + w), postExpectedg(x), postExpectedShrinkage(x2), R2(y), hits(z) {}	// initialize with hits

	modelInfo& operator=(const modelInfo& m); // assignment operator

	SEXP
	convert2list(double addLogMargLikConst,
                     double subLogPriorConst,
                     long double logNormConst,
                     const book& bookkeep) const;
};
	
struct modelPar{ // key: must have a strict weak ordering
	PowersVector fpPars; // vector of multisets
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

	SEXP
	convert2list(const fpInfo& currFp) const;
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
	
	typedef std::map<modelPar, modelInfo>::size_type NumberType;
	NumberType totalNumber; // cardinality of model space
	
	dataValues(const Matrix &x,
               const Matrix &xcentered,
               const ColumnVector &y,
               const double &totalNum);
};







struct model{
	modelPar par;
	modelInfo info;
	
	model(const modelPar& p, const modelInfo& i) : par(p), info(i) {} // initialize
	model(const model& m) : par(m.par), info(m.info) {}; // copy ctor
	
	model& operator=(const model& m); // assignment operator
	bool operator<(const model& m) const; // less		
	
	SEXP convert2list(const fpInfo& currFp,
	                  double addLogMargLikConst,
	                  double subLogPriorConst,
	                  long double normConst,
	                  const book& bookkeep) const; // model to list
};





// the model cache class.
// Caches the best models in a map of a given maximum size, and also stores the
// (unnormalized) log posterior probabilities in an ordered set, pointing to the models in the map.
class ModelCache {
public:

    // create a new ModelCache with given maximum size.
    ModelCache(int maxSize) :
        maxSize(maxSize),
        modelMap(),
        modelIterSet()
        {
        }

    // check if max size was reached
    bool
    isFull() const
    {
        return modelMap.size() == maxSize;
    }

    // return size of cache
    int
    size() const
    {
        return modelMap.size();
    }

    // insert model parameter and belonging model info into the cache.
    // returns false if not inserted (e.g. because the par was
    // already inside, or the model was not good enough)
    bool
    insert(const modelPar& par, const modelInfo& info);

    // search for the log marginal likelihood of a model config in the map,
    // and return NA if not found
    double
    getLogMargLik(const modelPar& par) const;

    // increment the sampling frequency for a model configuration
    // (of course, if this config is not cached nothing is done!)
    void
    incrementFrequency(const modelPar& par);

    // compute the log normalising constant from all cached models
    long double
    getLogNormConstant() const;

    // compute the inclusion probabilities from all cached models,
    // taking the log normalising constant and the total number of FPs / UC groups
    DoubleVector
    getInclusionProbs(long double logNormConstant, PosInt nFps, PosInt nUcs) const;

    // convert the best nModels from the cache into an R list
    SEXP
    getListOfBestModels(const fpInfo& currFp,
                        double addLogMargLikConst,
                        double subLogPriorConst,
                        long double logNormConst,
                        const book& bookkeep) const;


private:

    // the map type
    typedef std::map<modelPar, modelInfo> MapType;

    // define comparison function for iterators
    struct Compare_map_iterators
    {
        bool
        operator()(const MapType::iterator& first, const MapType::iterator& second) const
        {
            return (first->second.logPost) < (second->second.logPost);
        }
    };

    // the set type of ordered map iterators
    typedef std::set<MapType::iterator, Compare_map_iterators> SetType;

    // and finally the data members
    const MapType::size_type maxSize;
    MapType modelMap;
    SetType modelIterSet;
};

struct modelmcmc{ // all information needed in mcmc function
        modelPar modPar;
        std::set<unsigned int> freeCovs; // indices of free covs (starting from first fp with index 1 up to uc index = nFps + 1)
        std::set<unsigned int> presentCovs; // analogue
        std::set<int> freeUcs; // indices within uc groups, denoting the birthable ones
        unsigned int dim; // number of columns in this model's design matrix
        double birthprob, deathprob, moveprob; // move type probabilites, switchprob is 1-bprob-dprob-mprob.
        double logMargLik;
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
