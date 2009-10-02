// todo:
// correct normalized value of norm constant must be returned!

#include <RnewMat.h>
#include <combinatorics.h>
#include <dataStructure.h>
#include <hyperg.h>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <iostream>
#include <cassert>
#include <climits>

// using pretty much:
using std::map;
using std::set;
using std::vector;
using std::accumulate;
using std::find;
using std::set_difference;
using std::count;
using std::max_element;

typedef std::vector<long double>::size_type indexType;

// export to C interface ##########################################################################
extern "C"
{
SEXP exhaustiveGaussian(// declaration
                        SEXP R_x, // (not centered!) design matrix (with colnames)
                        SEXP R_xcentered, // centered design matrix
                        SEXP R_y, // response vector
                        SEXP R_fpmaxs, // vector of maximum fp degrees
                        SEXP R_fppos, // corresponding vector of fp column indices
                        SEXP R_fpcards, // corresponding vector of power set cardinalities
                        SEXP R_nFps, // number of fp terms
                        SEXP R_fpnames, // names of fp terms
                        SEXP R_ucIndices, // vector giving _unc_ertainty custer indices (column -> which group)
                        SEXP R_ucTermList, // list for group -> which columns mapping
                        SEXP R_nUcGroups, // number of uncertainty groups
                        SEXP R_totalNumber, // cardinality of model space,
                        SEXP R_hyperparam, // hyperparameter a for hyper-g prior
                        SEXP R_useSparsePrior, // use sparse model prior?
                        SEXP R_nModels, // number of best models to be returned
                        SEXP R_verbose); // should progress been displayed?

SEXP samplingGaussian(// declaration
                      SEXP R_x, // (not centered!) design matrix (with colnames)
                      SEXP R_xcentered, // centered design matrix
                      SEXP R_y, // response vector
                      SEXP R_fpmaxs, // vector of maximum fp degrees
                      SEXP R_fppos, // corresponding vector of fp column indices
                      SEXP R_fpcards, // corresponding vector of power set cardinalities
                      SEXP R_nFps, // number of fp terms
                      SEXP R_fpnames, // names of fp terms
                      SEXP R_ucIndices, // vector giving _unc_ertainty custer indices (column -> which group)
                      SEXP R_ucTermList, // list for group -> which columns mapping
                      SEXP R_nUcGroups, // number of uncertainty groups
                      SEXP R_hyperparam, // hyperparameter a for hyper-g prior
                      SEXP R_useSparsePrior, // use sparse model prior?
                      SEXP R_nModels, // number of best models to be returned
                      SEXP R_verbose, // should progress been displayed?
                      SEXP R_chainlength); // how many times should a jump been made?

SEXP logMargLik( //declaration
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha, // hyperparamater for hyper-g prior
                SEXP R_sst); // total sum of squares computed from y

SEXP postExpectedg( //declaration
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha); // hyperparamater for hyper-g prior

SEXP postExpectedShrinkage( //declaration
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha); // hyperparamater for hyper-g prior

} // extern "C"

// other functions ##########################################################################
void permPars(unsigned int pos, // current position in parameter vector, starting from 0
              const fpInfo &currFp,
              const int &nUcGroups,
              modelPar mod,
              set<model> &space,
              const hyperPriorPars &hyp,
              const dataValues &data,
              set<int>* const &ucTermList,
              const set<int> &fixedCols,
              book&);

set<int> getFreeUcs( // compute set of free uc group indices
                   const modelPar& mod,
                   const vector<unsigned int>& ucSizes,
                   const unsigned int& currDim,
                   const unsigned int& maxDim);

set<unsigned int> getFreeCovs( // compute set of free cov indices
                             const modelPar& mod,
                             const fpInfo& currFp,
                             const set<int>& freeUcs,
                             const unsigned int& currDim,
                             const unsigned int& maxDim);
set<unsigned int> getPresentCovs( // determine set of present cov indices
const modelPar& mod);

template <class T> T discreteUniform( // return random element of myset; should be enclosed in getRNGstate() etc.
const set<T>& myset);

template <class T> typename T::iterator dU( // return iterator of random element of myset; should be enclosed in getRNGstate() etc.
const T& myset);

int discreteUniform( // get random int x with lower <= x < upper; should be enclosed in getRNGstate() etc.
                    const int& lower,
                    const int& upper);

void computeModel(const modelPar &mod,
                  const hyperPriorPars &hyp,
                  const dataValues &data,
                  const fpInfo &currFp,
                  set<int>* const &ucTermList,
                  const int &nUcGroups,
                  const set<int> &fixedCols,
                  set<model> &space,
                  book&);

ReturnMatrix getDesignMatrix( // construct design matrix for the model
                             const modelPar &mod,
                             const dataValues &data,
                             const fpInfo &currFp,
                             set<int>* const &ucTermList,
                             const int &nUcGroups,
                             const set<int> &fixedCols);

double getR2( // compute coefficient of determination for the model
             const Matrix &design,
             const dataValues &data,
             const set<int> &fixedCols,
             const hyperPriorPars &hyp);

double getVarLogMargLik( // compute varying part of log marginal likelihood for specific model
                        const double &R2,
                        const int &n,
                        const int &dim,
                        const hyperPriorPars &hyp);

double getVarLogPrior( // compute logarithm of model prior
                      const modelPar &mod,
                      const fpInfo &currFp,
                      const hyperPriorPars &hyp);

ReturnMatrix getFpMatrix( // build Fp basis matrix from transformed cols and power indices
                         const vector<ColumnVector> &tcols,
                         const multiset<int> &powerinds,
                         const dataValues &data);

void pushInclusionProbs( // push back index into covGroupWisePosteriors-Array
                        const modelPar &mod,
                        const fpInfo &currFp,
                        const int &nUcGroups,
                        book &bookkeep);




// definitions ####################################################################################

SEXP exhaustiveGaussian(// definition
				SEXP R_x, 				// design matrix
				SEXP R_xcentered,       // centered design matrix
                SEXP R_y, 				// response vector
                SEXP R_fpmaxs, 			// vector of maximum fp degrees
                SEXP R_fppos,			// corresponding vector of fp column indices
	            SEXP R_fpcards,			// corresponding vector of power set cardinalities
                SEXP R_nFps, 			// number of fp terms
                SEXP R_fpnames,			// names of fp terms
                SEXP R_ucIndices, 		// vector giving _unc_ertainty custer indices (column -> which group)
                SEXP R_ucTermList,   	// list for (group -> which columns) one-to-many mapping
                SEXP R_nUcGroups, 		// number of uncertainty groups
                SEXP R_totalNumber,	 	// cardinality of model space
                SEXP R_hyperparam,		// hyperparameter a for hyper-g prior
                SEXP R_useSparsePrior, // use sparse model prior?
                SEXP R_nModels,			// number of best models to be returned
                SEXP R_verbose			// should progress been displayed?
				)
{

	unsigned int nProtect = 0;
	// unpack ###
	// data
	const Matrix x = getMatrix(R_x);
	const Matrix xcentered = getMatrix(R_xcentered);
	const ColumnVector y = vec2col(R_y);

	// prior specifications
	const double hyperparam = Rf_asReal(R_hyperparam);
	const bool useSparsePrior = Rf_asLogical(R_useSparsePrior);

	hyperPriorPars hyp(hyperparam,
	                   useSparsePrior);

	const double totalNumber = REAL(R_totalNumber)[0]; // cardinality of model space

	// constant information
	const dataValues data(x, xcentered, y, totalNumber);

	// fp info
	fpInfo currentFpInfo;
	currentFpInfo.nFps = INTEGER(R_nFps)[0];
	currentFpInfo.fpmaxs = INTEGER(R_fpmaxs);
	currentFpInfo.fppos = INTEGER(R_fppos);
	currentFpInfo.fpcards = INTEGER(R_fpcards);
	Rf_protect(currentFpInfo.fpnames = R_fpnames);
	nProtect++;

	const int* biggestMaxDegree = max_element(currentFpInfo.fpmaxs, currentFpInfo.fpmaxs + Rf_length(R_fpmaxs));
	double powerset[max(8, 5 + *biggestMaxDegree)];
	// corresponding indices        0   1     2  3    4  5  6  7
	const double fixedpowers[] = { -2, -1, -0.5, 0, 0.5, 1, 2, 3 }; // always in powerset
	copy(fixedpowers, fixedpowers + 8, powerset);
	for(int more = 3; more < *biggestMaxDegree; more++){ // additional powers
		*(powerset + 8 + (3 - more)) = more + 1;
	}
	currentFpInfo.powerset = powerset; // save maximum powerset

	vector<ColumnVector> transformedCols[currentFpInfo.nFps]; // build array of vectors of ColumnVectors holding the required
									// transformed values for the design matrices
	for (unsigned int i = 0; i != currentFpInfo.nFps; i++){ // for every fp term
		const int nCols = currentFpInfo.fpcards[i];
		const ColumnVector thisCol = x.Column(currentFpInfo.fppos[i]);
		vector<ColumnVector> thisFp;
		for (int j = 0; j != nCols; j++){ // for every possible power
			ColumnVector thisTransform = thisCol;
			double thisPower = currentFpInfo.powerset[j];
			if(thisPower){ // not 0
				for (int k = 0; k != thisTransform.Nrows(); k++){ // transform each element
					assert(thisTransform.element(k) > 0);
					thisTransform.element(k) = pow (thisTransform.element(k), thisPower);
					assert(! ISNAN(thisTransform.element(k)));
				}
			} else { // 0
				for (int k = 0; k != thisTransform.Nrows(); k++){
					assert(thisTransform.element(k) > 0);
					thisTransform.element(k) = log (thisTransform.element(k));
					assert(! ISNAN(thisTransform.element(k)));
				}
			}
			// do not! center the column. This is done inside getFpMatrix, because
			// the repeated powers case cannot be treated here!!

			// and put it into vector of columns
			thisFp.push_back(thisTransform);
		}
		transformedCols[i] = thisFp;
	}
	currentFpInfo.tcols = transformedCols;

	// uc info
	const int* ucIndicesArray = INTEGER(R_ucIndices);
	const vector<int> ucIndices(ucIndicesArray, ucIndicesArray + Rf_length(R_ucIndices));
	const int nUcGroups = INTEGER(R_nUcGroups)[0];
	set<int> ucTermList[nUcGroups]; // Array with length nUcGroups of int-vectors
	if(nUcGroups){ // catch case with no uc groups
		for(R_len_t i = 0; i != Rf_length(R_ucTermList); i++){
			SEXP temp = VECTOR_ELT(R_ucTermList, i);
			copy(INTEGER(temp), INTEGER(temp) + Rf_length(temp), inserter(ucTermList[i], ucTermList[i].begin()));
		}
	}

	set<int> nonfixedCols = set<int>(currentFpInfo.fppos, currentFpInfo.fppos + currentFpInfo.nFps); // fixed info
	set<int> allCols, fixedCols;
	for(vector<int>::size_type i = 0; i != ucIndices.size(); i++){
		allCols.insert(i + 1);
		if (ucIndices.at(i))
			nonfixedCols.insert(i + 1);
	}
	set_difference(allCols.begin(), allCols.end(), nonfixedCols.begin(), nonfixedCols.end(),
				   inserter(fixedCols, fixedCols.begin()));
	// now fixedCols contains indices of columns that are always present in the design matrix


	// no map needed for exhaustive search, a set is the right thing:
	set<model> orderedModels;

	// compute marginal likelihood for every model ###
	if (orderedModels.max_size() < totalNumber)
		Rf_error("\nmodel space is too large - cannot compute every model\n");

	// start model
	modelPar startModel(currentFpInfo.nFps, 0, 0);
	powervecType startFps(currentFpInfo.nFps); // allocate correct length of vector
	startModel.fpPars = startFps;

	// bookkeeping
	book bookkeep;
	bookkeep.verbose = LOGICAL(R_verbose)[0];

	// for computation of inclusion probs
	indexSafeSum cgwp[currentFpInfo.nFps + nUcGroups];
	bookkeep.covGroupWisePosteriors = cgwp;

	// how many models to return?
	if (Rf_length(R_nModels) == 0){
		bookkeep.nModels = max(static_cast<unsigned int>(floor(totalNumber * 1 / 100)), static_cast<unsigned int>(1));
	} else {
		bookkeep.nModels = INTEGER(R_nModels)[0];
	}

	// start computation
	permPars(0, currentFpInfo, nUcGroups, startModel, orderedModels, hyp, data, ucTermList, fixedCols, bookkeep);

	if (bookkeep.verbose){
		Rprintf("\nActual number of possible models:  %d ", bookkeep.modelCounter);
		Rprintf("\nNumber of non-identifiable models: %d", bookkeep.nanCounter);
		Rprintf("\nNumber of saved possible models:   %d\n", orderedModels.size());
	}

	// normalize posterior probabilities and correct log marg lik and log prior of the models to return
	const long double normConst = bookkeep.modelPropToPosteriors.sum();
	const double logMargLikConst = 	lgammafn((data.nObs - 1) / 2.0) -
									(data.nObs - 1) * sqrt(data.sumOfSquaresTotal) -
									(data.nObs - 1) * M_LN_SQRT_PI -
									0.5 * log(data.nObs);
	const double logPriorConst = nUcGroups * M_LN2; // note: M_LN2 = log(2)

	// inclusion probs
	SEXP inc;
	Rf_protect(inc = Rf_allocVector(REALSXP, currentFpInfo.nFps + nUcGroups));
	nProtect++;
	for (int i = 0; i != Rf_length(inc); i++)
		REAL(inc)[i] = bookkeep.covGroupWisePosteriors[i].sum(bookkeep.modelPropToPosteriors) / normConst;

	SEXP ret;
	Rf_protect(ret = Rf_allocVector(VECSXP, orderedModels.size()));
	nProtect++;

	unsigned int i = 0;
	for(set<model>::reverse_iterator j = orderedModels.rbegin(); j != orderedModels.rend(); j++){
		model modCopy = *j;
		SET_VECTOR_ELT(ret, i++, modCopy.convert2list(currentFpInfo, logMargLikConst, logPriorConst, normConst));
	}
	Rf_setAttrib(ret, Rf_install("numVisited"), Rf_ScalarReal(bookkeep.modelCounter));
	Rf_setAttrib(ret, Rf_install("inclusionProbs"), inc);
	Rf_setAttrib(ret, Rf_install("logNormConst"), Rf_ScalarReal(logl (normConst)));

	// return ###
	Rf_unprotect(nProtect);
	return ret;
}

// ***************************************************************************************************//

// recursion via:
void permPars(
			unsigned int pos, // current position in parameter vector, starting from 0 - copied.
			const fpInfo& currFp,
			const int &nUcGroups,
			modelPar mod,	// is copied every time! everything else is call by reference.
			set<model> &space,
			const hyperPriorPars &hyp,
			const dataValues &data,
			set<int>* const &ucTermList,
			const set<int> &fixedCols,
			book &bookkeep
			)
{
	if (pos != currFp.nFps){ // some fps are still left
		const int card = currFp.fpcards[pos]; // cardinality of this power set
		permPars(pos + 1, currFp, nUcGroups, mod, space, hyp, data, ucTermList, fixedCols, bookkeep); // degree 0
		for (int deg = 1; deg <= currFp.fpmaxs[pos]; deg++){ // different degrees for fp at pos
			mod.fpSize++; // increment sums of fp degrees
			int part[card]; // partition of deg into card parts
			bool more1 = false;
			int h(0), t(0); // internal variables for comp_next
			do {
				comp_next(deg, card, part, &more1, h, t);	// next partition of deg into card parts
				mod.fpPars[pos] = freqvec2multiset(part, card); // convert into multiset
				// and go on
				permPars(pos + 1, currFp, nUcGroups, mod, space, hyp, data, ucTermList, fixedCols, bookkeep);
			} while (more1);
		}
	} else { // no fps left
		computeModel(mod, hyp, data, currFp, ucTermList, nUcGroups, fixedCols, space, bookkeep);
		for (int deg = 1; deg <= nUcGroups; deg++){ // different number of uc groups
			mod.ucSize++; // increment number of uc groups present
			int subset[deg]; // partition of deg into card parts
			bool more2 = false;
			int m(0), m2(0); // internal variables for ksub_next
			do {
				ksub_next(nUcGroups, deg, subset, &more2, m, m2);	// next subset (positive integers)
				mod.ucPars = set<int>(subset, subset + deg); // convert into set
				computeModel(mod, hyp, data, currFp, ucTermList, nUcGroups, fixedCols, space, bookkeep);
			} while (more2);
		}
	}
}



// ***************************************************************************************************//

void computeModel(// compute (varying part of) marginal likelihood and prior of mod and insert into map
					const modelPar &mod,
					const hyperPriorPars &hyp,
					const dataValues &data,
					const fpInfo &currFp,
					set<int>* const &ucTermList,
					const int &nUcGroups,
					const set<int> &fixedCols,
					set<model> &space,
					book &bookkeep
				 )
{
	static set<model>::size_type compCounter = 0;

	// design matrix
	Matrix thisDesign = getDesignMatrix(mod, data, currFp, ucTermList, nUcGroups, fixedCols);

	// R2
	double thisR2 = getR2(thisDesign, data, fixedCols, hyp);

	if (! ISNAN(thisR2)){
		// log marginal likelihood
		double thisVarLogMargLik = getVarLogMargLik(thisR2, data.nObs, thisDesign.Ncols(), hyp);

		// log prior
		const double thisLogPrior = getVarLogPrior(mod, currFp, hyp);

		// posterior expected g
		double thisPostExpectedg = posteriorExpectedg_hyperg(thisR2, data.nObs, thisDesign.Ncols(), hyp.a, thisVarLogMargLik);

		// posterior expected shrinkage
		double thisPostExpectedShrinkage = posteriorExpectedShrinkage_hyperg(thisR2, data.nObs, thisDesign.Ncols(), hyp.a, thisVarLogMargLik);

		// put all this into the modelInfo
		modelInfo info(thisVarLogMargLik, thisLogPrior, thisPostExpectedg, thisPostExpectedShrinkage, thisR2);

		// altogether we have the model:
		model thisModel = model(mod, info);

		// and insert it into the model space
		if (space.size() >= bookkeep.nModels){
			set<model>::iterator it = space.begin();
			if (*it < thisModel){ // compare this model to the least probable model in the set
				space.erase(it);
				space.insert(thisModel); // exchange if it is better than this worst model in the set
			}
		} else {
			space.insert(thisModel);
		}

		const long double thisPropToPosterior = expl(thisVarLogMargLik + thisLogPrior);
		bookkeep.modelPropToPosteriors.add(thisPropToPosterior);

		pushInclusionProbs(mod, currFp, nUcGroups, bookkeep);
		bookkeep.modelCounter++;

	} else {
		bookkeep.nanCounter++;
	}
	// increase static vars
	if((++compCounter % max(data.totalNumber / 100, static_cast<dataValues::NumberType>(1)) == 0) && bookkeep.verbose)
		Rprintf("-"); // display computation progress at each percent
}

// ***************************************************************************************************//

SEXP samplingGaussian(// definition
				SEXP R_x, 				// design matrix (with colnames)
				SEXP R_xcentered, 		// centered design matrix
                SEXP R_y, 				// response vector
                SEXP R_fpmaxs, 			// vector of maximum fp degrees
                SEXP R_fppos,			// corresponding vector of fp column indices
	            SEXP R_fpcards,			// corresponding vector of power set cardinalities
                SEXP R_nFps, 			// number of fp terms
                SEXP R_fpnames,			// names of fp terms
                SEXP R_ucIndices, 		// vector giving _unc_ertainty custer indices (column -> which group)
                SEXP R_ucTermList,   	// list for group -> which columns mapping
                SEXP R_nUcGroups, 		// number of uncertainty groups
                SEXP R_hyperparam,		// hyperparameter a for hyper-g prior
                SEXP R_useSparsePrior, // use sparse model prior?
                SEXP R_nModels,			// number of best models to be returned
                SEXP R_verbose,			// should progress been displayed?
                SEXP R_chainlength		// how many times should a jump been made?
						)
{
	// important!!! We now assume that all elements of R_fpmaxs are identical!!!
	// It would be best to remove the option supporting different maximum FP degrees from the code,
	// to be inline with the paper.

	unsigned int nProtect = 0;
	// unpack ###
	// data
	const Matrix x = getMatrix(R_x);
	const Matrix xcentered = getMatrix(R_xcentered);
	const ColumnVector y = vec2col(R_y);

	dataValues data(x, xcentered, y, 0); // totalNumber is not needed

	// fp info
	fpInfo currentFpInfo;
	currentFpInfo.nFps = INTEGER(R_nFps)[0];
	currentFpInfo.fpmaxs = INTEGER(R_fpmaxs);
	currentFpInfo.fppos = INTEGER(R_fppos);
	currentFpInfo.fpcards = INTEGER(R_fpcards);
	Rf_protect(currentFpInfo.fpnames = R_fpnames);
	nProtect++;

	// the FP range
	const std::set<unsigned int> fpRange = constructSequence(currentFpInfo.nFps);

	// maximum fp dimension
	currentFpInfo.maxFpDim = accumulate(currentFpInfo.fpmaxs, currentFpInfo.fpmaxs + currentFpInfo.nFps, 0);

	// determine maximum powerset
	const int* biggestMaxDegree = max_element(currentFpInfo.fpmaxs, currentFpInfo.fpmaxs + currentFpInfo.nFps);
	double powerset[max(8, 5 + *biggestMaxDegree)];
//	const unsigned int fixedpowerIndices[] ={  0,  1,    2, 3,   4, 5, 6, 7 };
	const double fixedpowers[] = 			{ -2, -1, -0.5, 0, 0.5, 1, 2, 3 }; // always in powerset
	copy(fixedpowers, fixedpowers + 8, powerset);
	for(int more = 3; more < *biggestMaxDegree; more++){ // additional powers
		*(powerset + 8 + (3 - more)) = more + 1;
	}
	currentFpInfo.powerset = powerset; // save maximum powerset

	vector<ColumnVector> transformedCols[currentFpInfo.nFps]; // build array of vectors of ColumnVectors holding the required
									// transformed values for the design matrices
	for (unsigned int i = 0; i != currentFpInfo.nFps; i++){ // for every fp term
		const int nCols = currentFpInfo.fpcards[i];
		const ColumnVector thisCol = x.Column(currentFpInfo.fppos[i]);
		vector<ColumnVector> thisFp;
		for (int j = 0; j != nCols; j++){ // for every possible power
			ColumnVector thisTransform = thisCol;
			double thisPower = currentFpInfo.powerset[j];
			if(thisPower){ // not 0
				for (int k = 0; k != thisTransform.Nrows(); k++){ // transform each element
					assert(thisTransform.element(k) > 0);
					thisTransform.element(k) = pow (thisTransform.element(k), thisPower);
					assert(! ISNAN(thisTransform.element(k)));
				}
			} else { // 0
				for (int k = 0; k != thisTransform.Nrows(); k++){
					assert(thisTransform.element(k) > 0);
					thisTransform.element(k) = log (thisTransform.element(k));
					assert(! ISNAN(thisTransform.element(k)));
				}
			}
			// do not! center the column. This is done inside getFpMatrix, because
			// the repeated powers case cannot be treated here!!

			// and put it into vector of columns
			thisFp.push_back(thisTransform);
		}
		transformedCols[i] = thisFp;
	}
	currentFpInfo.tcols = transformedCols;

	// uc info
	const int* ucIndicesArray = INTEGER(R_ucIndices);
	const vector<int> ucIndices(ucIndicesArray, ucIndicesArray + Rf_length(R_ucIndices));
	const int nUcGroups = INTEGER(R_nUcGroups)[0];

	vector<unsigned int> ucSizes;
	for (int i = 1; i <= nUcGroups; i++){
		ucSizes.push_back(count(ucIndices.begin(), ucIndices.end(), i));
	}
	int maxUcDim = accumulate(ucSizes.begin(), ucSizes.end(), 0);


	set<int> ucTermList[nUcGroups]; // Array with length nUcGroups of int-vectors
	if(nUcGroups){ // catch case with no uc groups
		for(R_len_t i = 0; i != Rf_length(R_ucTermList); i++){
			SEXP temp = VECTOR_ELT(R_ucTermList, i);
			copy(INTEGER(temp), INTEGER(temp) + Rf_length(temp), inserter(ucTermList[i], ucTermList[i].begin()));
		}
	}

	// determine columns that are always in the design matrix
	set<int> nonfixedCols = set<int>(currentFpInfo.fppos, currentFpInfo.fppos + currentFpInfo.nFps); // fixed info
	set<int> allCols, fixedCols;
	for(vector<int>::size_type i = 0; i != ucIndices.size(); i++){
		allCols.insert(i + 1);
		if (ucIndices.at(i))
			nonfixedCols.insert(i + 1);
	}
	set_difference(allCols.begin(), allCols.end(), nonfixedCols.begin(), nonfixedCols.end(),
				   inserter(fixedCols, fixedCols.begin()));
	unsigned int fixedDim = fixedCols.size();

	// now fixedCols contains indices of columns that are always present in the design matrix

	// prior specifications
	const double hyperparam = Rf_asReal(R_hyperparam);
	const bool useSparsePrior = Rf_asLogical(R_useSparsePrior);

	hyperPriorPars hyp(hyperparam,
	                   useSparsePrior);

	// in case of a linear model the best models can be found during chain run and put in here:
	set<model> orderedModels;

	// bookkeeping:
	book bookkeep; // 0) initializes empty sum of prop to posteriors and modelCounter 0

	// a) length of chain
	double chainlength = REAL(R_chainlength)[0];
	if (ULLONG_MAX < chainlength){
		Rf_warning("\nchainlength too high - reducing to %d \n", ULLONG_MAX);
		bookkeep.chainlength = ULLONG_MAX;
	} else {
		bookkeep.chainlength = static_cast<unsigned long long int>(chainlength);
	}


	// b) verbose?
	bookkeep.verbose = LOGICAL(R_verbose)[0];

	// how many models to return?
	if (Rf_length(R_nModels) == 0){
		bookkeep.nModels = max(bookkeep.chainlength * 1 / 100, static_cast<unsigned long long int>(1));
	} else {
		bookkeep.nModels = INTEGER(R_nModels)[0];
	}

	// for computation of inclusion probs
	indexSafeSum cgwp[currentFpInfo.nFps + nUcGroups];
	bookkeep.covGroupWisePosteriors = cgwp;


	// models will be put in here:
	map<modelPar, modelInfo> modelSpace;

	// upper limit for num of columns
	unsigned int maxDim = min(static_cast<unsigned int>(data.nObs), fixedDim + currentFpInfo.maxFpDim + maxUcDim);

	// theta and theta':
	modelmcmc old, now;

	// start model
	modelPar startModel(currentFpInfo.nFps, 0, 0);
	powervecType startFps(currentFpInfo.nFps); // initialize empty vector of correct length
	startModel.fpPars = startFps;
	old.modPar = startModel;
	old.dim = fixedDim;
	old.freeUcs = getFreeUcs(old.modPar, ucSizes, old.dim, maxDim);
	old.freeCovs = getFreeCovs(old.modPar, currentFpInfo, old.freeUcs, old.dim, maxDim);
	old.presentCovs = getPresentCovs(old.modPar);
	old.birthprob = 1; old.deathprob = old.moveprob = 0;

	Matrix oldDesign = getDesignMatrix(old.modPar, data, currentFpInfo, ucTermList, nUcGroups, fixedCols);
	double oldR2 = getR2(oldDesign, data, fixedCols, hyp);

   // log marginal likelihood
	old.logMargLik = getVarLogMargLik(oldR2, data.nObs, oldDesign.Ncols(), hyp);

	// posterior expected g
	double oldPostExpectedg = posteriorExpectedg_hyperg(oldR2, data.nObs, oldDesign.Ncols(), hyp.a, old.logMargLik);

	// posterior expected shrinkage
	double oldPostExpectedShrinkage = posteriorExpectedShrinkage_hyperg(oldR2, data.nObs, oldDesign.Ncols(), hyp.a, old.logMargLik);

	// insert this model into map container
	double logPrior = getVarLogPrior(old.modPar, currentFpInfo, hyp);
	modelInfo startInfo(old.logMargLik, logPrior, oldPostExpectedg, oldPostExpectedShrinkage, oldR2, 1);
	modelSpace[old.modPar] = startInfo;
	old.mapPos = modelSpace.find(old.modPar);
	// and into the ordered models set
	orderedModels.insert(model(old.modPar, startInfo));

	// bookkeeping
	long double propToPosterior = expl(old.logMargLik + logPrior);
	bookkeep.modelPropToPosteriors.add(propToPosterior);

	pushInclusionProbs(old.modPar, currentFpInfo, nUcGroups, bookkeep);
	bookkeep.modelCounter++;

	now = old;

	// Start MCMC sampler***********************************************************//
	GetRNGstate(); // use R's random number generator
	for(unsigned long long int t = 0; t != bookkeep.chainlength; /* ++t explicitly at the end */){
		double logR; // log(prior times proposal ratio)
		// randomly select move type
		double u1 = unif_rand();
		if (u1 < old.birthprob){											// BIRTH
			unsigned int newCovInd = discreteUniform<unsigned int>(old.freeCovs);
			if (newCovInd <= currentFpInfo.nFps){ 					// some fp index
				int powerIndex = discreteUniform(0, currentFpInfo.fpcards[newCovInd-1]);
				now.modPar.fpPars.at(newCovInd-1).insert(powerIndex);
				now.modPar.fpSize++; // correct invariants
				now.dim++;
				unsigned int newPowersEqualPowerIndex = count(now.modPar.fpPars.at(newCovInd-1).begin(), now.modPar.fpPars.at(newCovInd-1).end(), powerIndex);
				unsigned int m = old.modPar.fpPars.at(newCovInd-1).size();
				logR = hyp.useSparsePrior ?
				        log(newPowersEqualPowerIndex) + log(currentFpInfo.fpcards[newCovInd-1]) - log(currentFpInfo.fpcards[newCovInd-1] + m)
				        : log(newPowersEqualPowerIndex) + log(currentFpInfo.fpcards[newCovInd-1]) - log1p(m);
			} else { 													// uc index
				int index = discreteUniform<int>(old.freeUcs);
				now.modPar.ucPars.insert(index);
				now.modPar.ucSize++;
				now.dim += ucSizes.at(index - 1);
				now.freeUcs = getFreeUcs(now.modPar, ucSizes, now.dim, maxDim);
				logR = log(old.freeUcs.size()) - log(now.modPar.ucSize);
			}
			now.presentCovs.insert(newCovInd);
			now.freeCovs = getFreeCovs(now.modPar, currentFpInfo, now.freeUcs, now.dim, maxDim);
			if (now.dim == maxDim){
				now.birthprob = 0; now.deathprob = now.moveprob = (now.modPar.fpSize > 0) ? 1.0 / 3 : 0.5;
			} else {
				now.birthprob = now.deathprob =	now.moveprob = (now.modPar.fpSize > 0) ? 0.25 : 1.0 / 3;
			}
			logR += log(now.deathprob) - log(old.birthprob) + log(old.freeCovs.size()) - log(now.presentCovs.size());
		} else if (u1 < old.birthprob + old.deathprob){					// DEATH
			unsigned int oldCovInd = discreteUniform<unsigned int>(old.presentCovs);
			if (oldCovInd <= currentFpInfo.nFps){ 					// some fp index
				multiset<int>::iterator powerIterator = dU<multiset<int> >(now.modPar.fpPars.at(oldCovInd-1));
				unsigned int oldPowersEqualPowerIndex = count(old.modPar.fpPars.at(oldCovInd-1).begin(), old.modPar.fpPars.at(oldCovInd-1).end(), *powerIterator);
				now.modPar.fpPars.at(oldCovInd-1).erase(powerIterator);
				now.modPar.fpSize--; // correct invariants
				now.dim--;
				logR = hyp.useSparsePrior ?
				        - log(oldPowersEqualPowerIndex) + log(currentFpInfo.fpcards[oldCovInd-1] + now.modPar.fpPars.at(oldCovInd-1).size()) - log(currentFpInfo.fpcards[oldCovInd-1])
				        : - log(oldPowersEqualPowerIndex) - log(currentFpInfo.fpcards[oldCovInd-1]) + log(old.modPar.fpPars.at(oldCovInd-1).size());
			} else { 													// uc index
				set<int>::iterator IndIterator = dU<set<int> >(now.modPar.ucPars);
				now.modPar.ucSize--;
				now.dim -= ucSizes.at(*IndIterator - 1);
				now.modPar.ucPars.erase(IndIterator);
				now.freeUcs = getFreeUcs(now.modPar, ucSizes, now.dim, maxDim);
				logR = log(old.modPar.ucSize) - log(now.freeUcs.size());
			}
			now.presentCovs = getPresentCovs(now.modPar);
			now.freeCovs = getFreeCovs(now.modPar, currentFpInfo, now.freeUcs, now.dim, maxDim);
			if (now.dim == fixedDim){
				now.birthprob = 1; now.deathprob = now.moveprob = 0;
			} else {
				now.birthprob = now.deathprob =	now.moveprob = (now.modPar.fpSize > 0) ? 0.25 : 1.0 / 3;
			}
			logR += log(now.birthprob) - log(old.deathprob) + log(old.presentCovs.size()) - log(now.freeCovs.size());

		} else if (u1 < old.birthprob + old.deathprob + old.moveprob){	 // MOVE
			unsigned int CovInd = discreteUniform<unsigned int>(old.presentCovs);
			if (CovInd <= currentFpInfo.nFps){ 						// some fp index
				multiset<int>::iterator powerIterator = dU<multiset<int> >(now.modPar.fpPars.at(CovInd-1));
				unsigned int oldPowersEqualPowerIndex = count(old.modPar.fpPars.at(CovInd-1).begin(), old.modPar.fpPars.at(CovInd-1).end(), *powerIterator);
				now.modPar.fpPars.at(CovInd-1).erase(powerIterator);
				int powerIndex = discreteUniform(0, currentFpInfo.fpcards[CovInd-1]);
				now.modPar.fpPars.at(CovInd-1).insert(powerIndex);
				unsigned int newPowersEqualPowerIndex = count(now.modPar.fpPars.at(CovInd-1).begin(), now.modPar.fpPars.at(CovInd-1).end(), powerIndex);
				// free, present Covs and move type probs are unchanged
				logR = log(newPowersEqualPowerIndex) - log(oldPowersEqualPowerIndex);

			} else { 													// uc index
				set<int>::iterator IndIterator = dU<set<int> >(now.modPar.ucPars);
				now.modPar.ucSize--;
				now.dim -= ucSizes.at(*IndIterator - 1);
				now.modPar.ucPars.erase(IndIterator);
				now.freeUcs = getFreeUcs(now.modPar, ucSizes, now.dim, maxDim);
				int index = discreteUniform<int>(now.freeUcs);
				now.modPar.ucPars.insert(index);
				now.modPar.ucSize++;
				now.dim += ucSizes.at(index - 1);
				now.freeUcs = getFreeUcs(now.modPar, ucSizes, now.dim, maxDim);
				// here something may change, therefore:
				now.freeCovs = getFreeCovs(now.modPar, currentFpInfo, now.freeUcs, now.dim, maxDim);
				if (now.dim == maxDim){
					now.birthprob = 0; now.deathprob = now.moveprob = (now.modPar.fpSize > 0) ? 1.0 / 3 : 0.5;
				} else {
					now.birthprob = now.deathprob =	now.moveprob = (now.modPar.fpSize > 0) ? 0.25 : 1.0 / 3;
				}
				logR = 0;
			}
		} else {													// SWITCH (of FP vectors)
			// select only the FP present covs
			std::set<unsigned int> presentFps = removeElement(old.presentCovs, currentFpInfo.nFps + 1);

			// so we have the first power vector:
			unsigned int firstFpInd = discreteUniform<unsigned int>(presentFps);
			powers first = now.modPar.fpPars.at(firstFpInd - 1);

			// the second power vector from all other FPs
			std::set<unsigned int> otherFps = removeElement(fpRange, firstFpInd);
			unsigned int secondFpInd = discreteUniform<unsigned int>(otherFps);
			powers second = now.modPar.fpPars.at(secondFpInd - 1);

			// save the first
			powers saveFirst = first;

			// copy second to first
			now.modPar.fpPars.at(firstFpInd - 1) = second;

			// and save to second
			now.modPar.fpPars.at(secondFpInd - 1) = saveFirst;

			// so now we have switched the power vectors.

			// move type probs are not changed, because the number of present FPs is unchanged,
			// as well as the dimension of the model.

			// but carefully update the information which covariates are free and which are present
			now.freeCovs = getFreeCovs(now.modPar, currentFpInfo, now.freeUcs, now.dim, maxDim);
			now.presentCovs = getPresentCovs(now.modPar);

			// and the prior, proposal ratios are 1, thus:
			logR = 0;
		}

		double u2 = unif_rand();
		now.mapPos = modelSpace.find(now.modPar); // search for proposed model

		if (now.mapPos == modelSpace.end()){ 						// "now" is a new model

			// construct design matrix and compute R^2
			Matrix nowDesign = getDesignMatrix(now.modPar, data, currentFpInfo, ucTermList, nUcGroups, fixedCols);
			double nowR2 = getR2(nowDesign, data, fixedCols, hyp);

			if(ISNAN(nowR2)){ // check if new model is OK, if not then nan
				modelSpace[now.modPar] = modelInfo(NA_REAL, logPrior, 0, 0, 0);
				bookkeep.nanCounter++;
				old.mapPos->second.hits++;
				now = old;
			} else { 					// OK: only then can it possibly be accepted

				// log marginal likelihood and log Bayes factor
				now.logMargLik = getVarLogMargLik(nowR2, data.nObs, nowDesign.Ncols(), hyp);
				double logBayesFactor = now.logMargLik - old.logMargLik;

				double logAcceptanceProb = logBayesFactor + logR;
				logPrior = getVarLogPrior(now.modPar, currentFpInfo, hyp);

				// posterior expected g
				double nowPostExpectedg = posteriorExpectedg_hyperg(nowR2, data.nObs, nowDesign.Ncols(), hyp.a, now.logMargLik);

				// posterior expected shrinkage
				double nowPostExpectedShrinkage = posteriorExpectedShrinkage_hyperg(nowR2, data.nObs, nowDesign.Ncols(), hyp.a, now.logMargLik);

				// construct model
				model thisModel = model(now.modPar, modelInfo(now.logMargLik, logPrior, nowPostExpectedg, nowPostExpectedShrinkage, nowR2));

				// save model posteriors
				propToPosterior = expl(now.logMargLik + logPrior);
				bookkeep.modelPropToPosteriors.add(propToPosterior);

				pushInclusionProbs(thisModel.par, currentFpInfo, nUcGroups, bookkeep);
				bookkeep.modelCounter++;

				// insert it into orderedModels
				if (orderedModels.size() >= bookkeep.nModels){  // insert into best set?
					set<model>::iterator it = orderedModels.begin();
					if (*it < thisModel){ // compare this model to the least probable model in the set
						orderedModels.erase(it);
						orderedModels.insert(thisModel); // exchange if it is better than this worst model in the set
					}
				} else {
					orderedModels.insert(thisModel);
				}
				// and modelSpace
				modelSpace[now.modPar] = thisModel.info;
				if (u2 <= exp(logAcceptanceProb)){ 	// acceptance
					now.mapPos = modelSpace.find(now.modPar);
					now.mapPos->second.hits++;
					old = now;
				} else {							//  rejection
					old.mapPos->second.hits++;
					now = old;
				}
			}
		} else {												// "now" was proposed before
			now.logMargLik = now.mapPos->second.logMargLik;
			if(ISNAN(now.logMargLik)){
				old.mapPos->second.hits++;
				now = old;
			} else {
				double logBayesFactor = now.logMargLik - old.logMargLik;
				double logAcceptanceProb = logBayesFactor + logR;
				if (u2 <= exp(logAcceptanceProb)){ 	// acceptance
					now.mapPos->second.hits++;
					old = now;
				} else {							//  rejection
					old.mapPos->second.hits++;
					now = old;
				}
			}
		}
		if((++t % max(bookkeep.chainlength / 100, static_cast<long long unsigned int>(1)) == 0) && bookkeep.verbose)
			Rprintf("-"); // display computation progress at each percent
	}
	PutRNGstate(); // no RNs required anymore

	if (bookkeep.verbose){
		Rprintf("\nNumber of models explored:                     %d", bookkeep.modelCounter);
		Rprintf("\nNumber of non-identifiable models encountered: %d", bookkeep.nanCounter);
		Rprintf("\nNumber of saved models:                        %d\n", orderedModels.size());
	}

	// normalize posterior probabilities and correct log marg lik and log prior
	const long double normConst = bookkeep.modelPropToPosteriors.sum();
	const double logMargLikConst = 	lgammafn((data.nObs - 1) / 2.0) -
									(data.nObs - 1) * sqrt(data.sumOfSquaresTotal) -
									(data.nObs - 1) * M_LN_SQRT_PI -
									0.5 * log(data.nObs);
	const double logPriorConst = nUcGroups * M_LN2; // note: M_LN2 = log(2)

	// return best models
	SEXP ret;
	Rf_protect(ret = Rf_allocVector(VECSXP, orderedModels.size()));
	nProtect++;

	// inclusion probs
	SEXP inc;
	Rf_protect(inc = Rf_allocVector(REALSXP, currentFpInfo.nFps + nUcGroups));
	nProtect++;
	for (int i = 0; i != Rf_length(inc); i++)
		REAL(inc)[i] = bookkeep.covGroupWisePosteriors[i].sum(bookkeep.modelPropToPosteriors) / normConst;

	// build return list
	unsigned long long int i = 0;
	for(set<model>::reverse_iterator j = orderedModels.rbegin(); j != orderedModels.rend(); j++){
		model modCopy = *j;
		map<modelPar, modelInfo>::const_iterator modPosInMap = modelSpace.find(modCopy.par); // plug in correct hit value
		modCopy.info.hits = modPosInMap->second.hits;
		SET_VECTOR_ELT(ret, i++, modCopy.convert2listMcmc(currentFpInfo, logMargLikConst, logPriorConst, normConst, bookkeep)); // and insert into return list
	}
	Rf_setAttrib(ret, Rf_install("numVisited"), Rf_ScalarReal(bookkeep.modelCounter));
	Rf_setAttrib(ret, Rf_install("inclusionProbs"), inc);
	Rf_setAttrib(ret, Rf_install("logNormConst"), Rf_ScalarReal(logl (normConst)));

	// return ###
	Rf_unprotect(nProtect);
	return ret;
}

// ***************************************************************************************************//

set<int> getFreeUcs(	// compute set of free uc group indices
					const modelPar& mod,
					const vector<unsigned int>& ucSizes,
					const unsigned int& currDim,
					const unsigned int& maxDim
					)
{
	set<int> ret;
	for (int i = 1; i <= static_cast<int>(ucSizes.size()); i++){ // for every uc index
		if ((find(mod.ucPars.begin(), mod.ucPars.end(), i) == mod.ucPars.end()) && (ucSizes.at(i-1) <= maxDim - currDim))
			ret.insert(i); // insert if not already in model and enough space in design matrix
	}
	return ret;
}


// ***************************************************************************************************//

set<unsigned int> getFreeCovs(					// compute set of free cov indices
				const modelPar& mod,
				const fpInfo& currFp,
				const set<int>& freeUcs,
				const unsigned int& currDim,
				const unsigned int& maxDim
					)
{
	set<unsigned int> ret;

	if (currDim == maxDim)
		return ret;

	for (unsigned int i = 0; i != mod.nFps; i++){
		if (mod.fpPars.at(i).size() < static_cast<unsigned int>(currFp.fpmaxs[i]))
			ret.insert(i + 1);
	}

	if (! freeUcs.empty())
		ret.insert(mod.nFps + 1);

	return ret;
}

// ***************************************************************************************************//


set<unsigned int> getPresentCovs	( // determine set of present cov indices
					const modelPar& mod
						)
{
	set<unsigned int> ret;
	for (unsigned int i = 0; i != mod.nFps; i++){
		if (! mod.fpPars.at(i).empty())
			ret.insert(i+1);
	}
	if (! mod.ucPars.empty())
		ret.insert(mod.nFps + 1);
	return ret;
}

// ***************************************************************************************************//


template <class T>
T discreteUniform (	// return random element of myset; should be enclosed in getRNGstate() etc.
				const set<T>& myset
					)
{
	if (myset.empty())
		Rf_error("\nmyset is empty!\n");

	double u = unif_rand();
	typename set<T>::size_type size = myset.size();
	typename set<T>::const_iterator i = myset.begin(); typename set<T>::size_type j = 1;
	while(u > 1.0 / size * j){
		i++; j++;
	}
	return *i;
}

// ***************************************************************************************************//

template <class T>
typename T::iterator dU (	// return iterator of random element of myset; should be enclosed in getRNGstate() etc.
				const T& container
					)
{
	if (container.empty())
		Rf_error("\ncontainer is empty!\n");

	double u = unif_rand();
	typename T::size_type size = container.size();
	typename T::iterator i = container.begin(); typename T::size_type j = 1;
	while(u > 1.0 / size * j){
		i++; j++;
	}
	return i;
}


// ***************************************************************************************************//

int discreteUniform ( // get random int x with lower <= x < upper; should be enclosed in getRNGstate() etc.
						const int& lower,
						const int& upper
					)
{
	if (lower >= upper)
		Rf_error("\nlower = %d >= %d = upper in discreteUniform call\n", lower, upper);

	int size = upper - lower;
	int ret = lower;
	double u = unif_rand();

	while(u > 1.0 / size * (ret - lower + 1)){
		ret++;
	}
	return ret;
}

// ***************************************************************************************************//

ReturnMatrix getDesignMatrix ( // construct centered design matrix including intercept for the model
					const modelPar &mod,
					const dataValues &data,
					const fpInfo &currFp,
					set<int>* const &ucTermList,
					const int &nUcGroups,
					const set<int> &fixedCols
					)
{
//
//	MATRIXSTORE(data.response, responseStore)
//
	// build design matrix B

	// intercept column
	Matrix B = data.onesVector;

	// centered fp matrices
	for (unsigned int i = 0; i != currFp.nFps; i++){
		multiset<int> powersi = mod.fpPars[i];
		if (! powersi.empty()){
			Matrix Fp = getFpMatrix(currFp.tcols[i], powersi, data); // this is centered
			B = B | Fp;
		}
	}

	// centered uc matrices
	for (int i = 0; i != nUcGroups; i++){
		set<int>::const_iterator ipos = find(mod.ucPars.begin(), mod.ucPars.end(), i + 1);
		if (ipos != mod.ucPars.end()){ // if mod.ucPars contains i
			Matrix Uc = getMultipleCols(data.centeredDesign, ucTermList[i]);	// this is centered
			B = B | Uc;
		}
	}

	B.release();
	return B;
}




// ***************************************************************************************************//

double getR2(	// compute coefficient of determination for the model
			const Matrix &design, // must be the centered design matrix!! (including the intercept)
			const dataValues &data,
			const set<int> &fixedCols,
			const hyperPriorPars &hyp
			)
{

	int dim = design.Ncols();
	if (dim - 1 >= data.nObs - 3 - hyp.a) return NA_REAL; // not a valid model
//
//	MATRIXSTORE(B, BStore)
//
	const int numFixedCols = fixedCols.size();

	if (dim == numFixedCols) { // then this is the null model
		return 0; // because SSE == SST in this case
	}
	else {	// else we need to work a bit

		// select nonfixed part of the centered design matrix, because the fixed part is not needed here
		Matrix X = design.Columns(numFixedCols + 1, dim);

//		double test;
//		MATRIXSTORE(X, XStore)

		// compute model-specific parts
		SymmetricMatrix XtX;
		XtX << X.t() * X;

		try // a cholesky decomposition of XtX
		{
		    LowerTriangularMatrix LeftRootOfXtX = Cholesky(XtX);

//		    MATRIXSTORE(LeftRootOfXtX, lroxtxStore)

		    // compute coefficient of determination R2
		    ColumnVector tmp = LeftRootOfXtX.i() * (X.t() * data.response);

//		    MATRIXSTORE(tmp, tmpStore)

		    double sumOfSquaresModel = tmp.sum_square();

			double R2 = sumOfSquaresModel / data.sumOfSquaresTotal;

			assert((R2 <= 1) & (R2 >= 0));

			return R2;
		}
		catch(NPDException) {return NA_REAL;} // if XtX is not p.d. then return NA
	}
}

// ***************************************************************************************************//

// compute varying part of log marginal likelihood for specific model
double getVarLogMargLik(const double &R2, const int &n, const int &dim, const hyperPriorPars &hyp)
{
	if(dim == 1){
		return 0;
	} else {

		double logBF = logBF_hyperg(R2, n, dim, hyp.a);
		assert(! ISNAN(logBF)); // each possible case has been treated above!

		return logBF;
	}
}

// ***************************************************************************************************//

double getVarLogPrior( // compute varying part of logarithm of model prior
                      const modelPar &mod,
                      const fpInfo &currFp,
                      const hyperPriorPars &hyp)
{
    if (hyp.useSparsePrior)
    {
        safeSum thisVarLogPrior;
        for (unsigned int i = 0; i != currFp.nFps; i++)
        { // for each fp covariate
            unsigned int degree = mod.fpPars[i].size();
            double thisVal = -lchoose(currFp.fpcards[i] - 1 + degree, degree)
                    -log1p(currFp.fpmaxs[i]);
            thisVarLogPrior.add(thisVal);
        }
        return thisVarLogPrior.sum();
    }
    else
    {
        return 0;
    }
}

// ***************************************************************************************************//

ReturnMatrix getFpMatrix( // build Fp basis matrix from vector, power indices and power set for one covariate
			const vector<ColumnVector> &tcols,
			const multiset<int> &powerinds,
			const dataValues &data
			)
{
	const int logInd = 3; // this index corresponds to power = 0, i.e. log.
	const int nrow = tcols.at(1).Nrows();

	Matrix ret(nrow, powerinds.size());

	// start recursion
	int lastInd = logInd;
	ColumnVector lastCol(nrow); lastCol = 1;

	// there is at least one power present
	multiset<int>::size_type cols = 1; // invariant: about to process column number cols

	for (multiset<int>::const_iterator now = powerinds.begin(); now != powerinds.end(); now++){
		if (*now == lastInd){ 	// repeated powers case
			lastCol = SP(lastCol, tcols[logInd]);
		} else {				// normal case
			lastInd = *now;
			lastCol = tcols[lastInd];
		}
		// center the column
		ret.Column(cols++) = lastCol - (lastCol.sum() / data.nObs) * data.onesVector;
	}

	ret.Release(); return ret;
}

// ***************************************************************************************************//

void pushInclusionProbs(	// push back index into covGroupWisePosteriors-Array
						const modelPar &mod,
						const fpInfo &currFp,
						const int &nUcGroups,
						book &bookkeep)
{
	for (unsigned int i = 0; i != currFp.nFps; i++){
		if (! mod.fpPars.at(i).empty())
			bookkeep.covGroupWisePosteriors[i].add(bookkeep.modelCounter);
	}

	for (int i = 1; i <= nUcGroups; i++){
		set<int>::const_iterator ipos = find(mod.ucPars.begin(), mod.ucPars.end(), i);
		if (ipos != mod.ucPars.end()){ // if mod.ucPars contains i
			bookkeep.covGroupWisePosteriors[i - 1 + currFp.nFps].add(bookkeep.modelCounter);
		}
	}
}



// ***************************************************************************************************//

// this is an interface for R to the log marg lik computation
SEXP logMargLik( //definition
				SEXP R_R2,				// coefficient of determination
				SEXP R_n,				// number of observations
				SEXP R_dim,				// number of columns of the design matrix
				SEXP R_alpha, 			// hyperparamater for hyper-g prior
				SEXP R_sst				// total sum of squares computed from y
				)
{
	unsigned int nProtect = 0;

	// unpack
	const double R2 = REAL(R_R2)[0];
	const int n = INTEGER(R_n)[0];
	const int dim = INTEGER(R_dim)[0];
	const double alpha = REAL(R_alpha)[0];
	const double sst = REAL(R_sst)[0];

	// compute
	hyperPriorPars hyp(alpha, 0); // value of useSparsePrior does not matter here
	double varLogMargLik = getVarLogMargLik(R2, n, dim, hyp);
	double logMargLikConst = lgammafn((n - 1) / 2.0) -
							(n - 1) * sqrt(sst) -
							(n - 1) * M_LN_SQRT_PI -
							0.5 * log(n);

	SEXP ret;
	Rf_protect(ret = Rf_ScalarReal(varLogMargLik + logMargLikConst));
	nProtect++;

	Rf_unprotect(nProtect);
	return(ret);
}

// ***************************************************************************************************//

// this is an interface for R to the computation of the posterior expected g
SEXP postExpectedg( // definition
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha) // hyperparamater for hyper-g prior
{
    unsigned int nProtect = 0;

    // unpack
    const double R2 = REAL(R_R2)[0];
    const int n = INTEGER(R_n)[0];
    const int dim = INTEGER(R_dim)[0];
    const double alpha = REAL(R_alpha)[0];

    // compute
    hyperPriorPars hyp(alpha, 0); // value of useSparsePrior does not matter here
    const double varLogMargLik = getVarLogMargLik(R2, n, dim, hyp);
    const double postExpectedg = posteriorExpectedg_hyperg(R2, n, dim, hyp.a, varLogMargLik);

    SEXP ret;
    Rf_protect(ret = Rf_ScalarReal(postExpectedg));
    nProtect++;

    Rf_unprotect(nProtect);
    return(ret);
}

// ***************************************************************************************************//

// this is an interface for R to the computation of the posterior expected shrinkage factor
SEXP postExpectedShrinkage( //definition
                SEXP R_R2, // coefficient of determination
                SEXP R_n, // number of observations
                SEXP R_dim, // number of columns of the design matrix
                SEXP R_alpha) // hyperparamater for hyper-g prior
{
    unsigned int nProtect = 0;

    // unpack
    const double R2 = REAL(R_R2)[0];
    const int n = INTEGER(R_n)[0];
    const int dim = INTEGER(R_dim)[0];
    const double alpha = REAL(R_alpha)[0];

    // compute
    hyperPriorPars hyp(alpha, 0); // value of useSparsePrior does not matter here
    const double varLogMargLik = getVarLogMargLik(R2, n, dim, hyp);
    const double postExpectedShrinkage = posteriorExpectedShrinkage_hyperg(R2, n, dim, hyp.a, varLogMargLik);

    SEXP ret;
    Rf_protect(ret = Rf_ScalarReal(postExpectedShrinkage));
    nProtect++;

    Rf_unprotect(nProtect);
    return(ret);
}

// ################################################################################################

int main() {} // dummy


// attention, when building shared library with R CMD SHLIB enter ALL source files that are entangled with
// this, too!! (RnewMat.cpp e.g.)
