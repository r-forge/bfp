#include "RnewMat.h"
#include "dataStructure.h"
#include <vector>
#include <numeric>
#include <algorithm>

using std::set;
using std::vector;
using std::accumulate;

ReturnMatrix getMatrix(const SEXP& m) // R-Matrix to Newmat-Matrix
{
	unsigned int mx, my;
	int *mDim;
	double *a;
	
	a = REAL(m); // read matrix m in array (columnwise)
	
	mDim = INTEGER(Rf_getAttrib(m, R_DimSymbol)); // get dimensions of m
	mx = mDim[0];
	my = mDim[1];
	
	Matrix M(my, mx); // copy array into matrix. This is done rowwise, so transpose necessary
	M << a;
	M = M.t();	
	
	M.Release(); return M;	
}

ColumnVector vec2col(const SEXP& v) // R-Vector to Newmat-column vector
{
	if (Rf_isMatrix(v))
		Rf_error("Argument of vec2col is a matrix\n");
		
	ColumnVector V(Rf_length(v));	
	double *a = REAL(v);
	V << a;
	
	return V;	
}

SEXP putMatrix(const Matrix& M) // Newmat-Matrix to R-Matrix
{
	unsigned int nProtected = 0, Mx, My, i, j;
	double *a;
	SEXP ret;
	
	Mx = M.Nrows(); // get dimensions of M
	My = M.Ncols();
	
	Rf_protect(ret = Rf_allocMatrix(REALSXP, Mx, My)); // allocate return matrix
	++nProtected;	
	a = REAL(ret);
	
	for (i = 0; i != Mx; i++){	// put values in return matrix
		for (j = 0; j != My; j++){
			a[i + j * Mx] = M.element(i,j);	
		}	
	}
		
	Rf_unprotect(nProtected); // unprotect everything
	return ret;
} 

ReturnMatrix getMultipleCols(const Matrix& M, const set<int>& s) // get different concatenated columns of matrix 
{
	
	
	// MATRIXSTORE(M, MStore)
	
	Matrix ret(M.Nrows(), s.size());

	set<int>::size_type cols = 1; // invariant: about to process column number cols
	for (set<int>::const_iterator i = s.begin(); i != s.end(); i++){
		ret.Column(cols++) = M.Column(*i);	
	}

	ret.Release(); return ret;
}

multiset<int> freqvec2multiset(int* const &vec, const int &vecLength) // convert frequency vector into multiset
{ 
	multiset<int> ret;
	for (int power = 0; power != vecLength; power++){ 
		for (int times = 0; times != vec[power]; times++)
			ret.insert(power);					
	}	
	return ret;
}

//double quadraticForm(const Matrix& M, const ColumnVector& v) // compute quadratic form v.t() * M * v where M is symmetric
//{	
//	int dim = M.Nrows();
//	if (dim != M.Ncols() || dim != v.Nrows())
//		Rf_error("quadratic got inputs with wrong dimensions\n");
//	
//	safeSum ret;
//		
//	for (int i = 0; i != dim; i++){
//		ret.add(pow(v.element(i), 2) * M.element(i, i));
//	}
//
//	for (int i = 0; i != dim; i++){
//		for (int j = 0; j < i; j++){
//			ret.add(2 * v.element(i) * v.element(j) * M.element(i,j));	
//		}	
//	}
//	
//	return(ret.sum()); 	 
//}
