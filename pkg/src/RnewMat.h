#ifndef RNEWMAT_H_
#define RNEWMAT_H_

#include <R.h> // in /usr/lib/R/include. For R
#include <Rinternals.h>
#include <Rmath.h>

#include <newmat.h>  // For Newmat
#include <newmatap.h> // for cholesky decomposition

#include <set>   

#define MATRIXSTORE(m, mName)  volatile double *mName = m.Store(); test=mName[0];

ReturnMatrix getMatrix(const SEXP&); 	// R-Matrix to Newmat-Matrix
ColumnVector vec2col(const SEXP&);	// R-Vector to Newmat-column vector
SEXP putMatrix(const Matrix&); 	// Newmat-Matrix to R-Matrix
ReturnMatrix getMultipleCols(const Matrix&, const std::set<int>&); // get different concatenated columns of matrix 
multiset<int> freqvec2multiset(int* const&, const int&); // convert frequency vector into multiset
double quadraticForm(const Matrix&, const ColumnVector&); // quadratic form of symmetric matrix
 
#endif /*RNEWMAT_H_*/
