#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _glmBfp_cpp_bfgs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmBfp_cpp_coxfit(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmBfp_cpp_evalZdensity(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmBfp_cpp_glmBayesMfp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmBfp_cpp_optimize(SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmBfp_cpp_sampleGlm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _glmBfp_predBMAcpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_glmBfp_cpp_bfgs",         (DL_FUNC) &_glmBfp_cpp_bfgs,         6},
    {"_glmBfp_cpp_coxfit",       (DL_FUNC) &_glmBfp_cpp_coxfit,       5},
    {"_glmBfp_cpp_evalZdensity", (DL_FUNC) &_glmBfp_cpp_evalZdensity, 7},
    {"_glmBfp_cpp_glmBayesMfp",  (DL_FUNC) &_glmBfp_cpp_glmBayesMfp,  7},
    {"_glmBfp_cpp_optimize",     (DL_FUNC) &_glmBfp_cpp_optimize,     4},
    {"_glmBfp_cpp_sampleGlm",    (DL_FUNC) &_glmBfp_cpp_sampleGlm,    9},
    {"_glmBfp_predBMAcpp",       (DL_FUNC) &_glmBfp_predBMAcpp,       3},
    {NULL, NULL, 0}
};

void R_init_glmBfp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
