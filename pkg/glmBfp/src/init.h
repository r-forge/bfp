#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Declare C++ functions

/* .Call calls */
RcppExport SEXP cpp_bfgs(SEXP);
RcppExport SEXP cpp_coxfit(SEXP);
RcppExport SEXP cpp_evalZdensity(SEXP);
RcppExport SEXP cpp_glmBayesMfp(SEXP);
RcppExport SEXP cpp_optimize(SEXP);
RcppExport SEXP cpp_sampleGlm(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"glmBfp_predBMAcpp", (DL_FUNC) &glmBfp_predBMAcpp, 3},
  {"cpp_bfgs",         (DL_FUNC) &cpp_bfgs,         6},
  {"cpp_coxfit",       (DL_FUNC) &cpp_coxfit,       5},
  {"cpp_evalZdensity", (DL_FUNC) &cpp_evalZdensity, 6},
  {"cpp_glmBayesMfp",  (DL_FUNC) &cpp_glmBayesMfp,  7},
  {"cpp_optimize",     (DL_FUNC) &cpp_optimize,     4},
  {"cpp_sampleGlm",    (DL_FUNC) &cpp_sampleGlm,    9},
  {NULL, NULL, 0}
};

void R_init_glmBfp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}