// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// predBMAcpp
NumericMatrix predBMAcpp(NumericMatrix SurvMat, NumericMatrix LpMat, NumericVector WtVec);
RcppExport SEXP glmBfp_predBMAcpp(SEXP SurvMatSEXP, SEXP LpMatSEXP, SEXP WtVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type SurvMat(SurvMatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type LpMat(LpMatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type WtVec(WtVecSEXP);
    __result = Rcpp::wrap(predBMAcpp(SurvMat, LpMat, WtVec));
    return __result;
END_RCPP
}
