// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// multiKernel
SEXP multiKernel(NumericMatrix Yr, NumericMatrix Zr, NumericMatrix Kr, double tau);
RcppExport SEXP multiKernel_multiKernel(SEXP YrSEXP, SEXP ZrSEXP, SEXP KrSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Yr(YrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Zr(ZrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Kr(KrSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    __result = Rcpp::wrap(multiKernel(Yr, Zr, Kr, tau));
    return __result;
END_RCPP
}
// multiKernel_noCon
SEXP multiKernel_noCon(NumericMatrix Yr, NumericMatrix Kr, double tau);
RcppExport SEXP multiKernel_multiKernel_noCon(SEXP YrSEXP, SEXP KrSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Yr(YrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Kr(KrSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    __result = Rcpp::wrap(multiKernel_noCon(Yr, Kr, tau));
    return __result;
END_RCPP
}
