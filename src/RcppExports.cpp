// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// GSEN
SEXP GSEN(NumericMatrix X, NumericVector b, NumericVector O, NumericVector xx, NumericVector e, double L, double A, int p);
RcppExport SEXP bWGR_GSEN(SEXP XSEXP, SEXP bSEXP, SEXP OSEXP, SEXP xxSEXP, SEXP eSEXP, SEXP LSEXP, SEXP ASEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type O(OSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    Rcpp::traits::input_parameter< double >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    __result = Rcpp::wrap(GSEN(X, b, O, xx, e, L, A, p));
    return __result;
END_RCPP
}
// KMUP
SEXP KMUP(NumericMatrix X, NumericVector b, NumericVector xx, NumericVector E, NumericVector L, int p, double Ve, double pi);
RcppExport SEXP bWGR_KMUP(SEXP XSEXP, SEXP bSEXP, SEXP xxSEXP, SEXP ESEXP, SEXP LSEXP, SEXP pSEXP, SEXP VeSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    __result = Rcpp::wrap(KMUP(X, b, xx, E, L, p, Ve, pi));
    return __result;
END_RCPP
}
// KMUP2
SEXP KMUP2(NumericMatrix X, NumericVector b, NumericVector xx, NumericVector E, NumericVector L, int p, double Ve, double pi, IntegerVector ro);
RcppExport SEXP bWGR_KMUP2(SEXP XSEXP, SEXP bSEXP, SEXP xxSEXP, SEXP ESEXP, SEXP LSEXP, SEXP pSEXP, SEXP VeSEXP, SEXP piSEXP, SEXP roSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ro(roSEXP);
    __result = Rcpp::wrap(KMUP2(X, b, xx, E, L, p, Ve, pi, ro));
    return __result;
END_RCPP
}
