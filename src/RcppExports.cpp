// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// softmaxC
NumericVector softmaxC(NumericVector x);
RcppExport SEXP _iplsv_softmaxC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(softmaxC(x));
    return rcpp_result_gen;
END_RCPP
}
// exp_nlpostC
double exp_nlpostC(List Z_exp, NumericMatrix PHI, NumericMatrix THETA, NumericMatrix PSI, List docs, NumericVector Ns, double eta, double gamma, double beta);
RcppExport SEXP _iplsv_exp_nlpostC(SEXP Z_expSEXP, SEXP PHISEXP, SEXP THETASEXP, SEXP PSISEXP, SEXP docsSEXP, SEXP NsSEXP, SEXP etaSEXP, SEXP gammaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type Z_exp(Z_expSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PHI(PHISEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PSI(PSISEXP);
    Rcpp::traits::input_parameter< List >::type docs(docsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ns(NsSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_nlpostC(Z_exp, PHI, THETA, PSI, docs, Ns, eta, gamma, beta));
    return rcpp_result_gen;
END_RCPP
}
// g_enlpC
List g_enlpC(List Z_exp, NumericMatrix PHI, NumericMatrix THETA, NumericMatrix PSI, List docs, NumericVector Ns, double eta, double gamma, double beta);
RcppExport SEXP _iplsv_g_enlpC(SEXP Z_expSEXP, SEXP PHISEXP, SEXP THETASEXP, SEXP PSISEXP, SEXP docsSEXP, SEXP NsSEXP, SEXP etaSEXP, SEXP gammaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type Z_exp(Z_expSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PHI(PHISEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type THETA(THETASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PSI(PSISEXP);
    Rcpp::traits::input_parameter< List >::type docs(docsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ns(NsSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(g_enlpC(Z_exp, PHI, THETA, PSI, docs, Ns, eta, gamma, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_iplsv_softmaxC", (DL_FUNC) &_iplsv_softmaxC, 1},
    {"_iplsv_exp_nlpostC", (DL_FUNC) &_iplsv_exp_nlpostC, 9},
    {"_iplsv_g_enlpC", (DL_FUNC) &_iplsv_g_enlpC, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_iplsv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
