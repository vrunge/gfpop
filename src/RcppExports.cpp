// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gfpopTransfer
List gfpopTransfer(NumericVector vectData, DataFrame mygraph, std::string type, NumericVector vectWeight, bool testMode);
RcppExport SEXP _gfpop_gfpopTransfer(SEXP vectDataSEXP, SEXP mygraphSEXP, SEXP typeSEXP, SEXP vectWeightSEXP, SEXP testModeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vectData(vectDataSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type mygraph(mygraphSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vectWeight(vectWeightSEXP);
    Rcpp::traits::input_parameter< bool >::type testMode(testModeSEXP);
    rcpp_result_gen = Rcpp::wrap(gfpopTransfer(vectData, mygraph, type, vectWeight, testMode));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gfpop_gfpopTransfer", (DL_FUNC) &_gfpop_gfpopTransfer, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_gfpop(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
