// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// altra
NumericVector altra(NumericVector v, int n, NumericMatrix ind, int nodes);
RcppExport SEXP _TreeGuidedRegression_altra(SEXP vSEXP, SEXP nSEXP, SEXP indSEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type nodes(nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(altra(v, n, ind, nodes));
    return rcpp_result_gen;
END_RCPP
}
// findLambdaMax
double findLambdaMax(NumericVector v, int n, NumericMatrix ind, int nodes);
RcppExport SEXP _TreeGuidedRegression_findLambdaMax(SEXP vSEXP, SEXP nSEXP, SEXP indSEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type nodes(nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(findLambdaMax(v, n, ind, nodes));
    return rcpp_result_gen;
END_RCPP
}
// treeNorm
double treeNorm(NumericVector x, int n, NumericMatrix ind, int nodes);
RcppExport SEXP _TreeGuidedRegression_treeNorm(SEXP xSEXP, SEXP nSEXP, SEXP indSEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type nodes(nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(treeNorm(x, n, ind, nodes));
    return rcpp_result_gen;
END_RCPP
}
// eplb
NumericVector eplb(NumericVector v, int n, double z, double lambda0);
RcppExport SEXP _TreeGuidedRegression_eplb(SEXP vSEXP, SEXP nSEXP, SEXP zSEXP, SEXP lambda0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    rcpp_result_gen = Rcpp::wrap(eplb(v, n, z, lambda0));
    return rcpp_result_gen;
END_RCPP
}
// general_altra
NumericVector general_altra(NumericVector v, int n, NumericVector G, NumericMatrix ind, int nodes);
RcppExport SEXP _TreeGuidedRegression_general_altra(SEXP vSEXP, SEXP nSEXP, SEXP GSEXP, SEXP indSEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type nodes(nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(general_altra(v, n, G, ind, nodes));
    return rcpp_result_gen;
END_RCPP
}
// general_findLambdaMax
double general_findLambdaMax(NumericVector v, int n, NumericVector G, NumericMatrix ind, int nodes);
RcppExport SEXP _TreeGuidedRegression_general_findLambdaMax(SEXP vSEXP, SEXP nSEXP, SEXP GSEXP, SEXP indSEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type nodes(nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(general_findLambdaMax(v, n, G, ind, nodes));
    return rcpp_result_gen;
END_RCPP
}
// general_treeNorm
double general_treeNorm(NumericVector x, int n, NumericVector G, NumericMatrix ind, int nodes);
RcppExport SEXP _TreeGuidedRegression_general_treeNorm(SEXP xSEXP, SEXP nSEXP, SEXP GSEXP, SEXP indSEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type nodes(nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(general_treeNorm(x, n, G, ind, nodes));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _TreeGuidedRegression_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TreeGuidedRegression_altra", (DL_FUNC) &_TreeGuidedRegression_altra, 4},
    {"_TreeGuidedRegression_findLambdaMax", (DL_FUNC) &_TreeGuidedRegression_findLambdaMax, 4},
    {"_TreeGuidedRegression_treeNorm", (DL_FUNC) &_TreeGuidedRegression_treeNorm, 4},
    {"_TreeGuidedRegression_eplb", (DL_FUNC) &_TreeGuidedRegression_eplb, 4},
    {"_TreeGuidedRegression_general_altra", (DL_FUNC) &_TreeGuidedRegression_general_altra, 5},
    {"_TreeGuidedRegression_general_findLambdaMax", (DL_FUNC) &_TreeGuidedRegression_general_findLambdaMax, 5},
    {"_TreeGuidedRegression_general_treeNorm", (DL_FUNC) &_TreeGuidedRegression_general_treeNorm, 5},
    {"_TreeGuidedRegression_rcpp_hello_world", (DL_FUNC) &_TreeGuidedRegression_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_TreeGuidedRegression(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}