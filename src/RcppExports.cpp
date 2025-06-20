// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calcSESimpleSlopes
arma::vec calcSESimpleSlopes(arma::mat const& X, arma::mat const& V);
RcppExport SEXP _modsem_calcSESimpleSlopes(SEXP XSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(calcSESimpleSlopes(X, V));
    return rcpp_result_gen;
END_RCPP
}
// muLmsCpp
arma::vec muLmsCpp(Rcpp::List model, arma::vec z);
RcppExport SEXP _modsem_muLmsCpp(SEXP modelSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(muLmsCpp(model, z));
    return rcpp_result_gen;
END_RCPP
}
// sigmaLmsCpp
arma::mat sigmaLmsCpp(Rcpp::List model, arma::vec z);
RcppExport SEXP _modsem_sigmaLmsCpp(SEXP modelSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmaLmsCpp(model, z));
    return rcpp_result_gen;
END_RCPP
}
// completeLogLikLmsCpp
double completeLogLikLmsCpp(Rcpp::List modelR, Rcpp::List P, Rcpp::List quad);
RcppExport SEXP _modsem_completeLogLikLmsCpp(SEXP modelRSEXP, SEXP PSEXP, SEXP quadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type modelR(modelRSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type P(PSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type quad(quadSEXP);
    rcpp_result_gen = Rcpp::wrap(completeLogLikLmsCpp(modelR, P, quad));
    return rcpp_result_gen;
END_RCPP
}
// gradLogLikLmsCpp
arma::vec gradLogLikLmsCpp(const Rcpp::List& modelR, const Rcpp::List& P, const arma::uvec& block, const arma::uvec& row, const arma::uvec& col, double eps);
RcppExport SEXP _modsem_gradLogLikLmsCpp(SEXP modelRSEXP, SEXP PSEXP, SEXP blockSEXP, SEXP rowSEXP, SEXP colSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type modelR(modelRSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type block(blockSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type row(rowSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type col(colSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(gradLogLikLmsCpp(modelR, P, block, row, col, eps));
    return rcpp_result_gen;
END_RCPP
}
// muQmlCpp
arma::mat muQmlCpp(Rcpp::List m, int t);
RcppExport SEXP _modsem_muQmlCpp(SEXP mSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(muQmlCpp(m, t));
    return rcpp_result_gen;
END_RCPP
}
// sigmaQmlCpp
arma::mat sigmaQmlCpp(Rcpp::List m, int t);
RcppExport SEXP _modsem_sigmaQmlCpp(SEXP mSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmaQmlCpp(m, t));
    return rcpp_result_gen;
END_RCPP
}
// calcKronXi
arma::mat calcKronXi(Rcpp::List m, int t);
RcppExport SEXP _modsem_calcKronXi(SEXP mSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(calcKronXi(m, t));
    return rcpp_result_gen;
END_RCPP
}
// calcBinvCpp
arma::mat calcBinvCpp(Rcpp::List m, int t);
RcppExport SEXP _modsem_calcBinvCpp(SEXP mSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(calcBinvCpp(m, t));
    return rcpp_result_gen;
END_RCPP
}
// dnormCpp
arma::vec dnormCpp(const arma::vec& x, const arma::vec& mu, const arma::vec& sigma);
RcppExport SEXP _modsem_dnormCpp(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(dnormCpp(x, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// varZCpp
arma::mat varZCpp(arma::mat Omega, arma::mat Sigma1, int numEta);
RcppExport SEXP _modsem_varZCpp(SEXP OmegaSEXP, SEXP Sigma1SEXP, SEXP numEtaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma1(Sigma1SEXP);
    Rcpp::traits::input_parameter< int >::type numEta(numEtaSEXP);
    rcpp_result_gen = Rcpp::wrap(varZCpp(Omega, Sigma1, numEta));
    return rcpp_result_gen;
END_RCPP
}
// multiplyIndicatorsCpp
NumericVector multiplyIndicatorsCpp(DataFrame df);
RcppExport SEXP _modsem_multiplyIndicatorsCpp(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(multiplyIndicatorsCpp(df));
    return rcpp_result_gen;
END_RCPP
}
// rep_dmvnorm
arma::vec rep_dmvnorm(const arma::mat& x, const arma::mat& expected, const arma::mat& sigma, const int t);
RcppExport SEXP _modsem_rep_dmvnorm(SEXP xSEXP, SEXP expectedSEXP, SEXP sigmaSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type expected(expectedSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(rep_dmvnorm(x, expected, sigma, t));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm_arma_mc
arma::vec dmvnrm_arma_mc(arma::mat const& x, arma::rowvec const& mean, arma::mat const& sigma, bool const logd);
RcppExport SEXP _modsem_dmvnrm_arma_mc(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec const& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool const >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_arma_mc(x, mean, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}
// totalDmvnWeightedCpp
double totalDmvnWeightedCpp(const arma::vec& mu, const arma::mat& sigma, const arma::vec& nu, const arma::mat& S, double tgamma, int n, int d);
RcppExport SEXP _modsem_totalDmvnWeightedCpp(SEXP muSEXP, SEXP sigmaSEXP, SEXP nuSEXP, SEXP SSEXP, SEXP tgammaSEXP, SEXP nSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type tgamma(tgammaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(totalDmvnWeightedCpp(mu, sigma, nu, S, tgamma, n, d));
    return rcpp_result_gen;
END_RCPP
}
// tracePathsNumericCpp
Rcpp::NumericVector tracePathsNumericCpp(Rcpp::CharacterVector x, Rcpp::CharacterVector y, Rcpp::DataFrame parTable, int maxlen);
RcppExport SEXP _modsem_tracePathsNumericCpp(SEXP xSEXP, SEXP ySEXP, SEXP parTableSEXP, SEXP maxlenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type parTable(parTableSEXP);
    Rcpp::traits::input_parameter< int >::type maxlen(maxlenSEXP);
    rcpp_result_gen = Rcpp::wrap(tracePathsNumericCpp(x, y, parTable, maxlen));
    return rcpp_result_gen;
END_RCPP
}
// tracePathsCharacterCpp
Rcpp::List tracePathsCharacterCpp(Rcpp::CharacterVector x, Rcpp::CharacterVector y, Rcpp::DataFrame parTable, std::string paramCol, int maxlen);
RcppExport SEXP _modsem_tracePathsCharacterCpp(SEXP xSEXP, SEXP ySEXP, SEXP parTableSEXP, SEXP paramColSEXP, SEXP maxlenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type parTable(parTableSEXP);
    Rcpp::traits::input_parameter< std::string >::type paramCol(paramColSEXP);
    Rcpp::traits::input_parameter< int >::type maxlen(maxlenSEXP);
    rcpp_result_gen = Rcpp::wrap(tracePathsCharacterCpp(x, y, parTable, paramCol, maxlen));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_modsem_calcSESimpleSlopes", (DL_FUNC) &_modsem_calcSESimpleSlopes, 2},
    {"_modsem_muLmsCpp", (DL_FUNC) &_modsem_muLmsCpp, 2},
    {"_modsem_sigmaLmsCpp", (DL_FUNC) &_modsem_sigmaLmsCpp, 2},
    {"_modsem_completeLogLikLmsCpp", (DL_FUNC) &_modsem_completeLogLikLmsCpp, 3},
    {"_modsem_gradLogLikLmsCpp", (DL_FUNC) &_modsem_gradLogLikLmsCpp, 6},
    {"_modsem_muQmlCpp", (DL_FUNC) &_modsem_muQmlCpp, 2},
    {"_modsem_sigmaQmlCpp", (DL_FUNC) &_modsem_sigmaQmlCpp, 2},
    {"_modsem_calcKronXi", (DL_FUNC) &_modsem_calcKronXi, 2},
    {"_modsem_calcBinvCpp", (DL_FUNC) &_modsem_calcBinvCpp, 2},
    {"_modsem_dnormCpp", (DL_FUNC) &_modsem_dnormCpp, 3},
    {"_modsem_varZCpp", (DL_FUNC) &_modsem_varZCpp, 3},
    {"_modsem_multiplyIndicatorsCpp", (DL_FUNC) &_modsem_multiplyIndicatorsCpp, 1},
    {"_modsem_rep_dmvnorm", (DL_FUNC) &_modsem_rep_dmvnorm, 4},
    {"_modsem_dmvnrm_arma_mc", (DL_FUNC) &_modsem_dmvnrm_arma_mc, 4},
    {"_modsem_totalDmvnWeightedCpp", (DL_FUNC) &_modsem_totalDmvnWeightedCpp, 7},
    {"_modsem_tracePathsNumericCpp", (DL_FUNC) &_modsem_tracePathsNumericCpp, 4},
    {"_modsem_tracePathsCharacterCpp", (DL_FUNC) &_modsem_tracePathsCharacterCpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_modsem(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
