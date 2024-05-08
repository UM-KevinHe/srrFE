// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// computeDirectExp
arma::vec computeDirectExp(const arma::vec& gamma_prov, const arma::vec& Z_beta, const int& threads);
RcppExport SEXP _srrFE_computeDirectExp(SEXP gamma_provSEXP, SEXP Z_betaSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type gamma_prov(gamma_provSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Z_beta(Z_betaSEXP);
    Rcpp::traits::input_parameter< const int& >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(computeDirectExp(gamma_prov, Z_beta, threads));
    return rcpp_result_gen;
END_RCPP
}
// logis_fe_prov
List logis_fe_prov(arma::vec& Y, arma::mat& Z, arma::vec& n_prov, arma::vec gamma, arma::vec beta, int backtrack, int max_iter, double bound, double tol, bool message);
RcppExport SEXP _srrFE_logis_fe_prov(SEXP YSEXP, SEXP ZSEXP, SEXP n_provSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP backtrackSEXP, SEXP max_iterSEXP, SEXP boundSEXP, SEXP tolSEXP, SEXP messageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type n_prov(n_provSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type backtrack(backtrackSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type bound(boundSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type message(messageSEXP);
    rcpp_result_gen = Rcpp::wrap(logis_fe_prov(Y, Z, n_prov, gamma, beta, backtrack, max_iter, bound, tol, message));
    return rcpp_result_gen;
END_RCPP
}
// logis_BIN_fe_prov
List logis_BIN_fe_prov(arma::vec& Y, arma::mat& Z, arma::vec& n_prov, arma::vec gamma, arma::vec beta, int parallel, int threads, double tol, int max_iter, double bound, bool message, bool backtrack);
RcppExport SEXP _srrFE_logis_BIN_fe_prov(SEXP YSEXP, SEXP ZSEXP, SEXP n_provSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP parallelSEXP, SEXP threadsSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP boundSEXP, SEXP messageSEXP, SEXP backtrackSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type n_prov(n_provSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type bound(boundSEXP);
    Rcpp::traits::input_parameter< bool >::type message(messageSEXP);
    Rcpp::traits::input_parameter< bool >::type backtrack(backtrackSEXP);
    rcpp_result_gen = Rcpp::wrap(logis_BIN_fe_prov(Y, Z, n_prov, gamma, beta, parallel, threads, tol, max_iter, bound, message, backtrack));
    return rcpp_result_gen;
END_RCPP
}
// logis_firth_prov
List logis_firth_prov(arma::vec& Y, arma::mat& Z, arma::vec& n_prov, arma::vec gamma, arma::vec beta, int n_obs, int m, int parallel, int threads, double tol, int max_iter, double bound, bool message, bool backtrack);
RcppExport SEXP _srrFE_logis_firth_prov(SEXP YSEXP, SEXP ZSEXP, SEXP n_provSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP n_obsSEXP, SEXP mSEXP, SEXP parallelSEXP, SEXP threadsSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP boundSEXP, SEXP messageSEXP, SEXP backtrackSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type n_prov(n_provSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type bound(boundSEXP);
    Rcpp::traits::input_parameter< bool >::type message(messageSEXP);
    Rcpp::traits::input_parameter< bool >::type backtrack(backtrackSEXP);
    rcpp_result_gen = Rcpp::wrap(logis_firth_prov(Y, Z, n_prov, gamma, beta, n_obs, m, parallel, threads, tol, max_iter, bound, message, backtrack));
    return rcpp_result_gen;
END_RCPP
}
// wald_covar
List wald_covar(arma::vec& Y, arma::mat& Z, arma::vec& n_prov, arma::vec& gamma, arma::vec& beta, arma::uvec& indices, double null, double alpha);
RcppExport SEXP _srrFE_wald_covar(SEXP YSEXP, SEXP ZSEXP, SEXP n_provSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP indicesSEXP, SEXP nullSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type n_prov(n_provSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< double >::type null(nullSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(wald_covar(Y, Z, n_prov, gamma, beta, indices, null, alpha));
    return rcpp_result_gen;
END_RCPP
}
// Modified_score
arma::vec Modified_score(arma::vec& Y, arma::mat& Z, arma::vec& n_prov, arma::vec gamma, arma::vec beta, double gamma_null, int m, arma::vec parm, int threads);
RcppExport SEXP _srrFE_Modified_score(SEXP YSEXP, SEXP ZSEXP, SEXP n_provSEXP, SEXP gammaSEXP, SEXP betaSEXP, SEXP gamma_nullSEXP, SEXP mSEXP, SEXP parmSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type n_prov(n_provSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_null(gamma_nullSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type parm(parmSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Modified_score(Y, Z, n_prov, gamma, beta, gamma_null, m, parm, threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_srrFE_computeDirectExp", (DL_FUNC) &_srrFE_computeDirectExp, 3},
    {"_srrFE_logis_fe_prov", (DL_FUNC) &_srrFE_logis_fe_prov, 10},
    {"_srrFE_logis_BIN_fe_prov", (DL_FUNC) &_srrFE_logis_BIN_fe_prov, 12},
    {"_srrFE_logis_firth_prov", (DL_FUNC) &_srrFE_logis_firth_prov, 14},
    {"_srrFE_wald_covar", (DL_FUNC) &_srrFE_wald_covar, 8},
    {"_srrFE_Modified_score", (DL_FUNC) &_srrFE_Modified_score, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_srrFE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
