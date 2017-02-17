// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Yik_cpp_arma
double Yik_cpp_arma(arma::mat x, int t, int k);
RcppExport SEXP MultSurvTests_Yik_cpp_arma(SEXP xSEXP, SEXP tSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(Yik_cpp_arma(x, t, k));
    return rcpp_result_gen;
END_RCPP
}
// summandG_cpp_arma
double summandG_cpp_arma(arma::mat x, arma::mat y, int n, int t, double t_delta, int k);
RcppExport SEXP MultSurvTests_summandG_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP nSEXP, SEXP tSEXP, SEXP t_deltaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type t_delta(t_deltaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(summandG_cpp_arma(x, y, n, t, t_delta, k));
    return rcpp_result_gen;
END_RCPP
}
// summandL_cpp_arma
double summandL_cpp_arma(arma::mat x, arma::mat y, int n, int t, double t_delta, int k);
RcppExport SEXP MultSurvTests_summandL_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP nSEXP, SEXP tSEXP, SEXP t_deltaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type t_delta(t_deltaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(summandL_cpp_arma(x, y, n, t, t_delta, k));
    return rcpp_result_gen;
END_RCPP
}
// GehanTest_cpp_arma
double GehanTest_cpp_arma(arma::mat x, arma::mat y, int n1, int n2, arma::mat delta_x, arma::mat delta_y, int k);
RcppExport SEXP MultSurvTests_GehanTest_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(GehanTest_cpp_arma(x, y, n1, n2, delta_x, delta_y, k));
    return rcpp_result_gen;
END_RCPP
}
// mvlogrankTest_cpp_arma
double mvlogrankTest_cpp_arma(arma::mat x, arma::mat y, int n1, int n2, arma::mat delta_x, arma::mat delta_y, int k);
RcppExport SEXP MultSurvTests_mvlogrankTest_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mvlogrankTest_cpp_arma(x, y, n1, n2, delta_x, delta_y, k));
    return rcpp_result_gen;
END_RCPP
}
// muG_cpp_arma
double muG_cpp_arma(arma::mat x, arma::mat y, int t, int k, int n);
RcppExport SEXP MultSurvTests_muG_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP tSEXP, SEXP kSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(muG_cpp_arma(x, y, t, k, n));
    return rcpp_result_gen;
END_RCPP
}
// muL_cpp_arma
double muL_cpp_arma(arma::mat x, arma::mat y, int t, int k);
RcppExport SEXP MultSurvTests_muL_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP tSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(muL_cpp_arma(x, y, t, k));
    return rcpp_result_gen;
END_RCPP
}
// psiG_cpp_arma
double psiG_cpp_arma(arma::mat x, arma::mat y, int t, arma::mat delta_x, int k, int n);
RcppExport SEXP MultSurvTests_psiG_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP tSEXP, SEXP delta_xSEXP, SEXP kSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(psiG_cpp_arma(x, y, t, delta_x, k, n));
    return rcpp_result_gen;
END_RCPP
}
// psiL_cpp_arma
double psiL_cpp_arma(arma::mat x, arma::mat y, int t, arma::mat delta_x, int k);
RcppExport SEXP MultSurvTests_psiL_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP tSEXP, SEXP delta_xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(psiL_cpp_arma(x, y, t, delta_x, k));
    return rcpp_result_gen;
END_RCPP
}
// sigma_iklG_cpp_arma
double sigma_iklG_cpp_arma(arma::mat x, arma::mat y, arma::mat delta_x, arma::mat delta_y, int k, int l, int n);
RcppExport SEXP MultSurvTests_sigma_iklG_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP kSEXP, SEXP lSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_iklG_cpp_arma(x, y, delta_x, delta_y, k, l, n));
    return rcpp_result_gen;
END_RCPP
}
// sigma_iklL_cpp_arma
double sigma_iklL_cpp_arma(arma::mat x, arma::mat y, arma::mat delta_x, arma::mat delta_y, int k, int l);
RcppExport SEXP MultSurvTests_sigma_iklL_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP kSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_iklL_cpp_arma(x, y, delta_x, delta_y, k, l));
    return rcpp_result_gen;
END_RCPP
}
// sigma_klG_cpp_arma
double sigma_klG_cpp_arma(arma::mat x, arma::mat y, arma::mat delta_x, arma::mat delta_y, int k, int l, int n1, int n2);
RcppExport SEXP MultSurvTests_sigma_klG_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP kSEXP, SEXP lSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_klG_cpp_arma(x, y, delta_x, delta_y, k, l, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// sigma_klL_cpp_arma
double sigma_klL_cpp_arma(arma::mat x, arma::mat y, arma::mat delta_x, arma::mat delta_y, int k, int l, int n1, int n2);
RcppExport SEXP MultSurvTests_sigma_klL_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP kSEXP, SEXP lSEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_klL_cpp_arma(x, y, delta_x, delta_y, k, l, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// sigma_G_cpp_arma
arma::mat sigma_G_cpp_arma(arma::mat x, arma::mat y, arma::mat delta_x, arma::mat delta_y, int k, int l, int n1, int n2, int p);
RcppExport SEXP MultSurvTests_sigma_G_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP kSEXP, SEXP lSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_G_cpp_arma(x, y, delta_x, delta_y, k, l, n1, n2, p));
    return rcpp_result_gen;
END_RCPP
}
// sigma_L_cpp_arma
arma::mat sigma_L_cpp_arma(arma::mat x, arma::mat y, arma::mat delta_x, arma::mat delta_y, int k, int l, int n1, int n2, int p);
RcppExport SEXP MultSurvTests_sigma_L_cpp_arma(SEXP xSEXP, SEXP ySEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP kSEXP, SEXP lSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_L_cpp_arma(x, y, delta_x, delta_y, k, l, n1, n2, p));
    return rcpp_result_gen;
END_RCPP
}
// gehan
arma::mat gehan(arma::mat x, arma::mat y, arma::mat delta_x, arma::mat delta_y, int n1, int n2, int p, int k, int l);
RcppExport SEXP MultSurvTests_gehan(SEXP xSEXP, SEXP ySEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP pSEXP, SEXP kSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(gehan(x, y, delta_x, delta_y, n1, n2, p, k, l));
    return rcpp_result_gen;
END_RCPP
}
// mvlogrank
arma::mat mvlogrank(arma::mat x, arma::mat y, arma::mat delta_x, arma::mat delta_y, int n1, int n2, int p, int k, int l);
RcppExport SEXP MultSurvTests_mvlogrank(SEXP xSEXP, SEXP ySEXP, SEXP delta_xSEXP, SEXP delta_ySEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP pSEXP, SEXP kSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type delta_y(delta_ySEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(mvlogrank(x, y, delta_x, delta_y, n1, n2, p, k, l));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP MultSurvTests_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP MultSurvTests_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP MultSurvTests_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP MultSurvTests_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
