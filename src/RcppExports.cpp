// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// nllh_bvn_censored_ogram
double nllh_bvn_censored_ogram(arma::vec pars, arma::vec x, arma::vec y, arma::vec ux, arma::vec uy, int freq);
RcppExport SEXP _deform_nllh_bvn_censored_ogram(SEXP parsSEXP, SEXP xSEXP, SEXP ySEXP, SEXP uxSEXP, SEXP uySEXP, SEXP freqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ux(uxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type uy(uySEXP);
    Rcpp::traits::input_parameter< int >::type freq(freqSEXP);
    rcpp_result_gen = Rcpp::wrap(nllh_bvn_censored_ogram(pars, x, y, ux, uy, freq));
    return rcpp_result_gen;
END_RCPP
}
// cnorml
double cnorml(arma::vec pars, arma::vec yvec, arma::vec lvec);
RcppExport SEXP _deform_cnorml(SEXP parsSEXP, SEXP yvecSEXP, SEXP lvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yvec(yvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lvec(lvecSEXP);
    rcpp_result_gen = Rcpp::wrap(cnorml(pars, yvec, lvec));
    return rcpp_result_gen;
END_RCPP
}
// cnormr
double cnormr(arma::vec pars, arma::vec yvec, arma::vec rvec);
RcppExport SEXP _deform_cnormr(SEXP parsSEXP, SEXP yvecSEXP, SEXP rvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yvec(yvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rvec(rvecSEXP);
    rcpp_result_gen = Rcpp::wrap(cnormr(pars, yvec, rvec));
    return rcpp_result_gen;
END_RCPP
}
// cnormlr
double cnormlr(arma::vec pars, arma::vec yvec, arma::vec lvec, arma::vec rvec);
RcppExport SEXP _deform_cnormlr(SEXP parsSEXP, SEXP yvecSEXP, SEXP lvecSEXP, SEXP rvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yvec(yvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lvec(lvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rvec(rvecSEXP);
    rcpp_result_gen = Rcpp::wrap(cnormlr(pars, yvec, lvec, rvec));
    return rcpp_result_gen;
END_RCPP
}
// coscovfn
arma::mat coscovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_coscovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(coscovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d1coscovfn
Rcpp::List d1coscovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d1coscovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d1coscovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d2coscovfn
Rcpp::List d2coscovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d2coscovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d2coscovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// covfn
arma::mat covfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_covfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(covfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d1covfn
Rcpp::List d1covfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d1covfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d1covfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d2covfn
Rcpp::List d2covfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d2covfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d2covfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// dampedcoscovfn
arma::mat dampedcoscovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_dampedcoscovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(dampedcoscovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d1dampedcoscovfn
Rcpp::List d1dampedcoscovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d1dampedcoscovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d1dampedcoscovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d2dampedcoscovfn
Rcpp::List d2dampedcoscovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d2dampedcoscovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d2dampedcoscovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// pcovfn
arma::mat pcovfn(Rcpp::List pars, arma::mat z, int n);
RcppExport SEXP _deform_pcovfn(SEXP parsSEXP, SEXP zSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(pcovfn(pars, z, n));
    return rcpp_result_gen;
END_RCPP
}
// pcoscovfn
arma::mat pcoscovfn(arma::vec pars, arma::mat z, int n);
RcppExport SEXP _deform_pcoscovfn(SEXP parsSEXP, SEXP zSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoscovfn(pars, z, n));
    return rcpp_result_gen;
END_RCPP
}
// punitcoscovfn
arma::mat punitcoscovfn(Rcpp::List pars, arma::mat z, int n);
RcppExport SEXP _deform_punitcoscovfn(SEXP parsSEXP, SEXP zSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(punitcoscovfn(pars, z, n));
    return rcpp_result_gen;
END_RCPP
}
// punitcovfn
arma::mat punitcovfn(Rcpp::List pars, arma::mat z, int n);
RcppExport SEXP _deform_punitcovfn(SEXP parsSEXP, SEXP zSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(punitcovfn(pars, z, n));
    return rcpp_result_gen;
END_RCPP
}
// pdampedcoscovfn
arma::mat pdampedcoscovfn(arma::vec pars, arma::mat z, int n);
RcppExport SEXP _deform_pdampedcoscovfn(SEXP parsSEXP, SEXP zSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(pdampedcoscovfn(pars, z, n));
    return rcpp_result_gen;
END_RCPP
}
// tnormr
double tnormr(arma::vec pars, arma::vec yvec, arma::vec rvec);
RcppExport SEXP _deform_tnormr(SEXP parsSEXP, SEXP yvecSEXP, SEXP rvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yvec(yvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rvec(rvecSEXP);
    rcpp_result_gen = Rcpp::wrap(tnormr(pars, yvec, rvec));
    return rcpp_result_gen;
END_RCPP
}
// tnorml
double tnorml(arma::vec pars, arma::vec yvec, arma::vec lvec);
RcppExport SEXP _deform_tnorml(SEXP parsSEXP, SEXP yvecSEXP, SEXP lvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yvec(yvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lvec(lvecSEXP);
    rcpp_result_gen = Rcpp::wrap(tnorml(pars, yvec, lvec));
    return rcpp_result_gen;
END_RCPP
}
// tnormlr
double tnormlr(arma::vec pars, arma::vec yvec, arma::vec lvec, arma::vec rvec);
RcppExport SEXP _deform_tnormlr(SEXP parsSEXP, SEXP yvecSEXP, SEXP lvecSEXP, SEXP rvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yvec(yvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lvec(lvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rvec(rvecSEXP);
    rcpp_result_gen = Rcpp::wrap(tnormlr(pars, yvec, lvec, rvec));
    return rcpp_result_gen;
END_RCPP
}
// unitcoscovfn
arma::mat unitcoscovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_unitcoscovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(unitcoscovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d1unitcoscovfn
Rcpp::List d1unitcoscovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d1unitcoscovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d1unitcoscovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d2unitcoscovfn
Rcpp::List d2unitcoscovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d2unitcoscovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d2unitcoscovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// unitcovfn
arma::mat unitcovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_unitcovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(unitcovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d1unitcovfn
Rcpp::List d1unitcovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d1unitcovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d1unitcovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}
// d2unitcovfn
Rcpp::List d2unitcovfn(Rcpp::List pars, Rcpp::List X, int n);
RcppExport SEXP _deform_d2unitcovfn(SEXP parsSEXP, SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(d2unitcovfn(pars, X, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_deform_nllh_bvn_censored_ogram", (DL_FUNC) &_deform_nllh_bvn_censored_ogram, 6},
    {"_deform_cnorml", (DL_FUNC) &_deform_cnorml, 3},
    {"_deform_cnormr", (DL_FUNC) &_deform_cnormr, 3},
    {"_deform_cnormlr", (DL_FUNC) &_deform_cnormlr, 4},
    {"_deform_coscovfn", (DL_FUNC) &_deform_coscovfn, 3},
    {"_deform_d1coscovfn", (DL_FUNC) &_deform_d1coscovfn, 3},
    {"_deform_d2coscovfn", (DL_FUNC) &_deform_d2coscovfn, 3},
    {"_deform_covfn", (DL_FUNC) &_deform_covfn, 3},
    {"_deform_d1covfn", (DL_FUNC) &_deform_d1covfn, 3},
    {"_deform_d2covfn", (DL_FUNC) &_deform_d2covfn, 3},
    {"_deform_dampedcoscovfn", (DL_FUNC) &_deform_dampedcoscovfn, 3},
    {"_deform_d1dampedcoscovfn", (DL_FUNC) &_deform_d1dampedcoscovfn, 3},
    {"_deform_d2dampedcoscovfn", (DL_FUNC) &_deform_d2dampedcoscovfn, 3},
    {"_deform_pcovfn", (DL_FUNC) &_deform_pcovfn, 3},
    {"_deform_pcoscovfn", (DL_FUNC) &_deform_pcoscovfn, 3},
    {"_deform_punitcoscovfn", (DL_FUNC) &_deform_punitcoscovfn, 3},
    {"_deform_punitcovfn", (DL_FUNC) &_deform_punitcovfn, 3},
    {"_deform_pdampedcoscovfn", (DL_FUNC) &_deform_pdampedcoscovfn, 3},
    {"_deform_tnormr", (DL_FUNC) &_deform_tnormr, 3},
    {"_deform_tnorml", (DL_FUNC) &_deform_tnorml, 3},
    {"_deform_tnormlr", (DL_FUNC) &_deform_tnormlr, 4},
    {"_deform_unitcoscovfn", (DL_FUNC) &_deform_unitcoscovfn, 3},
    {"_deform_d1unitcoscovfn", (DL_FUNC) &_deform_d1unitcoscovfn, 3},
    {"_deform_d2unitcoscovfn", (DL_FUNC) &_deform_d2unitcoscovfn, 3},
    {"_deform_unitcovfn", (DL_FUNC) &_deform_unitcovfn, 3},
    {"_deform_d1unitcovfn", (DL_FUNC) &_deform_d1unitcovfn, 3},
    {"_deform_d2unitcovfn", (DL_FUNC) &_deform_d2unitcovfn, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_deform(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
