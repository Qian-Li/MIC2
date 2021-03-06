// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// MIC_mcmc
List MIC_mcmc(Rcpp::List const& data, int const& K, int const& run, int const& thin);
RcppExport SEXP _MIC2_MIC_mcmc(SEXP dataSEXP, SEXP KSEXP, SEXP runSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List const& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int const& >::type K(KSEXP);
    Rcpp::traits::input_parameter< int const& >::type run(runSEXP);
    Rcpp::traits::input_parameter< int const& >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(MIC_mcmc(data, K, run, thin));
    return rcpp_result_gen;
END_RCPP
}
// specParzen
arma::vec specParzen(arma::vec const& ts, int lag, int const& maxf, int const& outn);
RcppExport SEXP _MIC2_specParzen(SEXP tsSEXP, SEXP lagSEXP, SEXP maxfSEXP, SEXP outnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec const& >::type ts(tsSEXP);
    Rcpp::traits::input_parameter< int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< int const& >::type maxf(maxfSEXP);
    Rcpp::traits::input_parameter< int const& >::type outn(outnSEXP);
    rcpp_result_gen = Rcpp::wrap(specParzen(ts, lag, maxf, outn));
    return rcpp_result_gen;
END_RCPP
}
// SpecSim
arma::cube SpecSim(arma::cube const& ts, int lag, int const& wn, int const& win, int const& overlap, int const& specN);
RcppExport SEXP _MIC2_SpecSim(SEXP tsSEXP, SEXP lagSEXP, SEXP wnSEXP, SEXP winSEXP, SEXP overlapSEXP, SEXP specNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube const& >::type ts(tsSEXP);
    Rcpp::traits::input_parameter< int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< int const& >::type wn(wnSEXP);
    Rcpp::traits::input_parameter< int const& >::type win(winSEXP);
    Rcpp::traits::input_parameter< int const& >::type overlap(overlapSEXP);
    Rcpp::traits::input_parameter< int const& >::type specN(specNSEXP);
    rcpp_result_gen = Rcpp::wrap(SpecSim(ts, lag, wn, win, overlap, specN));
    return rcpp_result_gen;
END_RCPP
}
// EigLap
arma::cube EigLap(arma::cube  const& data, int         const& D, bool        const& normal);
RcppExport SEXP _MIC2_EigLap(SEXP dataSEXP, SEXP DSEXP, SEXP normalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube  const& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int         const& >::type D(DSEXP);
    Rcpp::traits::input_parameter< bool        const& >::type normal(normalSEXP);
    rcpp_result_gen = Rcpp::wrap(EigLap(data, D, normal));
    return rcpp_result_gen;
END_RCPP
}
// SpecOnly
arma::cube SpecOnly(arma::cube const& ts, int lag, int const& wn, int const& win, int const& overlap, int const& specN);
RcppExport SEXP _MIC2_SpecOnly(SEXP tsSEXP, SEXP lagSEXP, SEXP wnSEXP, SEXP winSEXP, SEXP overlapSEXP, SEXP specNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube const& >::type ts(tsSEXP);
    Rcpp::traits::input_parameter< int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< int const& >::type wn(wnSEXP);
    Rcpp::traits::input_parameter< int const& >::type win(winSEXP);
    Rcpp::traits::input_parameter< int const& >::type overlap(overlapSEXP);
    Rcpp::traits::input_parameter< int const& >::type specN(specNSEXP);
    rcpp_result_gen = Rcpp::wrap(SpecOnly(ts, lag, wn, win, overlap, specN));
    return rcpp_result_gen;
END_RCPP
}
// clustalign
void clustalign(arma::mat& now, arma::mat const& ref);
RcppExport SEXP _MIC2_clustalign(SEXP nowSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type now(nowSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type ref(refSEXP);
    clustalign(now, ref);
    return R_NilValue;
END_RCPP
}
