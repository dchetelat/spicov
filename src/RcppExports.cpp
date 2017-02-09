// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// spiCov
NumericMatrix spiCov(arma::mat S, int n, SEXP r = R_NilValue, bool verbose = false);
RcppExport SEXP spicov_spiCov(SEXP SSEXP, SEXP nSEXP, SEXP rSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP );
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        Rcpp::traits::input_parameter< SEXP >::type r(rSEXP );
        Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP );
        NumericMatrix __result = spiCov(S, n, r, verbose);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
