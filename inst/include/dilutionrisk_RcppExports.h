// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_dilutionrisk_RCPPEXPORTS_H_GEN_
#define RCPP_dilutionrisk_RCPPEXPORTS_H_GEN_

#include <Rcpp.h>

namespace dilutionrisk {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("dilutionrisk", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("dilutionrisk", "_dilutionrisk_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in dilutionrisk");
            }
        }
    }

    inline NumericVector cpp_rtpois(const int& n, const NumericVector& lambda, const NumericVector& lower, const NumericVector& upper) {
        typedef SEXP(*Ptr_cpp_rtpois)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_cpp_rtpois p_cpp_rtpois = NULL;
        if (p_cpp_rtpois == NULL) {
            validateSignature("NumericVector(*cpp_rtpois)(const int&,const NumericVector&,const NumericVector&,const NumericVector&)");
            p_cpp_rtpois = (Ptr_cpp_rtpois)R_GetCCallable("dilutionrisk", "_dilutionrisk_cpp_rtpois");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cpp_rtpois(Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(lower)), Shield<SEXP>(Rcpp::wrap(upper)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

}

#endif // RCPP_dilutionrisk_RCPPEXPORTS_H_GEN_
