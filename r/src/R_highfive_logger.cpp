#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>

#include <string>
#include <highfive/H5Utility.hpp>

inline void R_highfive_logging_callback(HighFive::LogSeverity severity,
                                     const std::string& message,
                                     const std::string& file,
                                     int line) {
    Rcpp::Rcerr << file << ": " << line << " [" << HighFive::to_string(severity) << "] " << message << std::endl;
}

// [[Rcpp::export]]
void register_highfive_logging_callback() {
    static bool has_registered = false;
    if (!has_registered) {
        HighFive::register_logging_callback(&R_highfive_logging_callback);
        has_registered = true;
    }
}