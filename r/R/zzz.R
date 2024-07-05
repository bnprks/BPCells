
# Run all of the rlang `on_load()` hooks for registering 
# generics of suggested dependencies
.onLoad <- function(...) {
    # Register the logger to use with highfive so HDF5 errors go to Rcerr instead of std::cerr
    register_highfive_logging_callback()

    rlang::run_on_load()
}