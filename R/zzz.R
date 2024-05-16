
# Run all of the rlang `on_load()` hooks for registering 
# generics of suggested dependencies
.onLoad <- function(...) {
    rlang::run_on_load()
}