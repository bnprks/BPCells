#include <chrono>
#include <future>

#define RCPP_NO_RTTI
#define RCPP_NO_SUGAR
#include <Rcpp.h>

namespace {
inline void myCheckInterruptFn(void * /*dummy*/) { R_CheckUserInterrupt(); }
} // anonymous namespace

// Check for interrupts and throw the sentinel exception if one is pending
inline bool hasUserInterrupt() { return (R_ToplevelExec(myCheckInterruptFn, NULL) == FALSE); }

// Wrap a function call such that we will check for R user interrupts
// The function itself is run in a background thread, and the main thread checks for an
// R user interrupt every 100ms.
// The function must take a pointer to as std::atomic<bool> as the last argument, and exit
// early if the bool gets set to true
// NOTE: It is EXTREMELY IMPORTANT that no R objects are created/destroyed inside the spawned
// thread, which includes destructors that mess with R's GC protection.
template <class F, class... Args>
std::invoke_result_t<F, Args..., std::atomic<bool> *>
run_with_R_interrupt_check(F &&f, Args &&...args) {
    std::atomic<bool> interrupt(false);
    auto job =
        std::async(std::launch::async, std::forward<F>(f), std::forward<Args>(args)..., &interrupt);
    while (job.wait_for(std::chrono::milliseconds(100)) == std::future_status::timeout) {
        if (hasUserInterrupt()) {
            interrupt = true;
        }
    }
    if (interrupt) {
        throw Rcpp::internal::InterruptedException();
    }
    return job.get();
}
