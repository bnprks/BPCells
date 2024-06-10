// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <pybind11/pybind11.h>

#include <chrono>
#include <future>

namespace BPCells {


// Wrap a function call such that we will check for user interrupts
// The function itself is run in a background thread, and the main thread checks for a
// user interrupt every 100ms.
// The function must take a pointer to as std::atomic<bool> as the last argument, and exit
// early if the bool gets set to true
template <class F, class... Args>
std::invoke_result_t<F, Args..., std::atomic<bool> *>
run_with_py_interrupt_check(F &&f, Args &&...args) {
    std::atomic<bool> interrupt(false);
    auto job =
        std::async(std::launch::async, std::forward<F>(f), std::forward<Args>(args)..., &interrupt);
    while (job.wait_for(std::chrono::milliseconds(100)) == std::future_status::timeout) {
        if (PyErr_CheckSignals() != 0) {
            interrupt = true;
        }
    }
    if (interrupt) {
            throw pybind11::error_already_set();
    }
    return job.get();
}

}