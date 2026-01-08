# Highway Build Customization for BPCells

This directory contains scripts and configurations for building the `highway` SIMD library, vendored within the BPCells project.

The original `highway` library was modified to address CRAN check warnings related to `abort()` calls and `stderr` output. These modifications ensure compliance with CRAN policies, which disallow compiled code from terminating R or writing directly to standard error/output streams.

## Modifications:

1.  **`src/vendor/highway/manual-build/build_highway.sh`**:
    *   The `HWY_FLAGS` were updated to include `-DHWY_NO_ABORT`. This preprocessor definition prevents the `highway` library from including `abort()` calls in its compiled output when the relevant code is conditionally compiled.
    *   `hwy/nanobenchmark.cc` and `hwy/print.cc` were removed from `HWY_SOURCES`. These files contained `fprintf(stderr, ...)` and `printf(...)` calls, which are not permitted by CRAN. These components are not essential for BPCells.

2.  **`src/vendor/highway/hwy/targets.cc`**:
    *   The `Abort` function definition was wrapped with `#ifndef HWY_NO_ABORT ... #endif`. This ensures that the `Abort` function (which calls `abort()` and `fprintf(stderr, ...)`) is not compiled into `libhwy.a` when `-DHWY_NO_ABORT` is defined.
    *   The `fprintf(stderr, ...)` warning message within the `DetectTargets` function was commented out to prevent `stderr` output during runtime target detection.

3.  **`src/vendor/highway/hwy/base.h`**:
    *   Added conditional compilation logic (`#ifndef HWY_NO_ABORT`) around the `Abort` function declaration.
    *   Added alternative definitions for `HWY_ABORT` and `HWY_ASSERT` macros when `HWY_NO_ABORT` is defined, making them no-ops. This allows the code to compile without linking to the `Abort` function when the flag is set.

These changes are necessary to pass CRAN checks and should be reapplied if the vendored `highway` version is updated.