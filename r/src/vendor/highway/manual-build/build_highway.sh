#!/usr/bin/env bash

set -euo pipefail

# This is a fairly direct translation of the google/highway CMakeLists.txt from version 1.0.5
# https://github.com/google/highway/blob/1.0.5/CMakeLists.txt
# We ignore compatibility with MSVC and pre-2018 versions of GCC/clang

# Other omitted option: HWY_CMAKE_ARM7; HWY_RISCV

# Usage: CXX=/path/to/c++/compiler bash build_highway.sh [highway-source] [build-output-dir]
# Respects CXX environment variable
# Output dir will have lib/libhwy.a and include/hwy/*.h

SCRIPT_DIR=$(dirname $0)

SRC_DIR="$1"
OUT_DIR="$2"

# if (NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE RelWithDebInfo)
HWY_FLAGS=(
    -O2 
    -g 
    -DNDEBUG
)

# Skip the CONTRIB, since we don't need sorting or image libraries
HWY_SOURCES=(
    hwy/aligned_allocator.cc 
    hwy/nanobenchmark.cc
    hwy/per_target.cc
    hwy/print.cc
    hwy/targets.cc
)

HWY_HEADERS=(
    hwy/contrib/math/math-inl.h
    hwy/contrib/algo/copy-inl.h
    hwy/contrib/algo/find-inl.h
    hwy/contrib/algo/transform-inl.h

    hwy/aligned_allocator.h
    hwy/base.h
    hwy/cache_control.h
    hwy/detect_compiler_arch.h  # private
    hwy/detect_targets.h  # private
    hwy/foreach_target.h
    hwy/highway.h
    hwy/highway_export.h
    hwy/nanobenchmark.h
    hwy/ops/arm_neon-inl.h
    hwy/ops/arm_sve-inl.h
    hwy/ops/emu128-inl.h
    hwy/ops/generic_ops-inl.h
    hwy/ops/ppc_vsx-inl.h
    hwy/ops/rvv-inl.h
    hwy/ops/scalar-inl.h
    hwy/ops/set_macros-inl.h
    hwy/ops/shared-inl.h
    hwy/ops/wasm_128-inl.h
    hwy/ops/tuple-inl.h
    hwy/ops/x86_128-inl.h
    hwy/ops/x86_256-inl.h
    hwy/ops/x86_512-inl.h
    hwy/per_target.h
    hwy/print-inl.h
    hwy/print.h
    hwy/targets.h
    hwy/timer-inl.h
)


# Avoid changing binaries based on the current time and date.
HWY_FLAGS+=(
    -Wno-builtin-macro-redefined
    '-D__DATE__="redacted"'
    '-D__TIMESTAMP__="redacted"'
    '-D__TIME__="redacted"'
)

# Optimizations
HWY_FLAGS+=(
    -fmerge-all-constants
)

# Warnings
HWY_FLAGS+=(
    -Wall 
    -Wextra
)

# These are not included in Wall nor Wextra:
HWY_FLAGS+=(
    -Wconversion
    -Wsign-conversion
    -Wvla
    -Wnon-virtual-dtor
)

# if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
if bash $SCRIPT_DIR/check_macro_defined.sh __clang__; then
    HWY_FLAGS+=(
        -Wfloat-overflow-conversion
        -Wfloat-zero-conversion
        -Wfor-loop-analysis
        -Wgnu-redeclared-enum
        -Winfinite-recursion
        -Wself-assign
        -Wstring-conversion
        -Wtautological-overlap-compare
        -Wthread-safety-analysis
        -Wundefined-func-template

        -fno-cxx-exceptions
        -fno-slp-vectorize
        -fno-vectorize 
    )
    # Use color in messages
    HWY_FLAGS+=(-fdiagnostics-show-option -fcolor-diagnostics)
fi

# if (WIN32)
if bash $SCRIPT_DIR/check_macro_defined.sh __WIN32__; then
    # if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    if bash $SCRIPT_DIR/check_macro_defined.sh __clang__; then
        HWY_FLAGS+=(
            -Wno-global-constructors
            -Wno-language-extension-token
            -Wno-used-but-marked-unused
            -Wno-shadow-field-in-constructor
            -Wno-unused-member-function
            -Wno-unused-template
            -Wno-c++98-compat-pedantic
            -Wno-used-but-marked-unused
            -Wno-zero-as-null-pointer-constant
        )
    fi
    HWY_FLAGS+=(
        -Wno-cast-align
        -Wno-double-promotion
        -Wno-float-equal
        -Wno-format-nonliteral
        -Wno-shadow
        -Wno-sign-conversion
    )
else
    HWY_FLAGS+=(
        -fmath-errno
        -fno-exceptions
    )
fi


# Skip: if (HWY_CMAKE_ARM7)

# if(HWY_RISCV)
if bash $SCRIPT_DIR/check_macro_defined.sh __riscv; then
    HWY_FLAGS+=(-march=rv64gcv1p0)
    # if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    if bash $SCRIPT_DIR/check_macro_defined.sh __clang__; then
        HWY_FLAGS+=(-menable-experimental-extensions)
    fi
fi

# if (HWY_EMSCRIPTEN)
if bash $SCRIPT_DIR/check_macro_defined.sh __EMSCRIPTEN__; then
    HWY_FLAGS+=(-matomics)
fi

# check_include_file(asm/hwcap.h HAVE_ASM_HWCAP_H)
# if (NOT HAVE_ASM_HWCAP_H)
if ! bash $SCRIPT_DIR/check_include_file.sh asm/hwcap.h; then
    HWY_FLAGS+=(-DTOOLCHAIN_MISS_ASM_HWCAP_H)
fi

# if(NOT HAVE_SYS_AUXV_H)
if ! bash $SCRIPT_DIR/check_include_file.sh sys/auxv.h; then
    HWY_FLAGS+=(-DTOOLCHAIN_MISS_SYS_AUXV_H)
fi

# target_compile_definitions(hwy PUBLIC "${DLLEXPORT_TO_DEFINE}")
HWY_FLAGS+=(-DHWY_STATIC_DEFINE)

# set_property(TARGET hwy PROPERTY POSITION_INDEPENDENT_CODE ON)
HWY_FLAGS+=(-fPIC)

# target_compile_features(hwy PUBLIC cxx_std_11)
HWY_FLAGS+=(-std=c++11)

# Build the static library
mkdir -p $OUT_DIR/lib $OUT_DIR/include $OUT_DIR/build_tmp

for f in ${HWY_SOURCES[@]}; do
    # Get a .o path of the source file in the ouptput dir
    OBJ_PATH="$(basename $f)"
    OBJ_PATH="$OUT_DIR/build_tmp/${OBJ_PATH/%.cc/.o}"
    # Compile a .o file
    $CXX -c -o "$OBJ_PATH" "$SRC_DIR/$f" "${HWY_FLAGS[@]}" -I "$SRC_DIR"
done

ar rcs "$OUT_DIR/lib/libhwy.a" "$OUT_DIR"/build_tmp/*.o

rm -r "$OUT_DIR/build_tmp"

# Copy headers into the output
for f in "${HWY_HEADERS[@]}"; do
    mkdir -p "$OUT_DIR/include/$(dirname $f)"
    cp "$SRC_DIR/$f" "$OUT_DIR/include/$f"
done