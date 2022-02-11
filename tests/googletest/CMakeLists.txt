cmake_minimum_required(VERSION 3.14)
project(my_project)

option(TARGET_ARM "Whether to build SIMD tests for ARM architectures (rather than x86)")

set(CMAKE_CXX_STANDARD 17)

set(SRC ../../src)
add_library(
    bitpacking 
    ${SRC}/bitpacking/bp128.cpp  
    ${SRC}/bitpacking/simd_vec.cpp
)
target_include_directories(
    bitpacking
    PUBLIC
    ${SRC}
)

### TESTING ###
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip
)

FetchContent_MakeAvailable(googletest)
include(GoogleTest)
enable_testing()


# Build bitpacking tests for all available SIMD modes
if(TARGET_ARM)
    add_executable(
        bitpacking_test_fallback
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_fallback bitpacking gtest)
    target_compile_options(bitpacking_test_fallback PUBLIC -march=armv8-a+nofp)
    gtest_discover_tests(bitpacking_test_fallback EXTRA_ARGS --arch 0)

    add_executable(
        bitpacking_test_neon
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_neon bitpacking gtest)
    target_compile_options(bitpacking_test_neon PUBLIC -march=native)
    gtest_discover_tests(bitpacking_test_neon EXTRA_ARGS --arch 3)
else()
    add_executable(
        bitpacking_test_x86_full
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_x86_full bitpacking gtest)
    target_compile_options(bitpacking_test_x86_full PUBLIC -march=native)
    gtest_discover_tests(bitpacking_test_x86_full EXTRA_ARGS --arch 1)

    add_executable(
        bitpacking_test_x86
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_x86 bitpacking gtest)
    target_compile_options(bitpacking_test_x86 PUBLIC -mno-sse4)
    gtest_discover_tests(bitpacking_test_x86 EXTRA_ARGS --arch 2)

    add_executable(
        bitpacking_test_fallback
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_fallback bitpacking gtest)
    target_compile_options(bitpacking_test_fallback PUBLIC -mno-sse2)
    gtest_discover_tests(bitpacking_test_fallback EXTRA_ARGS --arch 0)
endif()
