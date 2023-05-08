cmake_minimum_required(VERSION 3.14)
project(BPCells)

find_package(HDF5 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

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

add_library(
    arrayIO
    ${SRC}/arrayIO/array_interfaces.cpp
    ${SRC}/arrayIO/binaryfile.cpp
    ${SRC}/arrayIO/hdf5.cpp
    ${SRC}/arrayIO/vector.cpp
    ${SRC}/arrayIO/bp128.cpp
)

target_link_libraries(
    arrayIO 
    ${HDF5_LIBRARIES}
    bitpacking
)
target_include_directories(
    arrayIO 
    PUBLIC 
    ${HDF5_INCLUDE_DIRS}
    ${SRC}
)

add_library(
    fragmentIterators
    ${SRC}/fragmentIterators/BedFragments.cpp
    ${SRC}/fragmentIterators/CellSelect.cpp
    ${SRC}/fragmentIterators/ChrSelect.cpp
    ${SRC}/fragmentIterators/FragmentIterator.cpp
    ${SRC}/fragmentIterators/LengthSelect.cpp
    ${SRC}/fragmentIterators/MergeFragments.cpp
    ${SRC}/fragmentIterators/RegionSelect.cpp
    ${SRC}/fragmentIterators/Rename.cpp
    ${SRC}/fragmentIterators/ShiftCoords.cpp
    ${SRC}/fragmentIterators/StoredFragments.cpp
    ${SRC}/fragmentUtils/InsertionIterator.cpp
)
target_link_libraries(
    fragmentIterators 
    arrayIO
    ${ZLIB_LIBRARIES}
)
target_include_directories(
    fragmentIterators 
    PUBLIC 
    ${ZLIB_INCLUDE_DIRS}
    ${SRC}
)

### TESTING ###
enable_testing()

# For newer cmake, suprress a cmake warning
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip
)

FetchContent_MakeAvailable(googletest)
include(GoogleTest)


# Build bitpacking tests for all available SIMD modes
if(TARGET_ARM)
    list(APPEND BP128_ARCH_FLAG -march=native -march=armv8-a+nosimd)
    list(APPEND BP128_ARCH_NAME neon fallback)
    list(APPEND BP128_ARCH_ENUM 3 0)
else()
    list(APPEND BP128_ARCH_FLAG -march=native -mno-sse4 -mno-sse2)
    list(APPEND BP128_ARCH_NAME x86_full x86 fallback)
    list(APPEND BP128_ARCH_ENUM 1 2 0)
endif()
list(LENGTH BP128_ARCH_FLAG targets_len)
math(EXPR targets_len "${targets_len} - 1")
foreach(arch_idx RANGE ${targets_len})
    list(GET BP128_ARCH_FLAG ${arch_idx} flag)
    list(GET BP128_ARCH_NAME ${arch_idx} label)
    list(GET BP128_ARCH_ENUM ${arch_idx} enum)
    add_executable(bitpacking_test_${label} test-bitpacking.cpp)
    target_link_libraries(bitpacking_test_${label} bitpacking gtest gmock)
    target_compile_options(bitpacking_test_${label} PUBLIC ${flag})
    gtest_discover_tests(bitpacking_test_${label} EXTRA_ARGS --arch ${enum})
endforeach()

add_executable(
    array_io_test
    test-arrayIO.cpp
)
target_link_libraries(array_io_test arrayIO gtest_main)
gtest_discover_tests(array_io_test)

add_executable(
    fragment_io_test
    test-fragmentIO.cpp
)
target_link_libraries(fragment_io_test fragmentIterators gtest)
gtest_discover_tests(fragment_io_test EXTRA_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../data/mini_fragments.tsv.gz)


add_executable(
    fragment_utils_test
    test-fragmentUtils.cpp
)
target_link_libraries(fragment_utils_test fragmentIterators gtest_main)
gtest_discover_tests(fragment_utils_test)

add_executable(
    matrix_io_test
    test-matrixIO.cpp
    ${SRC}/matrixIterators/ImportMatrixHDF5.cpp
)
target_link_libraries(matrix_io_test arrayIO Eigen3::Eigen gtest_main gmock)
gtest_discover_tests(matrix_io_test)

add_executable(
    matrix_transpose_test
    test-matrixTranspose.cpp
    ${SRC}/matrixIterators/ImportMatrixHDF5.cpp
)
target_link_libraries(matrix_transpose_test arrayIO Eigen3::Eigen gtest_main)
gtest_discover_tests(matrix_transpose_test)

# Build matrix math test for all available SIMD backends
if(TARGET_ARM)
    list(APPEND SLEEF_ARCH_FLAG -march=native -march=armv8-a+nosimd)
    list(APPEND SLEEF_ARCH_NAME NEON FALLBACK)
    list(APPEND SLEEF_ARCH_ENUM 4 0)
else()
    list(APPEND SLEEF_ARCH_FLAG -march=native -march=ivybridge -march=x86-64 -DBPCELLS_SLEEF_FALLBACK) # Warning: got rid of fallback implementation test
    list(APPEND SLEEF_ARCH_NAME AVX2 AVX SSE2 FALLBACK)
    list(APPEND SLEEF_ARCH_ENUM 3 2 1 0)
endif()

list(LENGTH SLEEF_ARCH_FLAG targets_len)
math(EXPR targets_len "${targets_len} - 1")
foreach(arch_idx RANGE ${targets_len})
    list(GET SLEEF_ARCH_FLAG ${arch_idx} flag)
    list(GET SLEEF_ARCH_NAME ${arch_idx} label)
    list(GET SLEEF_ARCH_ENUM ${arch_idx} enum)
    add_executable(matrix_math_test_${label} 
        test-matrixMath.cpp
        ${SRC}/matrixTransforms/MatrixTransform.cpp
        ${SRC}/matrixTransforms/Log1p.cpp
        ${SRC}/matrixTransforms/Pow.cpp
        ${SRC}/matrixTransforms/Min.cpp
        ${SRC}/matrixTransforms/Scale.cpp
        ${SRC}/matrixTransforms/Shift.cpp
    )
    
    target_include_directories(matrix_math_test_${label} PUBLIC ${SRC})
    target_link_libraries(matrix_math_test_${label} Eigen3::Eigen gtest gmock)
    target_compile_options(matrix_math_test_${label} PUBLIC ${flag})
    gtest_discover_tests(matrix_math_test_${label} EXTRA_ARGS --arch ${enum})
endforeach()


# Build sleef simd math test for all available SIMD backends
list(LENGTH SLEEF_ARCH_FLAG targets_len)
math(EXPR targets_len "${targets_len} - 1")
foreach(arch_idx RANGE ${targets_len})
    list(GET SLEEF_ARCH_FLAG ${arch_idx} flag)
    list(GET SLEEF_ARCH_NAME ${arch_idx} label)
    list(GET SLEEF_ARCH_ENUM ${arch_idx} enum)
    add_executable(simd_sleef_test_${label} test-simdSleef.cpp)
    
    target_include_directories(simd_sleef_test_${label} PUBLIC ${SRC})
    target_link_libraries(simd_sleef_test_${label} gtest gmock)
    target_compile_options(simd_sleef_test_${label} PUBLIC ${flag})
    gtest_discover_tests(simd_sleef_test_${label} EXTRA_ARGS --arch ${enum})
endforeach()

add_executable(
    peak_matrix_test
    test-peakMatrix.cpp
    ${SRC}/matrixIterators/PeakMatrix.cpp
    ${SRC}/matrixIterators/TileMatrix.cpp
)
target_link_libraries(peak_matrix_test 
    Eigen3::Eigen
    fragmentIterators 
    gtest_main)
gtest_discover_tests(peak_matrix_test)
