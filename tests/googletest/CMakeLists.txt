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

add_library(
    matrixTransforms
    ${SRC}/matrixTransforms/MatrixTransform.cpp
    ${SRC}/matrixTransforms/Log1p.cpp
    ${SRC}/matrixTransforms/Scale.cpp
    ${SRC}/matrixTransforms/Shift.cpp
)
target_link_libraries(
    matrixTransforms
    Eigen3::Eigen
)
target_include_directories(
    matrixTransforms
    PUBLIC
    ${SRC}
)

### TESTING ###
enable_testing()

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip
)

FetchContent_MakeAvailable(googletest)
include(GoogleTest)


# Build bitpacking tests for all available SIMD modes
if(TARGET_ARM)
    add_executable(
        bitpacking_test_fallback
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_fallback bitpacking gtest gmock)
    target_compile_options(bitpacking_test_fallback PUBLIC -march=armv8-a+nofp)
    gtest_discover_tests(bitpacking_test_fallback EXTRA_ARGS --arch 0)

    add_executable(
        bitpacking_test_neon
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_neon bitpacking gtest gmock)
    target_compile_options(bitpacking_test_neon PUBLIC -march=native)
    gtest_discover_tests(bitpacking_test_neon EXTRA_ARGS --arch 3)
else()
    add_executable(
        bitpacking_test_x86_full
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_x86_full bitpacking gtest gmock)
    target_compile_options(bitpacking_test_x86_full PUBLIC -march=native)
    gtest_discover_tests(bitpacking_test_x86_full EXTRA_ARGS --arch 1)

    add_executable(
        bitpacking_test_x86
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_x86 bitpacking gtest gmock)
    target_compile_options(bitpacking_test_x86 PUBLIC -mno-sse4)
    gtest_discover_tests(bitpacking_test_x86 EXTRA_ARGS --arch 2)

    add_executable(
        bitpacking_test_fallback
        test-bitpacking.cpp
    )
    target_link_libraries(bitpacking_test_fallback bitpacking gtest gmock)
    target_compile_options(bitpacking_test_fallback PUBLIC -mno-sse2)
    gtest_discover_tests(bitpacking_test_fallback EXTRA_ARGS --arch 0)
endif()

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
target_link_libraries(fragment_io_test fragmentIterators gtest_main)
gtest_discover_tests(fragment_io_test)


add_executable(
    fragment_utils_test
    test-fragmentUtils.cpp
)
target_link_libraries(fragment_utils_test fragmentIterators gtest_main)
gtest_discover_tests(fragment_utils_test)

add_executable(
    matrix_io_test
    test-matrixIO.cpp
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

add_executable(
    matrix_math_test
    test-matrixMath.cpp
)
target_link_libraries(matrix_math_test matrixTransforms gtest_main gmock)
gtest_discover_tests(matrix_math_test)

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
