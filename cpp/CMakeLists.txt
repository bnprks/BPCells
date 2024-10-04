cmake_minimum_required(VERSION 3.14)
project(BPCells)

find_package(HDF5 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(HWY REQUIRED NO_MODULE)

option(TARGET_ARM "Whether to build SIMD tests for ARM architectures (rather than x86)")

set(CMAKE_CXX_STANDARD 17)

# if (NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE RelWithDebInfo)
# endif()

# add_compile_options(-fsanitize=address -fsanitize=undefined)
# add_link_options(-fsanitize=address -fsanitize=undefined)

set(SRC bpcells-cpp)
add_library(
    bitpacking 
    ${SRC}/simd/bp128/d1_maxbits.cpp
    ${SRC}/simd/bp128/d1_pack.cpp
    ${SRC}/simd/bp128/d1_unpack.cpp
    ${SRC}/simd/bp128/d1z_maxbits.cpp
    ${SRC}/simd/bp128/d1z_pack.cpp
    ${SRC}/simd/bp128/d1z_unpack.cpp
    ${SRC}/simd/bp128/diff_maxbits.cpp
    ${SRC}/simd/bp128/diff_pack.cpp
    ${SRC}/simd/bp128/diff_unpack.cpp
    ${SRC}/simd/bp128/for_maxbits.cpp
    ${SRC}/simd/bp128/for_pack.cpp
    ${SRC}/simd/bp128/for_unpack.cpp
    ${SRC}/simd/bp128/vanilla_maxbits.cpp
    ${SRC}/simd/bp128/vanilla_pack.cpp
    ${SRC}/simd/bp128/vanilla_unpack.cpp
)
target_link_libraries(
    bitpacking
    hwy::hwy
)
target_include_directories(
    bitpacking
    PUBLIC
    ${HWY_INCLUDE_DIRS}
    ${SRC}
)

add_library(
    simd_math
    ${SRC}/simd/math.cpp
)
target_link_libraries(
    simd_math
    hwy::hwy
)
target_include_directories(
    simd_math
    PUBLIC
    ${HWY_INCLUDE_DIRS}
    ${SRC}
)

add_library(
    simd_overlaps
    ${SRC}/simd/overlaps.cpp
)
target_link_libraries(
    simd_overlaps
    hwy::hwy
)
target_include_directories(
    simd_overlaps
    PUBLIC
    ${HWY_INCLUDE_DIRS}
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
    simd_math
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
  URL https://github.com/google/googletest/archive/40412d85124f7c6f3d88454583c4633e5e10fc8c.zip
)

FetchContent_MakeAvailable(googletest)
include(GoogleTest)


# Build bitpacking tests for all available SIMD modes
add_executable(
    test-bp128
    tests/test-bp128.cpp
)
target_link_libraries(test-bp128 bitpacking gtest_main gmock)
gtest_discover_tests(test-bp128)

add_executable(
    test-arrayIO
    tests/test-arrayIO.cpp
)
target_link_libraries(test-arrayIO arrayIO gtest_main)
gtest_discover_tests(test-arrayIO)

add_executable(
    test-fragmentIO
    tests/test-fragmentIO.cpp
)
target_link_libraries(test-fragmentIO fragmentIterators gtest)
gtest_discover_tests(test-fragmentIO EXTRA_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../r/tests/data/mini_fragments.tsv.gz)


add_executable(
    test-fragmentUtils
    tests/test-fragmentUtils.cpp
)
target_link_libraries(test-fragmentUtils fragmentIterators gtest_main)
gtest_discover_tests(test-fragmentUtils)

add_executable(
    test-matrixIO
    tests/test-matrixIO.cpp
    ${SRC}/matrixIterators/ImportMatrixHDF5.cpp
)
target_link_libraries(test-matrixIO arrayIO Eigen3::Eigen gtest_main gmock)
gtest_discover_tests(test-matrixIO)

add_executable(
    test-matrixTranspose
    tests/test-matrixTranspose.cpp
    ${SRC}/matrixIterators/ImportMatrixHDF5.cpp
)
target_link_libraries(test-matrixTranspose arrayIO Eigen3::Eigen gtest_main)
gtest_discover_tests(test-matrixTranspose)

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


add_executable(test-matrixMath
    tests/test-matrixMath.cpp
    ${SRC}/matrixTransforms/MatrixTransform.cpp
    ${SRC}/matrixTransforms/Log1p.cpp
    ${SRC}/matrixTransforms/Pow.cpp
    ${SRC}/matrixTransforms/Min.cpp
    ${SRC}/matrixTransforms/Scale.cpp
    ${SRC}/matrixTransforms/Shift.cpp
)

target_include_directories(test-matrixMath PUBLIC ${SRC})
target_link_libraries(test-matrixMath simd_math Eigen3::Eigen gtest_main gmock)
gtest_discover_tests(test-matrixMath)

add_executable(
    test-peakMatrix
    tests/test-peakMatrix.cpp
    ${SRC}/matrixIterators/PeakMatrix.cpp
    ${SRC}/matrixIterators/TileMatrix.cpp
    ${SRC}/simd/overlaps.cpp
)
target_link_libraries(test-peakMatrix 
    Eigen3::Eigen
    fragmentIterators 
    simd_overlaps
    gtest_main)
gtest_discover_tests(test-peakMatrix)
