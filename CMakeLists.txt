cmake_minimum_required(VERSION 3.20)
project(BPCells)

set(CMAKE_CXX_STANDARD 14)

include_directories(src/bitpacking)
include_directories(src/fragmentIterators)
include_directories(src/lib)
include_directories(src/matrixIterators)
include_directories(src/matrixTransforms)


add_executable(BPCells
        src/fragmentIterators/BedFragments.cpp
        src/fragmentIterators/BedFragments.h
        src/fragmentIterators/FragmentsIterator.cpp
        src/fragmentIterators/CellSelect.h
        src/fragmentIterators/CellSelect.cpp
        tests/debug_test.cpp)

target_link_libraries(BPCells z)
target_include_directories(BPCells PRIVATE /Users/ben/Downloads/eigen-3.4.0/Eigen/)
target_compile_definitions(BPCells PRIVATE REAL_EIGEN)