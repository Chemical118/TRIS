cmake_minimum_required(VERSION 3.12)

project(TRIS)

set(CMAKE_CXX_STANDARD 17)

find_package(argparse REQUIRED)
find_package(Boost REQUIRED)
find_package(ZLIB REQUIRED)
find_package(TBB REQUIRED)

add_executable(tris
        src/tris.cpp
        src/kmer.cpp
        src/kmer.h
)

add_executable(tris_test
        test/test.cpp
        src/kmer.cpp
        src/kmer.h
)


target_link_libraries(tris PRIVATE ZLIB::ZLIB TBB::tbb Boost::boost argparse::argparse)
target_link_libraries(tris_test PRIVATE ZLIB::ZLIB TBB::tbb Boost::boost)