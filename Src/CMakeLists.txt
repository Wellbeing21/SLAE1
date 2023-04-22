cmake_minimum_required(VERSION 3.24)
project(SLAE)

set(CMAKE_CXX_STANDARD 20)
file(GLOB_RECURSE source source *.h *.hpp *.cpp )
#add_library(TRIDIAGONAL Tridiagonal_matrix/tridiagonal.cpp)
add_library(SLAE ${source})
#target_link_libraries(SLAE TRIDIAGONAL)


set(SRC_FILE Tridiagonal_matrix/tridiagonal.cpp Tridiagonal_matrix/tridiagonal.h solvers/solve_tridiagonal.hpp)


enable_testing()
#target_link_libraries(SLAE gtest gmock pthread)
#target_link_libraries(SLAE gtest_main gmock_main)

#add_subdirectory(gtest)\