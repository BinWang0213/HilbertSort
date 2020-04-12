cmake_minimum_required(VERSION 2.8.12)
project(pyHilbertSort)

find_package(Eigen3)
message(STATUS "EIGEN_DIR: " ${Eigen3_DIR})
set (EIGEN3_INCLUDE_DIR  "${Eigen3_DIR}/../../../include/eigen3")
message(STATUS "EIGEN_INCLUDE: " ${EIGEN3_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(pybind11)
pybind11_add_module(pyHilbertSort src/hilbertSort.cpp)
