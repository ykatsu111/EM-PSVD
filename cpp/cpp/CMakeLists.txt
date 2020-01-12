cmake_minimum_required(VERSION 3.1)
project(EM-PSVD CXX)

set(CMAKE_BUILD_TYPE "Release")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

if(BUILD_SHARED_LIBS)
  add_library(empsvd SHARED EmpsvdCore.cpp)
else()
  add_library(empsvd STATIC EmpsvdCore.cpp)
endif()
# target_compile_features(empsvd PRIVATE cxx_std_11)
# find_package(Eigen3 REQUIRED)
# target_include_directories(empsvd INTERFACE Eigen3_INCLUDE_DIRS)

add_executable(sample sample.cpp)
target_link_libraries(empsvd)
