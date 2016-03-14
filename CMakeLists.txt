cmake_minimum_required(VERSION 2.8)
include(CheckCXXCompilerFlag)

add_definitions("-g")
add_definitions("-DNDEBUG")

# link_directories(~/gitkawailab/cluster-expand/spglib/src/.libs)

# set(CMAKE_CXX_COMPILER icpc)
# set(CMAKE_CXX_COMPILER_ID Intel)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  add_definitions("-fast -DEIGEN_NO_DEBUG")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()
message("${CMAKE_CXX_COMPILER_ID} Compiler ")

CHECK_CXX_COMPILER_FLAG("-std=c++11"  COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)

add_executable(cluster
	./src/parsePoscar.cpp
	./src/cluster.cpp
)

add_executable(wang-landau
	./src/wang-landau.cpp
	./src/parser.cpp
	./src/parsePoscar.cpp
	./src/conf2corr.cpp
)
