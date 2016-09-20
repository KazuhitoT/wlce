cmake_minimum_required(VERSION 2.8)
include(CheckCXXCompilerFlag)
project(wlce)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  add_definitions("-fast -DEIGEN_NO_DEBUG")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()
message("-- Compiler ${CMAKE_CXX_COMPILER}")

CHECK_CXX_COMPILER_FLAG("-std=c++11"  COMPILER_SUPPORTS_CXX11)
if(!COMPILER_SUPPORTS_CXX11)
	message(STATUS "WLCE code requires a compiler which supports C++11.")
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-write-strings")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -s -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g3 -Og -pg")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -s -DNDEBUG")

# seed fix
# add_definitions("-DSEED")

set(CMAKE_CXX_COMPILER mpic++)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)

add_library(mpiparser STATIC ../src/parser.cpp)
add_library(mpiconf STATIC ../src/conf2corr.cpp ../src/metroconf.cpp ../src/wlconf.cpp)

add_executable(remc
	./src/remc.cpp
)
target_link_libraries(remc mpiparser mpiconf)

set(CMAKE_CXX_COMPILER mpic++)
add_executable(remc-gs-search
	./src/remc-gs-search.cpp
)
target_link_libraries(remc-gs-search mpiparser mpiconf)
