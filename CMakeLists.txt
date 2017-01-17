cmake_minimum_required(VERSION 2.8)
include(CheckCXXCompilerFlag)
project(wlce)

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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ./lib)

add_library(common STATIC ./src/input.cpp ./src/parser.cpp)
add_library(conf STATIC ./src/conf2corr.cpp ./src/metroconf.cpp ./src/wlconf.cpp)
add_library(cluster-lib STATIC ./src/site.cpp ./src/myindex.cpp ./src/search_equiv_sites.cpp)
#
add_executable(cluster ./src/cluster.cpp)
target_link_libraries(cluster common cluster-lib)

add_executable(getconf ./src/getconf.cpp)
 target_link_libraries(getconf common conf cluster-lib)

add_executable(metropolis ./src/metropolis.cpp)
target_link_libraries(metropolis common conf)

add_executable(wang-landau ./src/wang-landau.cpp)
target_link_libraries(wang-landau common conf)

add_executable(directZ ./src/directZ.cpp ./src/wlconfT.cpp)
target_link_libraries(directZ common conf)


find_package(Boost 1.4 COMPONENTS unit_test_framework)
if(Boost_FOUND)
	add_executable(test-main
		tests/test-main.cpp
		tests/test-conf2corr.cpp
		tests/test-metroconf.cpp
	)
	target_link_libraries(test-main common conf ${Boost_LIBRARIES})
	enable_testing()
	add_test(
		NAME MyTest
		COMMAND $<TARGET_FILE:test-main>
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
	)
endif()
