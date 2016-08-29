cmake_minimum_required(VERSION 2.8)
include(CheckCXXCompilerFlag)
SET(BUILD_SHARED_LIBS ON)
project(wlce)

set (WLCE_VERSION_MAJOR 0)
set (WLCE_VERSION_MINOR 7)

# seed fix
# add_definitions("-DSEED")

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
add_definitions("-g -DNDEBUG -Wno-write-strings")

add_executable(cluster
	./src/parser.cpp
	./src/cluster.cpp
)
#
add_executable(getconf
	./src/getconf.cpp
	./src/parser.cpp
	./src/conf2corr.cpp
	./src/metroconf.cpp
)

add_executable(wang-landau
	./src/wang-landau.cpp
	./src/parser.cpp
	./src/conf2corr.cpp
	./src/metroconf.cpp
	./src/wlconf.cpp
	./src/input.cpp
)

add_executable(directZ
	./src/directZ.cpp
	./src/parser.cpp
	./src/conf2corr.cpp
	./src/metroconf.cpp
	./src/wlconfT.cpp
	./src/input.cpp
)

add_executable(metropolis
	./src/metropolis.cpp
	./src/parser.cpp
	./src/conf2corr.cpp
	./src/metroconf.cpp
)

find_package(Boost 1.4 COMPONENTS unit_test_framework-mt)
if(Boost_FOUND)
	add_executable(test-main
		tests/test-main.cpp
		tests/test-conf2corr.cpp
		./src/parser.cpp
		./src/conf2corr.cpp
	)
	target_link_libraries(test-main
		boost_unit_test_framework-mt
	)
	enable_testing()
	add_test(
		NAME MyTest
		COMMAND $<TARGET_FILE:test-main>
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
	)
endif()

set(CMAKE_CXX_COMPILER mpic++)
add_executable(remc
	./src/remc.cpp
	./src/parser.cpp
	./src/conf2corr.cpp
	./src/metroconf.cpp
)

set(CMAKE_CXX_COMPILER mpic++)
add_executable(remc-gs-search
	./src/remc-gs-search.cpp
	./src/parser.cpp
	./src/conf2corr.cpp
	./src/metroconf.cpp
)
