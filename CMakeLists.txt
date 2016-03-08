cmake_minimum_required(VERSION 2.8)
include(CheckCXXCompilerFlag)

add_definitions("-g")

CHECK_CXX_COMPILER_FLAG("-std=c++11"  COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)
# 
# add_executable(corr
# 	./src/parsePoscar.cpp
# 	./src/corr.cpp
# )

add_executable(corr-sym
	./src/parsePoscar.cpp
	./src/corr-sym.cpp
)

# link_directories(~/gitkawailab/cluster-expand/spglib/src/.libs)
# target_link_libraries(corr
# 	symspg
# 	m
# )

target_link_libraries(corr-sym
	symspg
	m
)
