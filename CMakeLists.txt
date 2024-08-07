cmake_minimum_required(VERSION 3.6)
project(eos)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-unknown-pragmas -Wextra")

file(GLOB HEADERS *.hpp)
file(GLOB SOURCE source/*.cpp)

add_library(eos ${HEADERS} ${SOURCE})

set_target_properties(eos PROPERTIES
        SUFFIX -${CMAKE_BUILD_TYPE}.a
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

include(../matplotlib-cpp/matplotlibcpp.cmake)

target_link_libraries(eos PRIVATE ../gsl_wrapper/libgsl_wrapper-${CMAKE_BUILD_TYPE}.a)