cmake_minimum_required(VERSION 3.10)
project(Desbordante)

# By default select Debug build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

# compiler and platform-dependent settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/target")
if (MSVC)
    add_compile_options(/MT /MTd /EHsc)
    add_compile_options("$<$<CONFIG:Release>:/O2>")

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  "${CMAKE_BINARY_DIR}/target")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  "${CMAKE_BINARY_DIR}/target")
else()
    # -DELPP_THREAD_SAFE -- for easylogging++ thread safety
    add_compile_options("$<$<CONFIG:Debug>:-O0;-DELPP_THREAD_SAFE;-g;-Wall;-Wextra;-Werror;-fno-omit-frame-pointer;-fsanitize=address>")
    add_link_options("$<$<CONFIG:Debug>:-fsanitize=address>")
    add_compile_options("$<$<CONFIG:Release>:-O3;-DELPP_THREAD_SAFE>")
endif()

# configuring boost
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.72.0 REQUIRED COMPONENTS container program_options thread)
include_directories(${Boost_INCLUDE_DIRS})
message(${Boost_INCLUDE_DIRS})

# providing subdirectories for header inclusion
include_directories(
    "src/algorithms"
    "src/algorithms/DFD"
    "src/algorithms/association-rules"
    "src/algorithms/depminer"
    "src/algorithms/Metric_FD"
    "src/model"
    "src/model/types"
    "src/util"
    "src/parser"
    "src/parser/json"
    "src/tests"
    "src/core"
    "src/caching"
    "src/logging"
    "src"

    "lib/easyloggingpp/src"
    "lib/better-enums/"
    )

# adding submodules
add_subdirectory("lib/googletest")

option(build_static_lib "Build easyloggingpp as a static library" ON)
add_subdirectory("lib/easyloggingpp")

add_subdirectory("src")
add_subdirectory("tests")

add_subdirectory("datasets")
add_subdirectory("cfg")
