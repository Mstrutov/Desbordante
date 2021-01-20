cmake_minimum_required(VERSION 3.10)
project(fdtester)

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
    add_compile_options("$<$<CONFIG:Debug>:-Wall;-Wextra;-fno-omit-frame-pointer>")
endif()

# configuring boost
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.65.1 REQUIRED COMPONENTS container program_options)
include_directories(${Boost_INCLUDE_DIRS})

# providing subdirectories for header inclusion
include_directories(
    "src/algorithms"
    "src/model"
    "src/util"
    "src/parser"
    "src/parser/json"
    "src/tests"
    "src/core"
    "src/caching"
    "src/logging"
    "src"
    )

# adding submodules
add_subdirectory("src")
add_subdirectory("lib/googletest")
add_subdirectory("tests")

#target_link_libraries(metanome_cpp PRIVATE "stdc++fs" gmock)
#target_link_libraries(metanome_cpp PRIVATE "stdc++fs" Boost::container)
#install(TARGETS ${PROJECT_NAME} DESTINATION bin)
#install(FILES ${metanome_cpp_src} DESTINATION src)
