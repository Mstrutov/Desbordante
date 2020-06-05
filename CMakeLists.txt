cmake_minimum_required(VERSION 3.10)
project(metanome_cpp)

# compiler and platform-dependent settings
set(CMAKE_CXX_STANDARD 17)
if (MSVC)
        set(additional_flags "/MT /MTd /EHsc")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${additional_flags}")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${additional_flags}")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${additional_flags} /O2")

        set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  "${CMAKE_BINARY_DIR}/target" )
        set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  "${CMAKE_BINARY_DIR}/target" )
else()
        set(CMAKE_CXX_FLAGS_DEBUG  "-O0")
        set(CMAKE_CXX_FLAGS_RELEASE  "-O3")

        set( CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/target" )
endif()

# configuring boost
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.65.1 COMPONENTS container)
if (NOT Boost_FOUND)
        message(SEND_ERROR "Failed to find boost::contatiner.")
        return()
else()
        include_directories(${Boost_INCLUDE_DIRS})
endif()

# providing subdirectories for header inclusion
include_directories(
        "scr/algorithms"
        "src/model"
        "src/util"
        "src/parser"
        "src/parser/json"
        "src/tests"
        "src/core"
        "src/caching"
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