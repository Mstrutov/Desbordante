cmake_minimum_required(VERSION 3.10)
project(metanome_cpp)
enable_testing()
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "-O0")
SET(CMAKE_C_FLAGS "-O0")

# configuring boost
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.65.1 COMPONENTS container)
if (NOT Boost_FOUND)
        message(SEND_ERROR "Failed to find boost::contatiner.")
        return()
else()
        include_directories(${Boost_INCLUDE_DIRS})
endif()

# getting GTest
find_package(GTest REQUIRED)
#find_package(GMock REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

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

#надо копировать csvшники
#разные прописать цели здесь
file(GLOB_RECURSE metanome_cpp_src "src/*.h*" "src/*.cpp*")
add_executable(metanome_cpp ${metanome_cpp_src})

target_link_libraries(metanome_cpp PRIVATE "stdc++fs" gtest)
#target_link_libraries(metanome_cpp PRIVATE "stdc++fs" gmock)
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" pthread)
#target_link_libraries(metanome_cpp PRIVATE "stdc++fs" Boost::container)

#install(TARGETS ${PROJECT_NAME} DESTINATION bin)
#install(FILES ${metanome_cpp_src} DESTINATION src)