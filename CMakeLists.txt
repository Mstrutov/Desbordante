# Only for code highlight in CLion

cmake_minimum_required(VERSION 3.10)
project(metanome_cpp)
enable_testing()
set(CMAKE_CXX_STANDARD 17)

SET(BOOST_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "../Software/boost_1_70_0/boost/include")
SET(BOOST_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "../Software/boost_1_70_0/boost/lib")
find_package(Boost 1.70.0 COMPONENTS filesystem regex container)
find_package(GTest REQUIRED)
#find_package(GMock REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(src)
include_directories(/usr/include)
link_directories(/usr/lib)

#if(Boost_FOUND)
#include_directories(container)
    #add_executable(progname file1.cxx file2.cxx)
    #target_link_libraries(metanome_cpp ${Boost_LIBRARIES})
#ADD_DEFINITIONS( "-DHAS_BOOST" )
#endif()

file(GLOB_RECURSE metanome_cpp_src "src/*.h*" "src/*.cpp*")

add_executable(metanome_cpp ${metanome_cpp_src})
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" gtest)
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" gmock)
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" pthread)
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" Boost::container)