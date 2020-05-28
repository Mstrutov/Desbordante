cmake_minimum_required(VERSION 3.10)
project(metanome_cpp)
enable_testing()
set(CMAKE_CXX_STANDARD 17)

#SET(BOOST_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "../Software/boost_1_70_0/boost/include")
#SET(BOOST_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "../Software/boost_1_70_0/boost/lib")
SET(CMAKE_CXX_FLAGS "-O0")
SET(CMAKE_C_FLAGS "-O0")
find_package(Boost 1.65.1 COMPONENTS container)
find_package(GTest REQUIRED)
#find_package(GMock REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(src)
include_directories(/usr/include)
link_directories(/usr/lib)

file(GLOB_RECURSE metanome_cpp_src "src/*.h*" "src/*.cpp*")

#надо копировать csvшники
#разные прописать цели здесь
add_executable(metanome_cpp ${metanome_cpp_src})
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" gtest)
#target_link_libraries(metanome_cpp PRIVATE "stdc++fs" gmock)
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" pthread)
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" Boost::container)