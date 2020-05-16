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
add_executable(metanome_cpp ${metanome_cpp_src} src/core/VerticalInfo.h src/core/DependencyCandidate.h src/util/ConfidenceInterval.h src/core/DependencyCandidate.cpp src/util/AgreeSetSample.h src/util/AgreeSetSample.cpp src/util/ListAgreeSetSample.h src/util/test.h src/util/ListAgreeSetSample.cpp src/util/AgreeSetSample_impl.h src/core/DependencyStrategy.h src/core/Configuration.h src/util/VerticalMap.h src/util/VerticalMap.cpp src/core/ProfilingContext.h src/util/PLICache.h src/caching/CachingMethod.h src/caching/CacheEvictionMethod.h src/core/DependencyConsumer.h src/model/PartialFD.h src/model/PartialKey.h src/core/DependencyConsumer.cpp src/util/PLICache.cpp)
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" gtest)
#target_link_libraries(metanome_cpp PRIVATE "stdc++fs" gmock)
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" pthread)
target_link_libraries(metanome_cpp PRIVATE "stdc++fs" Boost::container)