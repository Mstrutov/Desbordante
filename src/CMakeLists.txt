set(BINARY ${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE lib_sources "*/*.h*" "*/*.cpp*" "*/*.cc*")

add_library(${BINARY}_lib STATIC ${lib_sources})

set(run_sources "main.cpp")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
add_executable(${BINARY}_run ${run_sources})
target_link_libraries(${BINARY}_run PUBLIC ${CMAKE_PROJECT_NAME}_lib)
target_link_libraries(${BINARY}_run LINK_PUBLIC ${Boost_LIBRARIES} Threads::Threads)
#target_link_libraries(${BINARY}_run )

