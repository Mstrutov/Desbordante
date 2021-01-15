set(BINARY ${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE lib_sources "*/*.h*" "*/*.cpp*" "*/*.cc*")

add_library(${BINARY}_lib STATIC ${lib_sources})

set(run_sources "main.cpp")

add_executable(${BINARY}_run ${run_sources})
target_link_libraries(${BINARY}_run PUBLIC ${CMAKE_PROJECT_NAME}_lib)
