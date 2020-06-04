set(BINARY ${CMAKE_PROJECT_NAME})

file(GLOB_RECURSE sources "*.h*" "*.cpp*")

add_executable(${BINARY}_run ${sources})
add_library(${BINARY}_lib STATIC ${sources})