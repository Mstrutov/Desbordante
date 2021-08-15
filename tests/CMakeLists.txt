set(BINARY ${CMAKE_PROJECT_NAME}_test)


# getting GTest
enable_testing()
#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.65.1 REQUIRED COMPONENTS thread)
include_directories(${Boost_INCLUDE_DIRS})

# building tests
file(GLOB_RECURSE test_sources "*.h*" "*.cpp*")
add_executable(${BINARY} ${test_sources})
add_test(NAME ${BINARY} COMMAND ${BINARY})


# linking with gtest and implemented classes
target_link_libraries(${BINARY} PUBLIC gtest gmock ${CMAKE_PROJECT_NAME}_lib ${Boost_LIBRARIES})

# copying sample csv's for testing
add_custom_target(copy-files ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/tests/inputData
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/inputData
        )


