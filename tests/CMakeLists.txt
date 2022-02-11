set(BINARY ${CMAKE_PROJECT_NAME}_test)

# getting GTest
enable_testing()

# building tests
file(GLOB_RECURSE test_sources "*.h*" "*.cpp*")
add_executable(${BINARY} ${test_sources})
add_test(NAME ${BINARY} COMMAND ${BINARY})

# linking with gtest and implemented classes
target_link_libraries(${BINARY} PUBLIC ${CMAKE_PROJECT_NAME}_lib)
target_link_libraries(${BINARY} LINK_PUBLIC gtest gmock Boost::thread easyloggingpp)

# copying sample csv's for testing
add_custom_target(copy-files ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/tests/inputData
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/inputData
        )


