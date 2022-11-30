include(FetchContent)

FetchContent_Declare(
        integration
        GIT_REPOSITORY https://github.com/Mstrutov/desbordante-integration.git
        GIT_TAG        f13e712166b2a5af7ea19cf18c077db3cb14e149
        GIT_PROGRESS   true
        DOWNLOAD_DIR   ${CMAKE_SOURCE_DIR}/external/
)

FetchContent_MakeAvailable(integration)

add_custom_target( fetchEclat ALL )
add_custom_command(
        TARGET fetchEclat
        COMMAND ${CMAKE_COMMAND} -E copy
            ${integration_SOURCE_DIR}/borgelt-executable/eclat
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
        DEPENDS
            ${integration_SOURCE_DIR}/external/borgelt-executable
        COMMENT "Copying Borgelt"
        VERBATIM
)

add_custom_target( fetchFpgrowth ALL )
add_custom_command(
        TARGET fetchFpgrowth
        COMMAND ${CMAKE_COMMAND} -E copy
            ${integration_SOURCE_DIR}/borgelt-executable/fpgrowth
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
        DEPENDS
            ${integration_SOURCE_DIR}/external/borgelt-executable
        COMMENT "Copying Borgelt"
        VERBATIM
)
