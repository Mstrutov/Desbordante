add_custom_target(
	copy-config ALL
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_SOURCE_DIR}/cfg/logging.conf
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
)
