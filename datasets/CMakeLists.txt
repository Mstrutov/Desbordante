# unzip heavy datasets - TODO: not a cross-platform command, works only for Linux
add_custom_target( unZip ALL ) 
add_custom_command(
	TARGET unZip
	COMMAND ${CMAKE_COMMAND} -E tar xzf
		${CMAKE_SOURCE_DIR}/datasets/datasets.zip
	DEPENDS 
		${CMAKE_SOURCE_DIR}/datasets/datasets.zip
	WORKING_DIRECTORY 
		${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/inputData/
	COMMENT "Unpacking datasets.zip"
	VERBATIM
)
