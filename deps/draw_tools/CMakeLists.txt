cmake_minimum_required(VERSION 3.6)
project(Draw_tools C)

option(DRAW_TOOLS_BUILD_EXAMPLES "Build the Draw_tools example programs" ON)

add_library(draw_tools src/draw_tools.c)

set_target_properties(draw_tools PROPERTIES
                      C_STANDARD 99
                      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
target_include_directories(draw_tools PUBLIC include)
target_include_directories(draw_tools PRIVATE src)
target_include_directories(draw_tools PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/shaderIncludes)

add_subdirectory(deps/glfw)
add_subdirectory(deps/glad)

# transform shaders in a valid C string contained in a C header file
file(GLOB_RECURSE shadersGLSL ${CMAKE_CURRENT_SOURCE_DIR}/shaders/*.glsl)
foreach(fileGLSL ${shadersGLSL})
	get_filename_component(basename ${fileGLSL} NAME_WE)
	set(fileH "${CMAKE_CURRENT_BINARY_DIR}/shaderIncludes/${basename}.h")
	list(APPEND shadersH ${fileH})

	file(RELATIVE_PATH printFileH ${CMAKE_SOURCE_DIR} ${fileH})
	file(RELATIVE_PATH printFileGLSL ${CMAKE_SOURCE_DIR} ${fileGLSL})

	add_custom_command(OUTPUT ${fileH}
	                   COMMAND cmake -DIN=${fileGLSL} -DVAR=${basename} -DOUT=${fileH} -P ${CMAKE_CURRENT_SOURCE_DIR}/shaders/shaderToString.cmake
	                   DEPENDS ${fileGLSL}
	                   COMMENT "generating ${printFileH} from shader ${printFileGLSL}")
endforeach(fileGLSL)
add_custom_target(shaderIncludes
                  DEPENDS ${shadersH}
                  SOURCES ${shadersH})
add_dependencies(draw_tools shaderIncludes)

# add dependencies: glfw, glad, -lm
target_link_libraries(draw_tools
                      PUBLIC glad glfw)



if(DRAW_TOOLS_BUILD_EXAMPLES)
	add_executable(animation ${CMAKE_CURRENT_SOURCE_DIR}/examples/animation.c)
	add_executable(character_table ${CMAKE_CURRENT_SOURCE_DIR}/examples/character_table.c)
	add_executable(hello_world ${CMAKE_CURRENT_SOURCE_DIR}/examples/hello_world.c)
	add_executable(koch_snowflake ${CMAKE_CURRENT_SOURCE_DIR}/examples/koch_snowflake.c)
	add_executable(markers ${CMAKE_CURRENT_SOURCE_DIR}/examples/markers.c)
	add_executable(points_and_lines ${CMAKE_CURRENT_SOURCE_DIR}/examples/points_and_lines.c)
	add_executable(text_style ${CMAKE_CURRENT_SOURCE_DIR}/examples/text_style.c)

	target_link_libraries(animation draw_tools)
	target_link_libraries(character_table draw_tools)
	target_link_libraries(hello_world draw_tools)
	target_link_libraries(koch_snowflake draw_tools)
	target_link_libraries(markers draw_tools)
	target_link_libraries(points_and_lines draw_tools)
	target_link_libraries(text_style draw_tools)

	set_target_properties(animation
	                      character_table
	                      koch_snowflake
	                      markers
	                      points_and_lines
	                      text_style
	                      hello_world
	                      PROPERTIES
	                      C_STANDARD 99
	                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/examples")

endif(DRAW_TOOLS_BUILD_EXAMPLES)