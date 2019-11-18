cmake_minimum_required(VERSION 3.9)
project(NGP C)

file(GLOB_RECURSE sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
add_executable(project ${sources})

target_include_directories(project PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set_target_properties(project PROPERTIES
                      C_STANDARD 99
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

# add dependency to draw_tools and the libm (math library)
# set(DRAW_TOOLS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) #do not build examples
add_subdirectory(deps/draw_tools)
target_link_libraries(project draw_tools)

# set project as the startup project in visual studio
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT project)