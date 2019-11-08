cmake_minimum_required(VERSION 3.9)
project(Project C)

# TODO: add a for each loop for every file in a test directory
add_executable(exec src/main.c)
set_target_properties(exec PROPERTIES C_STANDARD 99)

# add dependency to draw_tools and the libm (math library)
add_subdirectory(deps/draw_tools)
target_link_libraries(exec draw_tools)

# set exec as the startup project in visual studio
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT exec)