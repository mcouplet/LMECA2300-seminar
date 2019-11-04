cmake_minimum_required(VERSION 3.6)
project(Draw_tools C)

add_library(draw_tools src/draw_tools.c)

set_target_properties(draw_tools PROPERTIES C_STANDARD 99)
target_include_directories(draw_tools PUBLIC include)
target_include_directories(draw_tools PRIVATE src)
target_include_directories(draw_tools PRIVATE shaders/include)

add_subdirectory(deps/glfw)
add_subdirectory(deps/glad)

# transform shaders in a valid C string contained in a C header file
file(GLOB_RECURSE shadersGLSL ${CMAKE_CURRENT_SOURCE_DIR}/shaders/src/*.glsl)
foreach(fileGLSL ${shadersGLSL})
    get_filename_component(basename ${fileGLSL} NAME_WE)
    set(fileH "${CMAKE_CURRENT_SOURCE_DIR}/shaders/include/${basename}.h")
    list(APPEND shadersH ${fileH})

    file(RELATIVE_PATH printFileH ${CMAKE_SOURCE_DIR} ${fileH})
    file(RELATIVE_PATH printFileGLSL ${CMAKE_SOURCE_DIR} ${fileGLSL})

    add_custom_command(OUTPUT ${fileH}
                       COMMAND cmake -DIN=${fileGLSL} -DVAR=${basename} -DOUT=${fileH} -P ${CMAKE_CURRENT_SOURCE_DIR}/shaders/shaderToString.cmake
                       DEPENDS ${fileGLSL}
                       COMMENT "generating include file ${printFileH} from shader ${printFileGLSL}")
endforeach(fileGLSL)
add_custom_target(shaderIncludes DEPENDS ${shadersH})
add_dependencies(draw_tools shaderIncludes)

# add dependencies: glfw, glad, -lm
target_link_libraries(draw_tools PUBLIC glad glfw)