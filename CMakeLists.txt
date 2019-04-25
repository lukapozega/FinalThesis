cmake_minimum_required (VERSION 3.2)
project(assembly)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#The version number.
set (assembly_VERSION_MAJOR 1)
set (assembly_VERSION_MINOR 0)

set(CMAKE_CXX_STANDARD 11)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/Config.h.in"
  "${PROJECT_BINARY_DIR}/Config.h"
)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_library(repeats_parser ${PROJECT_SOURCE_DIR}/src/repeats_parser.cpp)

add_executable(assembly src/main.cpp)
target_link_libraries(assembly bioparser repeats_parser)