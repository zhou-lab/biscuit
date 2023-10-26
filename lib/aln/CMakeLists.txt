# Minimum CMake version
cmake_minimum_required(VERSION 3.21)

# Source and header files
file(GLOB sources *.c)
file(GLOB headers *.h)

# Add directories for included libraries
include_directories(../utils)

# Setup static library for index and align
add_library(aln_static STATIC ${sources} ${headers})
set_target_properties(aln_static PROPERTIES OUTPUT_NAME "aln")
target_include_directories(aln_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
