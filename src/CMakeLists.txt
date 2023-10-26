# Minimum CMake version
cmake_minimum_required(VERSION 3.21)

# Source and header files
file(GLOB sources *.c)
file(GLOB headers *.h)

# main.c is what executable is compiled from
list(REMOVE_ITEM sources main.c)
list(REMOVE_ITEM sources test.c)

# Add directories for included libraries
include_directories(../lib/utils)
include_directories(../lib/sgsl)
include_directories(../lib/htslib/htslib)

# Need to create a core library that gets linked with the executable
add_library(biscuit_core ${sources} ${headers})
target_include_directories(biscuit_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Build biscuit executable
add_executable(biscuit main.c)

target_link_libraries(biscuit
    biscuit_core
    m
    ${LIB_DIR}/utils/libutils.a
    ${LIB_DIR}/sgsl/libgsl.a
    ${LIB_DIR}/htslib/libhts.a
    )

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    target_link_libraries(biscuit ${ZLIB_LIBRARIES})
else()
    message(FATAL_ERROR "zlib not found. Required for compilation")
endif(ZLIB_FOUND)

find_package(Threads REQUIRED)
if(Threads_FOUND AND CMAKE_USE_PTHREADS_INIT)
    include_directories(${Threads_INCLUDE_DIRS})
    target_link_libraries(biscuit ${Threads_LIBRARIES})
else()
    message(FATAL_ERROR "pthread not found. Required for compilation")
endif(Threads_FOUND AND CMAKE_USE_PTHREADS_INIT)

find_package(Curses REQUIRED)
if(Curses_FOUND)
    include_directories(${CURSES_INCLUDE_DIRS})
    target_link_libraries(biscuit ${CURSES_LIBRARIES})
else()
    message(FATAL_ERROR "ncurses not found. Required for compilation")
endif(Curses_FOUND)
