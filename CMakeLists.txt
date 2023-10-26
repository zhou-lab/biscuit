# Minimum CMake version
cmake_minimum_required(VERSION 3.21)

# Define BISCUIT version from biscuit.h
file(READ "src/biscuit.h" VER_FILE)

string(REGEX MATCH "#define BISCUIT_VERSION_MAJOR ([0-9]*)" _ ${VER_FILE})
set(VER_MAJOR ${CMAKE_MATCH_1})

string(REGEX MATCH "#define BISCUIT_VERSION_MINOR ([0-9]*)" _ ${VER_FILE})
set(VER_MINOR ${CMAKE_MATCH_1})

string(REGEX MATCH "#define BISCUIT_VERSION_PATCH ([0-9]*)" _ ${VER_FILE})
set(VER_PATCH ${CMAKE_MATCH_1})

# Setup project
project(biscuit
    VERSION "${VER_MAJOR}.${VER_MINOR}.${VER_PATCH}"
    DESCRIPTION "BISulfite-seq Command line User Interface Toolkit"
    HOMEPAGE_URL "https://github.com/huishenlab/biscuit"
    LANGUAGES C
    )

# Set C standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Standard installation paths
include(GNUInstallDirs)

# Setup directories used
set(LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(SCRIPT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/scripts)

# Default build in release mode
# Can be overridden with: cmake -DCMAKE_BUILD_TYPE=Debug
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# External libraries
include(ExternalProject)

# Utils
ExternalProject_Add(utils
    GIT_REPOSITORY https://github.com/huishenlab/utils
    GIT_TAG 5f4aeabb42229440a99a5286f8052d3e5dc665cb
    DOWNLOAD_DIR ${LIB_DIR}
    SOURCE_DIR ${LIB_DIR}/utils
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make libutils.a
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    )

# GSL subset
ExternalProject_Add(sgsl
    GIT_REPOSITORY https://github.com/huishenlab/sgsl
    GIT_TAG 6533277fd490b68b0c664efc168a1f05de54b010
    DOWNLOAD_DIR ${LIB_DIR}
    SOURCE_DIR ${LIB_DIR}/sgsl
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make libgsl.a
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    )

# htslib
find_package(CURL REQUIRED)
set(SHA_CHECK ${SCRIPT_DIR}/confirm_download.sh)
set(htslib_version "1.18")
set(htslib_base "htslib-${htslib_version}")

# NOTE: Will make CRAM files potentially unreadable.
#       If planning to add CRAM functionality in the future, will need to remove --disable-bz2 and --disable-lzma
ExternalProject_Add(htslib
    DOWNLOAD_DIR ${LIB_DIR}
    DOWNLOAD_COMMAND curl --silent --show-error -k -L https://github.com/samtools/htslib/releases/download/${htslib_version}/${htslib_base}.tar.bz2 -o ${htslib_base}.tar.bz2 &&
        ${SHA_CHECK} ${htslib_base}.tar.bz2 f1ab53a593a2320a1bfadf4ef915dae784006c5b5c922c8a8174d7530a9af18f &&
        mkdir -p ${htslib_base} &&
        tar -jxf ${htslib_base}.tar.bz2 --strip-components=1 -C ${htslib_base} &&
        rm -rf htslib &&
        mv -f ${htslib_base} htslib
    SOURCE_DIR ${LIB_DIR}/htslib
    INSTALL_DIR ${LIB_DIR}/install
    CONFIGURE_COMMAND ${LIB_DIR}/htslib/configure --prefix=<INSTALL_DIR> --disable-bz2 --disable-lzma --disable-libcurl
    BUILD_COMMAND make lib-static
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    )

# Include directories
include_directories(${utils_PREFIX} ${sgsl_PREFIX} ${htslib_PREFIX}/htslib)

# Build subcommands
add_subdirectory(src)
include_directories(${LIB_DIR})

# Install
install(TARGETS biscuit DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(PROGRAMS
    ${SCRIPT_DIR}/QC.sh
    ${SCRIPT_DIR}/build_biscuit_QC_assets.pl
    ${SCRIPT_DIR}/flip_pbat_strands.sh
    TYPE BIN
    )
