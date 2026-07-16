# ./CMakeLists.txt
cmake_minimum_required(VERSION 3.30)

project(mrock VERSION 3.0.0 LANGUAGES CXX)
set(MROCK_VERSION ${PROJECT_VERSION})

include(GNUInstallDirs)
include(CTest)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

#
# Default project options
#

set(MROCK_EXTRA_INCLUDE_DIRS
    ""
    CACHE STRING
    "Additional include directories for mrock (semicolon-separated)"
)

set(MROCK_ARCH
    ""
    CACHE STRING
    "Architecture passed to compiler as -march=<arch> (empty disables)"
)

option(MROCK_BUILD_UTILITY
    "Build the utility component"
    ON
)

option(MROCK_BUILD_SYMBOLIC_OPERATORS
    "Build the symbolic operators component"
    ON
)

option(MROCK_BUILD_IEOM
    "Build the iEoM component"
    ON
)

#
# Cmake includes
#

include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock_message.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock_install.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock_build_options.cmake)

#
# Custom user config
#

set(MROCK_USER_CONFIG
    "${CMAKE_SOURCE_DIR}/cmake/mrock-user.cmake"
    CACHE FILEPATH
    "Path to user-specific mrock CMake configuration"
)

if(EXISTS "${MROCK_USER_CONFIG}")
    mrock_message("Loading user configuration: ${MROCK_USER_CONFIG}")
    include("${MROCK_USER_CONFIG}")
else()
    mrock_message("No user configuration at ${MROCK_USER_CONFIG}")
endif()

#
# Core interface target
#

add_library(mrock INTERFACE)
add_library(mrock::mrock ALIAS mrock)

target_compile_features(mrock INTERFACE cxx_std_20)

mrock_message("EXTRA INCLUDES: ${MROCK_EXTRA_INCLUDE_DIRS}")

target_include_directories(mrock INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    ${MROCK_EXTRA_INCLUDE_DIRS}
)

#
# Install core package
#

mrock_install(mrock)
install(
    FILES 
        "${CMAKE_CURRENT_LIST_DIR}/cmake/mrock-information.cmake"
        "${CMAKE_CURRENT_LIST_DIR}/cmake/info.h.in"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mrock"
)

#
# Tests
#

if(BUILD_TESTING)
    set(UTIL_TESTS ON)
    set(SYM_OP_TESTS ON)
    set(IEOM_TESTS ON)
endif()

#
# Components
#

if(MROCK_BUILD_UTILITY)
    add_subdirectory(utility)
endif()

if(MROCK_BUILD_SYMBOLIC_OPERATORS)
    add_subdirectory(symbolic_operators)
endif()

if(MROCK_BUILD_IEOM)
    add_subdirectory(iEoM)
endif()

#
# Version information
#

include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock-information.cmake)
mrock_generate_information_header(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/mrock/info.h"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)

install(
    FILES ${MROCK_INFO_HEADER}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mrock
)