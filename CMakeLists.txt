# ./CMakeLists.txt
cmake_minimum_required(VERSION 3.30)

project(mrock VERSION 3.0.0 LANGUAGES CXX)
set(MROCK_VERSION ${PROJECT_VERSION})

get_filename_component(
    MROCK_DEFAULT_INSTALL_PREFIX
    "${CMAKE_CURRENT_LIST_DIR}/../.mrock"
    ABSOLUTE
)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(
        CMAKE_INSTALL_PREFIX
        "${MROCK_DEFAULT_INSTALL_PREFIX}"
        CACHE PATH
        "Install prefix"
        FORCE
    )
endif()

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

set(MROCK_EXTRA_INCLUDE_DIRS "~/usr/local/include" CACHE STRING
    "Additional include directories for mrock, semicolon-separated")
# Do not use FORCE if the user should be able to override it.

set(MROCK_EXTRA_INCLUDE_DIRS_NORMALIZED "")

foreach(dir IN LISTS MROCK_EXTRA_INCLUDE_DIRS)
    file(REAL_PATH "${dir}" dir_normalized
         BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
         EXPAND_TILDE)

    list(APPEND MROCK_EXTRA_INCLUDE_DIRS_NORMALIZED "${dir_normalized}")
endforeach()

# Assert that a empty prefix path will still work (defaults to ../.mrock/)
if(CMAKE_INSTALL_PREFIX STREQUAL "")
    set(
        CMAKE_INSTALL_PREFIX
        "${MROCK_DEFAULT_INSTALL_PREFIX}"
        CACHE PATH
        "Install prefix"
        FORCE
    )
endif()

mrock_message("Install prefix is '${CMAKE_INSTALL_PREFIX}'.")

#
# Set the default build type as Release
#
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
mrock_message("Build type is '${CMAKE_BUILD_TYPE}'.")

#
# Core interface target
#

add_library(mrock INTERFACE)
add_library(mrock::mrock ALIAS mrock)

target_compile_features(mrock INTERFACE cxx_std_20)

mrock_message("EXTRA INCLUDES: ${MROCK_EXTRA_INCLUDE_DIRS_NORMALIZED}")

target_include_directories(mrock INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    ${MROCK_EXTRA_INCLUDE_DIRS_NORMALIZED}
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