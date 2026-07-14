cmake_minimum_required(VERSION 3.16)

project(mrock VERSION 3.0.0 LANGUAGES CXX)
set(MROCK_VERSION ${PROJECT_VERSION})

include(GNUInstallDirs)
include(CTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(MROCK_ENABLE_NATIVE_ARCH
    "Enable architecture-specific optimization flags such as -march=native"
    OFF
)

set(MROCK_EXTRA_INCLUDE_DIRS "" CACHE STRING
    "Additional include directories for mrock"
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
# Core interface target
#
set(MROCK_INSTALL_INCLUDEDIR
    ${CMAKE_INSTALL_INCLUDEDIR}
)

add_library(mrock INTERFACE)
add_library(mrock::mrock ALIAS mrock)

target_include_directories(mrock INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:${MROCK_INSTALL_INCLUDEDIR}>
    ${MROCK_EXTRA_INCLUDE_DIRS}
)

#
# Install core package
#

include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock_install.cmake)

mrock_install(mrock)

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

install(
    FILES ${MROCK_INFO_HEADER}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mrock
)