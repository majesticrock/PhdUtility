cmake_minimum_required(VERSION 3.16)
project(mrock VERSION 3.0.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CTest)

set(MROCK_VERSION ${PROJECT_VERSION})
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(BUILD_TESTS "Build the unit tests" ON)
option(MROCK_ENABLE_NATIVE_ARCH "Enable architecture-specific optimization flags such as -march=native" OFF)

add_library(mrock INTERFACE)
add_library(mrock::mrock ALIAS mrock)

target_include_directories(mrock INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock_install.cmake)
mrock_install(mrock)

if(BUILD_TESTS)
    set(UTIL_TESTS ON)
    set(SYM_OP_TESTS ON)
    set(IEOM_TESTS ON)
endif()

add_subdirectory(utility)
add_subdirectory(symbolic_operators)
add_subdirectory(iEoM)

include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock-information.cmake)
install(FILES ${MROCK_INFO_HEADER} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mrock)