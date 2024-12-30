cmake_minimum_required(VERSION 3.16)
project(mrock LANGUAGES CXX VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(mrock INTERFACE)

target_include_directories(mrock INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock_install.cmake)
mrock_install(mrock)

add_subdirectory(utility)
add_subdirectory(symbolic_operators)