cmake_minimum_required(VERSION 3.16)
project(PhdUtility LANGUAGES CXX)

# Set the C++ standard globally
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Global Installation Paths
#set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")
#set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)

add_subdirectory(Utility)
add_subdirectory(SymbolicOperators)

# Write Global Config File
#include(CMakePackageConfigHelpers)
#configure_package_config_file(
#    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/mrockConfig.cmake.in
#    ${CMAKE_BINARY_DIR}/mrockConfig.cmake
#    INSTALL_DESTINATION lib/cmake/mrock
#)
#
#write_basic_package_version_file(
#    "${CMAKE_BINARY_DIR}/mrockConfigVersion.cmake"
#    VERSION 1.0.0
#    COMPATIBILITY SameMajorVersion
#)
#
## Install Global Configs
#install(FILES
#    "${CMAKE_BINARY_DIR}/mrockConfig.cmake"
#    "${CMAKE_BINARY_DIR}/mrockConfigVersion.cmake"
#    DESTINATION lib/cmake/mrock
#)