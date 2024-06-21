cmake_minimum_required(VERSION 3.16)

# Set the project name for the Utility library
project(Utility LANGUAGES CXX)

# Add the include directory as an interface library
add_library(Utility INTERFACE)

# Specify include directories for the interface library
target_include_directories(Utility INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Optionally install the header-only library (for use with `find_package`)
install(DIRECTORY include/ DESTINATION include)
install(TARGETS Utility EXPORT UtilityConfig)
install(EXPORT UtilityConfig NAMESPACE Utility:: DESTINATION lib/cmake/Utility)
