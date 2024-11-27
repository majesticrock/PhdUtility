cmake_minimum_required(VERSION 3.16)

project(Utility LANGUAGES CXX)

include(CMakeFindDependencyMacro)
find_dependency(ZLIB REQUIRED)
find_dependency(Boost REQUIRED COMPONENTS serialization iostreams)

add_library(Utility INTERFACE)
target_link_libraries(Utility INTERFACE ZLIB::ZLIB ${Boost_LIBRARIES})

target_include_directories(Utility INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${Boost_INCLUDE_DIRS}
)

add_subdirectory(tests)

install(TARGETS Utility
    EXPORT UtilityTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(EXPORT UtilityTargets
    FILE UtilityTargets.cmake
    NAMESPACE Utility::
    DESTINATION lib/cmake/Utility
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/UtilityConfigVersion.cmake"
    VERSION 1.0
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/UtilityConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/UtilityConfig.cmake"
    INSTALL_DESTINATION lib/cmake/Utility
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/UtilityConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/UtilityConfigVersion.cmake"
    DESTINATION lib/cmake/Utility
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/UtilityConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/UtilityConfig.cmake"
    @ONLY
)

install(DIRECTORY include/ DESTINATION include)
