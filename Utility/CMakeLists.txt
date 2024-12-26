cmake_minimum_required(VERSION 3.16)

project(Utility LANGUAGES CXX VERSION 1.0)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(UTIL_DIR mrock/Utility)

# --- Dependencies ---
include(CMakeFindDependencyMacro)
find_dependency(ZLIB REQUIRED)
find_dependency(Boost REQUIRED COMPONENTS serialization iostreams)

# --- Library Definition ---
add_library(Utility INTERFACE)
add_library(mrock::Utility ALIAS Utility)

target_include_directories(Utility INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(Utility INTERFACE ZLIB::ZLIB ${Boost_LIBRARIES})

include(${CMAKE_CURRENT_LIST_DIR}/cmake/UtilityInformation.cmake)

# --- Installation & Export ---
install(TARGETS Utility
    EXPORT UtilityTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(EXPORT UtilityTargets
    FILE UtilityTargets.cmake
    NAMESPACE mrock::
    DESTINATION lib/cmake/${UTIL_DIR}
)

# --- Config and Version Files ---
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/UtilityConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_LIST_DIR}/cmake/UtilityConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/UtilityConfig.cmake"
    INSTALL_DESTINATION lib/cmake/${UTIL_DIR}
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/UtilityConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/UtilityConfigVersion.cmake"
    DESTINATION lib/cmake/${UTIL_DIR}
)

install(DIRECTORY include/ DESTINATION include)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/utility_info.h" DESTINATION include/${UTIL_DIR})

if(UTIL_TESTS)
    add_subdirectory(tests)
endif()