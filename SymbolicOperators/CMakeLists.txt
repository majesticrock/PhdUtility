cmake_minimum_required(VERSION 3.16)

# Define the project name and version
project(SymbolicOperators VERSION 1.0.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

set(SYM_OP_SOURCES
    src/Coefficient.cpp
    src/IndexWrapper.cpp
    src/Momentum.cpp
    src/MomentumList.cpp
    src/Operator.cpp
    src/OperatorType.cpp
    src/Term.cpp
    src/TermLoader.cpp
    src/Wick.cpp
    src/WickOperator.cpp
    src/WickOperatorTemplate.cpp
    src/WickSymmetry.cpp
    src/WickTerm.cpp
)

# Add the library target
add_library(SymbolicOperators STATIC ${SYM_OP_SOURCES})

# Specify include directories for the library
target_include_directories(SymbolicOperators PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

find_package(Utility REQUIRED)
if(NOT Utility_FOUND)
    message(FATAL_ERROR "Utility library not found")
endif()
target_link_libraries(SymbolicOperators PRIVATE Utility::Utility)
target_include_directories(SymbolicOperators PRIVATE ${Utility_INCLUDE_DIRS})

# Installation settings
install(TARGETS SymbolicOperators
    EXPORT SymbolicOperatorsTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/SymbolicOperators DESTINATION include)

# Export the package for use with find_package
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/SymbolicOperatorsConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

export(EXPORT SymbolicOperatorsTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/SymbolicOperatorsTargets.cmake"
    NAMESPACE SymbolicOperators::
)

configure_file(SymbolicOperatorsConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/SymbolicOperatorsConfig.cmake"
    @ONLY
)

# Install the export
install(EXPORT SymbolicOperatorsTargets
    FILE SymbolicOperatorsTargets.cmake
    NAMESPACE SymbolicOperators::
    DESTINATION lib/cmake/SymbolicOperators
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SymbolicOperatorsConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SymbolicOperatorsConfigVersion.cmake"
    DESTINATION lib/cmake/SymbolicOperators
)
