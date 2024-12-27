cmake_minimum_required(VERSION 3.16)

project(SymbolicOperators VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SYM_OP_DIR mrock/SymbolicOperators)

add_library(SymbolicOperators STATIC
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
add_library(mrock::SymbolicOperators ALIAS SymbolicOperators)

# Compiler flags
if(NOT WIN32)
    set(BoldCyan "\\[1;36m")
    set(ColourReset "\\[m")
endif()
if(CLUSTER_BUILD)
    include(${CMAKE_CURRENT_LIST_DIR}/cmake/ClusterCompilerFlags.cmake)
    message("${BoldCyan}Building for cascadelake!${ColourReset}")
else()
    include(${CMAKE_CURRENT_LIST_DIR}/cmake/DefaultCompilerFlags.cmake)
    message("${BoldCyan}Building for the local machine!${ColourReset}")
endif()
SET_COMPILER_FLAGS(SymbolicOperators)

if(ONLY_SYM_OP)
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../Utility ${CMAKE_BINARY_DIR}/Utility)
endif()
target_link_libraries(SymbolicOperators PUBLIC mrock::Utility)

target_include_directories(SymbolicOperators PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${UTILITY_INCLUDE_DIRS}>
)

# Installation settings
install(TARGETS SymbolicOperators
    EXPORT SymbolicOperatorsTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY include/ DESTINATION include)

# Package Configuration
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/SymbolicOperatorsConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_LIST_DIR}/cmake/SymbolicOperatorsConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/SymbolicOperatorsConfig.cmake"
    INSTALL_DESTINATION lib/cmake/${SYM_OP_DIR}
)

# Export Package
install(EXPORT SymbolicOperatorsTargets
    FILE SymbolicOperatorsTargets.cmake
    NAMESPACE mrock::
    DESTINATION lib/cmake/${SYM_OP_DIR}
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SymbolicOperatorsConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SymbolicOperatorsConfigVersion.cmake"
    DESTINATION lib/cmake/${SYM_OP_DIR}
)

if(NOT CLUSTER_BUILD)
    if(SYM_OP_TESTS)
        add_subdirectory(tests)
    endif()
endif()