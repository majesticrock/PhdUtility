cmake_minimum_required(VERSION 3.16)

project(
    mrock_symbolic_operators
    VERSION ${MROCK_VERSION}
    LANGUAGES CXX
)

include(CMakeFindDependencyMacro)

find_dependency(
    Boost
    REQUIRED
    COMPONENTS
        serialization
        iostreams
)

add_library(
    mrock_symbolic_operators
    STATIC

    sources/Coefficient.cpp
    sources/IndexWrapper.cpp
    sources/Momentum.cpp
    sources/MomentumList.cpp
    sources/Operator.cpp
    sources/OperatorType.cpp
    sources/SumContainer.cpp
    sources/Term.cpp
    sources/TermLoader.cpp
    sources/Wick.cpp
    sources/WickOperator.cpp
    sources/WickOperatorTemplate.cpp
    sources/WickSymmetry.cpp
    sources/WickTerm.cpp
)
add_library(
    mrock::symbolic_operators
    ALIAS
    mrock_symbolic_operators
)
set_target_properties(mrock_symbolic_operators PROPERTIES
    EXPORT_NAME symbolic_operators
)

target_link_libraries(
    mrock_symbolic_operators
    PUBLIC
        mrock::mrock
        Boost::serialization
        Boost::iostreams
)

target_include_directories(
    mrock_symbolic_operators
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
        $<INSTALL_INTERFACE:${MROCK_INSTALL_INCLUDEDIR}>
)

include(
    ${CMAKE_CURRENT_LIST_DIR}/../cmake/mrock_install.cmake
)

mrock_install(mrock_symbolic_operators)

if(SYM_OP_TESTS)
    add_subdirectory(tests)
endif()