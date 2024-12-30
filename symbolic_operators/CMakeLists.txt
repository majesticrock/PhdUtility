cmake_minimum_required(VERSION 3.16)

project(mrock_symbolic_operators VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(mrock_symbolic_operators STATIC
    sources/Coefficient.cpp
    sources/IndexWrapper.cpp
    sources/Momentum.cpp
    sources/MomentumList.cpp
    sources/Operator.cpp
    sources/OperatorType.cpp
    sources/Term.cpp
    sources/TermLoader.cpp
    sources/Wick.cpp
    sources/WickOperator.cpp
    sources/WickOperatorTemplate.cpp
    sources/WickSymmetry.cpp
    sources/WickTerm.cpp
)
add_library(mrock::symbolic_operators ALIAS mrock_symbolic_operators)

# Compiler flags
if(NOT WIN32)
    string(ASCII 27 Esc)
    set(BoldCyan "${Esc}[1;36m")
    set(ColourReset "${Esc}[m")
endif()
if(CLUSTER_BUILD)
    include(${CMAKE_CURRENT_LIST_DIR}/cmake/ClusterCompilerFlags.cmake)
    message("${BoldCyan}Building for cascadelake!${ColourReset}")
else()
    include(${CMAKE_CURRENT_LIST_DIR}/cmake/DefaultCompilerFlags.cmake)
    message("${BoldCyan}Building for the local machine!${ColourReset}")
endif()
SET_COMPILER_FLAGS(mrock_symbolic_operators)

if(ONLY_SYM_OP)
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../utility ${CMAKE_BINARY_DIR}/utility)
endif()
target_link_libraries(mrock_symbolic_operators PUBLIC mrock::utility)

target_include_directories(mrock_symbolic_operators PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/mrock_install.cmake)
mrock_install(mrock_symbolic_operators)

if(NOT CLUSTER_BUILD)
    if(SYM_OP_TESTS)
        add_subdirectory(tests)
    endif()
endif()