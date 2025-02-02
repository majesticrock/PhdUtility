cmake_minimum_required(VERSION 3.16)
project(mrock VERSION 2.2.0 LANGUAGES CXX)

set(MROCK_VERSION ${PROJECT_VERSION})
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(mrock INTERFACE)

target_include_directories(mrock INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock_install.cmake)
mrock_install(mrock)

if(BUILD_TESTS)
    SET(UTIL_TESTS ON)
    set(SYM_OP_TESTS ON)
endif()

add_subdirectory(utility)
add_subdirectory(symbolic_operators)

include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock-information.cmake)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/info.h" DESTINATION include/mrock)

if(BUILD_TESTS)
    include(CTest)
    add_test(
        NAME run_utility_tests
        COMMAND ${CMAKE_CTEST_COMMAND}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/utility/tests
    )
    add_test(
        NAME run_symbolic_operators_tests
        COMMAND ${CMAKE_CTEST_COMMAND}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/symbolic_operators/tests
    )
endif()