cmake_minimum_required(VERSION 3.16)
project(iEoM_tests)

function(SET_COMPILER_FLAGS TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -march=native -O3)
    else()
        message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only GCC is supported.")
    endif()
endfunction()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)
include_directories("~/usr/local/include")

add_executable(block_matrix_test block_matrix.cpp)
SET_COMPILER_FLAGS(block_matrix_test)
target_include_directories(block_matrix_test PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../../iEoM/include)

add_executable(ieom_bcs_test ieom_bcs.cpp)
SET_COMPILER_FLAGS(ieom_bcs_test)
target_include_directories(ieom_bcs_test PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../../iEoM/include)

# Enable CTest
enable_testing()

# Add a test
add_test(NAME block_matrix_test COMMAND block_matrix_test)
add_test(NAME ieom_bcs_test COMMAND ieom_bcs_test)