cmake_minimum_required(VERSION 3.16)
project(UtilityTests)

function(SET_COMPILER_FLAGS TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.0)
            message(FATAL_ERROR "GCC version ${CMAKE_CXX_COMPILER_VERSION} is not supported. GCC 11.0 or newer is required for C++20 support.")
        endif()

        target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -fopenmp -march=native -O3)
    else()
        message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only GCC is supported.")
    endif()
endfunction()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)
include_directories("~/usr/local/include")

find_package(OpenMP REQUIRED)

add_executable(minimization_test minimization.cpp)
SET_COMPILER_FLAGS(minimization_test)

add_executable(block_matrix_test block_matrix.cpp)
SET_COMPILER_FLAGS(block_matrix_test)

target_link_libraries(block_matrix_test PRIVATE OpenMP::OpenMP_CXX)

# Enable CTest
enable_testing()

# Add a test
add_test(NAME minimization_test COMMAND minimization_test)
add_test(NAME block_matrix_test COMMAND block_matrix_test)