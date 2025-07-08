cmake_minimum_required(VERSION 3.16)
project(utilityTests)

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
add_executable(function_time_test function_time.cpp)
SET_COMPILER_FLAGS(function_time_test)
add_executable(block_matrix_test block_matrix.cpp)
SET_COMPILER_FLAGS(block_matrix_test)
add_executable(cauchy_principal_value_test cauchy_principal_value.cpp)
SET_COMPILER_FLAGS(cauchy_principal_value_test)
add_executable(string_utility_test string_utility.cpp)
SET_COMPILER_FLAGS(string_utility_test)
add_executable(trapezoidal_test trapezoidal_test.cpp)
SET_COMPILER_FLAGS(trapezoidal_test)
add_executable(hypergeometric_test hypergeometric.cpp)
SET_COMPILER_FLAGS(hypergeometric_test)

target_link_libraries(block_matrix_test PRIVATE OpenMP::OpenMP_CXX)

# Enable CTest
enable_testing()

# Add a test
add_test(NAME minimization_test COMMAND minimization_test)
add_test(NAME function_time_test COMMAND function_time_test)
add_test(NAME block_matrix_test COMMAND block_matrix_test)
add_test(NAME cauchy_principal_value_test COMMAND cauchy_principal_value_test)
add_test(NAME string_utility_test COMMAND string_utility_test)
add_test(NAME trapezoidal_test COMMAND trapezoidal_test)
add_test(NAME hypergeometric_test COMMAND hypergeometric_test)