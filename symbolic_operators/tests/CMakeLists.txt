cmake_minimum_required(VERSION 3.16)
project(symbolic_operatorsTests)

function(SET_COMPILER_FLAGS TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.0)
            message(FATAL_ERROR "GCC version ${CMAKE_CXX_COMPILER_VERSION} is not supported. GCC 11.0 or newer is required for C++20 support.")
        endif()

        target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -march=native -O3)
    else()
        message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only GCC is supported.")
    endif()
endfunction()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(continuum_test continuum.cpp)
SET_COMPILER_FLAGS(continuum_test)
target_link_libraries(continuum_test PRIVATE mrock_symbolic_operators)
target_include_directories(continuum_test PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../../symbolic_operators/include)
target_include_directories(continuum_test PRIVATE ${CMAKE_CURRENT_LIST_DIR}/../../utility/include)

enable_testing()

add_test(NAME continuum_test COMMAND continuum_test)