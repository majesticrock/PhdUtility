cmake_minimum_required(VERSION 3.16)
project(symbolic_operatorsTests)

function(SET_COMPILER_FLAGS TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -O3)
        if(MROCK_ENABLE_NATIVE_ARCH)
            target_compile_options(${TARGET} PRIVATE -march=native)
        endif()
    else()
        message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only GCC/Clang are supported.")
    endif()
endfunction()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(continuum_test continuum.cpp)
SET_COMPILER_FLAGS(continuum_test)
target_link_libraries(continuum_test PRIVATE mrock::symbolic_operators)

add_executable(bosons_test bosons.cpp)
SET_COMPILER_FLAGS(bosons_test)
target_link_libraries(bosons_test PRIVATE mrock::symbolic_operators)

add_executable(bcs_test bcs.cpp)
SET_COMPILER_FLAGS(bcs_test)
target_link_libraries(bcs_test PRIVATE mrock::symbolic_operators)

add_executable(string_helper_test string_helper.cpp)
SET_COMPILER_FLAGS(string_helper_test)

enable_testing()

add_test(NAME continuum_test COMMAND continuum_test)
add_test(NAME bosons_test COMMAND bosons_test)
add_test(NAME bcs_test COMMAND bcs_test)
add_test(NAME string_helper_test COMMAND string_helper_test)