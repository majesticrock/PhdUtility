cmake_minimum_required(VERSION 3.16)
project(iEoM_tests)

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

add_executable(block_matrix_test block_matrix.cpp)
SET_COMPILER_FLAGS(block_matrix_test)
target_link_libraries(block_matrix_test PRIVATE mrock::iEoM)

add_executable(ieom_bcs_test ieom_bcs.cpp)
SET_COMPILER_FLAGS(ieom_bcs_test)
target_link_libraries(ieom_bcs_test PRIVATE mrock::iEoM)

# Enable CTest
enable_testing()

# Add a test
add_test(NAME block_matrix_test COMMAND block_matrix_test)
add_test(NAME ieom_bcs_test COMMAND ieom_bcs_test)