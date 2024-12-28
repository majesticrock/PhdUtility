# ClusterCompilerFlags.cmake

function(SET_COMPILER_FLAGS TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.1)
            message(FATAL_ERROR "GCC version ${CMAKE_CXX_COMPILER_VERSION} is not supported. GCC 12.1 or newer is required for C++20 support.")
        endif()

        target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -fopenmp -march=cascadelake -O3)
    else()
        message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only GCC is supported.")
    endif()
endfunction()