# DefaultCompilerFlags.cmake

function(SET_COMPILER_FLAGS TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -O3)
        if(MROCK_ENABLE_NATIVE_ARCH)
            target_compile_options(${TARGET} PRIVATE -march=native)
        endif()
    else()
        message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only gcc/Clang are supported.")
    endif()
endfunction()