# DefaultCompilerFlags.cmake

function(SET_COMPILER_FLAGS TARGET)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${TARGET} PRIVATE -Wall -Wno-sign-compare -O3)

        if(NOT DEFINED MROCK_ARCH)
            set(MROCK_ARCH "native")
        endif()
        target_compile_options(${TARGET} PRIVATE -march=${MROCK_ARCH})
        
    else()
        message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only gcc/Clang are supported.")
    endif()
endfunction()