# mrock_build_options.cmake
function(mrock_set_build_options target)

    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wextra>
    )

    if(MROCK_ARCH)
        target_compile_options(${target} PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-march=${MROCK_ARCH}>
        )
    endif()

endfunction()