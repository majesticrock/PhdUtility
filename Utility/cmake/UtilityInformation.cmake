function(run_command RESULT_VAR)
    execute_process(
        COMMAND ${ARGN}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE ${RESULT_VAR}
        ERROR_VARIABLE ERROR_VAR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE RETURN_CODE
    )
    if(NOT RETURN_CODE EQUAL 0)
        message(WARNING "Error running command: ${ARGN}")
        message(WARNING "Error: ${ERROR_VAR}")
        message(WARNING "Command failed with return code: ${RETURN_CODE}")
    endif()
    set(${RESULT_VAR} "${${RESULT_VAR}}" PARENT_SCOPE)
endfunction()

set(INFO_HEADER ${CMAKE_BINARY_DIR}/utility_info.h)

run_command(GIT_COMMIT_HASH git rev-parse HEAD)
run_command(GIT_COMMIT_NAME git log -1 --format=%B)
run_command(GIT_COMMIT_DATE git log -1 --format=%cd)
run_command(MAKE_DATE "date")

message("GIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
message("GIT_COMMIT_NAME=${GIT_COMMIT_NAME}")
message("GIT_COMMIT_DATE=${GIT_COMMIT_DATE}")
message("MAKE_DATE=${MAKE_DATE}")

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/utility_info.h.in
    ${INFO_HEADER}
    @ONLY
)