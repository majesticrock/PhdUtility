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

set(INFO_HEADER ${CMAKE_CURRENT_BINARY_DIR}/info.h)

run_command(MROCK_GIT_COMMIT_VERSION git describe --always --dirty)
run_command(MROCK_GIT_COMMIT_NAME git log -1 --format=%B)
run_command(MROCK_GIT_COMMIT_DATE git log -1 --format=%cd)
run_command(MROCK_MAKE_DATE "date")
run_command(MROCK_HOSTNAME "hostname")

run_command(_MROCK_COMPILER_VERSION ${CMAKE_CXX_COMPILER} --version | sed 1q)
string(REPLACE "\n" ";" _MROCK_COMPILER_VERSION_LINES "${_MROCK_COMPILER_VERSION}")
list(GET _MROCK_COMPILER_VERSION_LINES 0 MROCK_COMPILER_VERSION)
unset(_MROCK_COMPILER_VERSION)
unset(_MROCK_COMPILER_VERSION_LINES)

configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/info.h.in
    ${INFO_HEADER}
    @ONLY
)

include(${CMAKE_CURRENT_LIST_DIR}/mrock-message.cmake)
mrock_message("Generated metadata header file.")