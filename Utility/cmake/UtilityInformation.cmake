function(get_output TO_EXECUTE RESULT_VAR)
    execute_process(
        COMMAND ${TO_EXECUTE}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE ${RESULT_VAR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endfunction()

set(INFO_HEADER ${CMAKE_BINARY_DIR}/utility_info.h)

get_output("git rev-parse HEAD" GIT_COMMIT_HASH)
get_output("git log -1 --format=%B" GIT_COMMIT_NAME)
get_output("git log -1 --format=%cd" GIT_COMMIT_DATE)
get_output("date" COMPILE_DATE)
message("${GIT_COMMIT_HASH}")

configure_file(
    ${CMAKE_SOURCE_DIR}/utility_info.h.in
    ${INFO_HEADER}
    @ONLY
)