# cmake/mrock-information.cmake
include_guard(GLOBAL)

function(_mrock_run_command RESULT_VAR WORKING_DIRECTORY)
    execute_process(
        COMMAND ${ARGN}
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_VARIABLE _mrock_output
        ERROR_VARIABLE _mrock_error
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _mrock_return_code
    )

    if(NOT _mrock_return_code EQUAL 0)
        message(WARNING "Error running command: ${ARGN}")
        message(WARNING "Error: ${_mrock_error}")
        message(WARNING "Command failed with return code: ${_mrock_return_code}")
    endif()

    set(${RESULT_VAR} "${_mrock_output}" PARENT_SCOPE)
endfunction()


function(mrock_generate_information_header)
    set(options)
    set(one_value_args
        OUTPUT
        TEMPLATE
        WORKING_DIRECTORY
        OUT_HEADER
        OUT_INCLUDE_DIR
    )
    set(multi_value_args)

    cmake_parse_arguments(
        MROCK_INFO
        "${options}"
        "${one_value_args}"
        "${multi_value_args}"
        ${ARGN}
    )

    if(NOT MROCK_INFO_WORKING_DIRECTORY)
        set(MROCK_INFO_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
    endif()

    if(NOT MROCK_INFO_OUTPUT)
        set(MROCK_INFO_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/include/mrock/info.h")
    endif()

    if(NOT MROCK_INFO_TEMPLATE)
        set(MROCK_INFO_TEMPLATE "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/info.h.in")
    endif()

    if(NOT EXISTS "${MROCK_INFO_TEMPLATE}")
        message(FATAL_ERROR
            "mrock information header template not found: ${MROCK_INFO_TEMPLATE}"
        )
    endif()

    get_filename_component(_mrock_info_output_dir "${MROCK_INFO_OUTPUT}" DIRECTORY)
    file(MAKE_DIRECTORY "${_mrock_info_output_dir}")

    get_filename_component(_mrock_generated_include_dir "${MROCK_INFO_OUTPUT}" DIRECTORY)
    get_filename_component(_mrock_generated_include_dir "${_mrock_generated_include_dir}" DIRECTORY)

    _mrock_run_command(
        MROCK_GIT_COMMIT_VERSION
        "${MROCK_INFO_WORKING_DIRECTORY}"
        git describe --always --dirty
    )

    _mrock_run_command(
        MROCK_GIT_COMMIT_NAME
        "${MROCK_INFO_WORKING_DIRECTORY}"
        git log -1 --format=%s
    )

    _mrock_run_command(
        MROCK_GIT_COMMIT_DATE
        "${MROCK_INFO_WORKING_DIRECTORY}"
        git log -1 --format=%cd
    )

    string(TIMESTAMP MROCK_MAKE_DATE "%Y-%m-%d %H:%M:%S %z")

    cmake_host_system_information(
        RESULT MROCK_HOSTNAME
        QUERY HOSTNAME
    )

    if(CMAKE_CXX_COMPILER)
        _mrock_run_command(
            _MROCK_COMPILER_VERSION
            "${MROCK_INFO_WORKING_DIRECTORY}"
            "${CMAKE_CXX_COMPILER}" --version
        )

        string(REPLACE "\n" ";" _MROCK_COMPILER_VERSION_LINES "${_MROCK_COMPILER_VERSION}")
        list(GET _MROCK_COMPILER_VERSION_LINES 0 MROCK_COMPILER_VERSION)
    else()
        set(MROCK_COMPILER_VERSION "unknown")
    endif()

    configure_file(
        "${MROCK_INFO_TEMPLATE}"
        "${MROCK_INFO_OUTPUT}"
        @ONLY
    )

    set(MROCK_INFO_HEADER "${MROCK_INFO_OUTPUT}" PARENT_SCOPE)
    set(MROCK_INFO_INCLUDE_DIR "${_mrock_generated_include_dir}" PARENT_SCOPE)

    if(MROCK_INFO_OUT_HEADER)
        set(${MROCK_INFO_OUT_HEADER} "${MROCK_INFO_OUTPUT}" PARENT_SCOPE)
    endif()

    if(MROCK_INFO_OUT_INCLUDE_DIR)
        set(${MROCK_INFO_OUT_INCLUDE_DIR} "${_mrock_generated_include_dir}" PARENT_SCOPE)
    endif()

    message(STATUS "Generated mrock metadata header: ${MROCK_INFO_OUTPUT}")
endfunction()