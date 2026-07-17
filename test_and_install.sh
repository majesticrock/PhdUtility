#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${ROOT_DIR}/build"
USER_CONFIG="${ROOT_DIR}/cmake/mrock-user.cmake"

BEGIN_MARKER="# >>> mrock installer managed settings >>>"
END_MARKER="# <<< mrock installer managed settings <<<"

###############################################################################
# Helper
###############################################################################

is_cmake_true() {
    local value="${1:-OFF}"
    value="${value^^}"

    case "${value}" in
        ON|YES|TRUE|Y|1)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

normalize_bool() {
    if is_cmake_true "${1:-OFF}"; then
        echo "ON"
    else
        echo "OFF"
    fi
}

default_letter() {
    if is_cmake_true "${1:-OFF}"; then
        echo "Y"
    else
        echo "N"
    fi
}

ask_yes_no() {
    local prompt="$1"
    local default="${2:-Y}"
    local reply

    while true; do
        if [[ "$default" == "Y" ]]; then
            read -rp "${prompt} [Y/n]: " reply
            reply="${reply:-Y}"
        else
            read -rp "${prompt} [y/N]: " reply
            reply="${reply:-N}"
        fi

        case "${reply}" in
            [Yy]*) return 0 ;;
            [Nn]*) return 1 ;;
        esac

        echo "Please answer yes or no."
    done
}

cmake_escape_string() {
    local value="$1"

    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"

    printf '%s' "${value}"
}

###############################################################################
# Read current CMake/user config defaults
###############################################################################

CURRENT_MROCK_BUILD_UTILITY="ON"
CURRENT_MROCK_BUILD_SYMBOLIC_OPERATORS="ON"
CURRENT_MROCK_BUILD_IEOM="ON"
CURRENT_BUILD_TESTING="ON"
CURRENT_MROCK_EXTRA_INCLUDE_DIRS=""
CURRENT_CMAKE_INSTALL_PREFIX=""

read_user_config_defaults() {
    if [[ ! -f "${USER_CONFIG}" ]]; then
        return
    fi

    local tmp
    tmp="$(mktemp)"

    local root_escaped
    local user_config_escaped

    root_escaped="$(cmake_escape_string "${ROOT_DIR}")"
    user_config_escaped="$(cmake_escape_string "${USER_CONFIG}")"

    cat > "${tmp}" <<EOF
set(CMAKE_SOURCE_DIR "${root_escaped}")
set(PROJECT_SOURCE_DIR "${root_escaped}")

set(MROCK_BUILD_UTILITY ON CACHE BOOL "" FORCE)
set(MROCK_BUILD_SYMBOLIC_OPERATORS ON CACHE BOOL "" FORCE)
set(MROCK_BUILD_IEOM ON CACHE BOOL "" FORCE)
set(BUILD_TESTING ON CACHE BOOL "" FORCE)
set(MROCK_EXTRA_INCLUDE_DIRS "" CACHE STRING "" FORCE)
set(CMAKE_INSTALL_PREFIX "" CACHE PATH "" FORCE)

include("${user_config_escaped}")

foreach(var IN ITEMS
    MROCK_BUILD_UTILITY
    MROCK_BUILD_SYMBOLIC_OPERATORS
    MROCK_BUILD_IEOM
    BUILD_TESTING
    MROCK_EXTRA_INCLUDE_DIRS
    CMAKE_INSTALL_PREFIX
)
    if(DEFINED \${var})
        message("MROCK_DUMP:\${var}=[\${\${var}}]")
    endif()
endforeach()
EOF

local dump_regex
dump_regex='MROCK_DUMP:([^=]+)=
$$
(.*)
$$
'

while IFS= read -r line; do
    if [[ "${line}" =~ ${dump_regex} ]]; then
        local var="${BASH_REMATCH[1]}"
        local value="${BASH_REMATCH[2]}"

        case "${var}" in
            MROCK_BUILD_UTILITY)
                CURRENT_MROCK_BUILD_UTILITY="${value}"
                ;;
            MROCK_BUILD_SYMBOLIC_OPERATORS)
                CURRENT_MROCK_BUILD_SYMBOLIC_OPERATORS="${value}"
                ;;
            MROCK_BUILD_IEOM)
                CURRENT_MROCK_BUILD_IEOM="${value}"
                ;;
            BUILD_TESTING)
                CURRENT_BUILD_TESTING="${value}"
                ;;
            MROCK_EXTRA_INCLUDE_DIRS)
                CURRENT_MROCK_EXTRA_INCLUDE_DIRS="${value}"
                ;;
            CMAKE_INSTALL_PREFIX)
                CURRENT_CMAKE_INSTALL_PREFIX="${value}"
                ;;
        esac
    fi
done < <(cmake -P "${tmp}" 2>&1)

    rm -f "${tmp}"

    CURRENT_MROCK_BUILD_UTILITY="$(normalize_bool "${CURRENT_MROCK_BUILD_UTILITY}")"
    CURRENT_MROCK_BUILD_SYMBOLIC_OPERATORS="$(normalize_bool "${CURRENT_MROCK_BUILD_SYMBOLIC_OPERATORS}")"
    CURRENT_MROCK_BUILD_IEOM="$(normalize_bool "${CURRENT_MROCK_BUILD_IEOM}")"
    CURRENT_BUILD_TESTING="$(normalize_bool "${CURRENT_BUILD_TESTING}")"
}

###############################################################################
# Generate/update user config
###############################################################################

generate_user_config_block() {
    local extra_include_dirs_escaped
    extra_include_dirs_escaped="$(cmake_escape_string "${MROCK_EXTRA_INCLUDE_DIRS}")"

    local install_prefix_escaped
    install_prefix_escaped="$(cmake_escape_string "${INSTALL_PREFIX}")"

    cat <<EOF
${BEGIN_MARKER}
# This block is managed by the mrock installer script.
# You may edit this file manually, but changes inside this block may be
# overwritten by the installer.

set(MROCK_BUILD_UTILITY ${BUILD_UTILITY} CACHE BOOL "Build the utility component" FORCE)
set(MROCK_BUILD_SYMBOLIC_OPERATORS ${BUILD_SYMBOLIC_OPERATORS} CACHE BOOL "Build the symbolic operators component" FORCE)
set(MROCK_BUILD_IEOM ${BUILD_IEOM} CACHE BOOL "Build the iEoM component" FORCE)

set(BUILD_TESTING ${BUILD_TESTING} CACHE BOOL "Build tests" FORCE)

set(CMAKE_INSTALL_PREFIX "${install_prefix_escaped}" CACHE PATH "Install prefix" FORCE)

set(MROCK_EXTRA_INCLUDE_DIRS "${extra_include_dirs_escaped}" CACHE STRING "Additional include directories for mrock, semicolon-separated" FORCE)
${END_MARKER}
EOF
}

write_user_config_if_needed() {
    local changed="OFF"

    if [[ "${BUILD_UTILITY}" != "${CURRENT_MROCK_BUILD_UTILITY}" ]]; then
        changed="ON"
    fi

    if [[ "${BUILD_SYMBOLIC_OPERATORS}" != "${CURRENT_MROCK_BUILD_SYMBOLIC_OPERATORS}" ]]; then
        changed="ON"
    fi

    if [[ "${BUILD_IEOM}" != "${CURRENT_MROCK_BUILD_IEOM}" ]]; then
        changed="ON"
    fi

    if [[ "${BUILD_TESTING}" != "${CURRENT_BUILD_TESTING}" ]]; then
        changed="ON"
    fi

    if [[ "${MROCK_EXTRA_INCLUDE_DIRS}" != "${CURRENT_MROCK_EXTRA_INCLUDE_DIRS}" ]]; then
        changed="ON"
    fi

    if [[ "${INSTALL_PREFIX}" != "${CURRENT_CMAKE_INSTALL_PREFIX}" ]]; then
        changed="ON"
    fi

    if [[ -f "${USER_CONFIG}" && "${changed}" == "OFF" ]]; then
        echo "User config unchanged: ${USER_CONFIG}"
        return
    fi

    mkdir -p "$(dirname "${USER_CONFIG}")"

    local new_block
    new_block="$(generate_user_config_block)"

    if [[ ! -f "${USER_CONFIG}" ]]; then
        printf '%s\n' "${new_block}" > "${USER_CONFIG}"
        echo "Created user config: ${USER_CONFIG}"
        return
    fi

    local tmp
    tmp="$(mktemp)"

    if grep -Fq "${BEGIN_MARKER}" "${USER_CONFIG}" &&
       grep -Fq "${END_MARKER}" "${USER_CONFIG}"; then

        awk \
            -v begin="${BEGIN_MARKER}" \
            -v end="${END_MARKER}" \
            -v block="${new_block}" '
                $0 == begin {
                    print block
                    in_block = 1
                    next
                }

                $0 == end {
                    in_block = 0
                    next
                }

                !in_block {
                    print
                }
            ' "${USER_CONFIG}" > "${tmp}"

        mv "${tmp}" "${USER_CONFIG}"
        echo "Updated managed block in user config: ${USER_CONFIG}"
    else
        {
            cat "${USER_CONFIG}"
            echo
            printf '%s\n' "${new_block}"
        } > "${tmp}"

        mv "${tmp}" "${USER_CONFIG}"
        echo "Appended managed block to user config: ${USER_CONFIG}"
    fi
}

###############################################################################
# Main
###############################################################################

read_user_config_defaults

echo
echo "=================================="
echo "        mrock installer"
echo "=================================="
echo
echo "User config: ${USER_CONFIG}"
echo

###############################################################################
# Components
###############################################################################

if ask_yes_no "Build utility?" "$(default_letter "${CURRENT_MROCK_BUILD_UTILITY}")"; then
    BUILD_UTILITY="ON"
else
    BUILD_UTILITY="OFF"
fi

if ask_yes_no "Build symbolic_operators?" "$(default_letter "${CURRENT_MROCK_BUILD_SYMBOLIC_OPERATORS}")"; then
    BUILD_SYMBOLIC_OPERATORS="ON"
else
    BUILD_SYMBOLIC_OPERATORS="OFF"
fi

if ask_yes_no "Build iEoM?" "$(default_letter "${CURRENT_MROCK_BUILD_IEOM}")"; then
    BUILD_IEOM="ON"
else
    BUILD_IEOM="OFF"
fi

if [[ "${BUILD_UTILITY}" == "OFF" &&
      "${BUILD_SYMBOLIC_OPERATORS}" == "OFF" &&
      "${BUILD_IEOM}" == "OFF" ]]; then
    echo
    echo "Nothing selected."
    exit 1
fi

###############################################################################
# Tests
###############################################################################

if ask_yes_no "Build tests?" "$(default_letter "${CURRENT_BUILD_TESTING}")"; then
    BUILD_TESTING="ON"
else
    BUILD_TESTING="OFF"
fi

RUN_TESTS="OFF"

if [[ "${BUILD_TESTING}" == "ON" ]]; then
    if ask_yes_no "Run tests after build?" Y; then
        RUN_TESTS="ON"
    fi
else
    echo "Tests are disabled, so ctest will not be run."
fi

###############################################################################
# Install
###############################################################################

INSTALL_PREFIX="${CURRENT_CMAKE_INSTALL_PREFIX}"

echo
if [[ -n "${CURRENT_CMAKE_INSTALL_PREFIX}" ]]; then
    echo "Current install prefix: ${CURRENT_CMAKE_INSTALL_PREFIX}"
else
    echo "Current install prefix: <CMake default>"
fi

if ask_yes_no "Change install prefix?" N; then
    read -rp "Install prefix, leave empty for CMake default: " INSTALL_PREFIX
fi

DO_INSTALL="OFF"

if ask_yes_no "Install after build?" Y; then
    DO_INSTALL="ON"
fi

###############################################################################
# Extra include directories
###############################################################################

MROCK_EXTRA_INCLUDE_DIRS="${CURRENT_MROCK_EXTRA_INCLUDE_DIRS}"

echo
if [[ -n "${CURRENT_MROCK_EXTRA_INCLUDE_DIRS}" ]]; then
    echo "Current custom include directories:"
    echo "  ${CURRENT_MROCK_EXTRA_INCLUDE_DIRS}"
else
    echo "Current custom include directories: <none>"
fi

if ask_yes_no "Replace custom include directories?" N; then
    MROCK_EXTRA_INCLUDE_DIRS=""

    echo
    echo "Enter include directories one by one."
    echo "Leave empty to finish."

    while true; do
        read -rp "Include directory: " INCLUDE_DIR

        if [[ -z "${INCLUDE_DIR}" ]]; then
            break
        fi

        if [[ -n "${MROCK_EXTRA_INCLUDE_DIRS}" ]]; then
            MROCK_EXTRA_INCLUDE_DIRS="${MROCK_EXTRA_INCLUDE_DIRS};${INCLUDE_DIR}"
        else
            MROCK_EXTRA_INCLUDE_DIRS="${INCLUDE_DIR}"
        fi
    done
fi

###############################################################################
# Write user config
###############################################################################

write_user_config_if_needed

###############################################################################
# Configure
###############################################################################

echo
echo "Configuring..."

CMAKE_ARGS=(
    -S "${ROOT_DIR}"
    -B "${BUILD_DIR}"
    -DMROCK_USER_CONFIG="${USER_CONFIG}"
    -DMROCK_BUILD_UTILITY="${BUILD_UTILITY}"
    -DMROCK_BUILD_SYMBOLIC_OPERATORS="${BUILD_SYMBOLIC_OPERATORS}"
    -DMROCK_BUILD_IEOM="${BUILD_IEOM}"
    -DBUILD_TESTING="${BUILD_TESTING}"
    "-DMROCK_EXTRA_INCLUDE_DIRS=${MROCK_EXTRA_INCLUDE_DIRS}"
    "-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}"
)

cmake "${CMAKE_ARGS[@]}"

###############################################################################
# Build
###############################################################################

echo
echo "Building..."

cmake --build "${BUILD_DIR}" --parallel

###############################################################################
# Test
###############################################################################

if [[ "${RUN_TESTS}" == "ON" ]]; then
    echo
    echo "Running tests..."

    ctest \
        --test-dir "${BUILD_DIR}" \
        --output-on-failure
fi

###############################################################################
# Install
###############################################################################

if [[ "${DO_INSTALL}" == "ON" ]]; then
    echo
    echo "Installing..."

    cmake --install "${BUILD_DIR}"
fi

echo
echo "Done."