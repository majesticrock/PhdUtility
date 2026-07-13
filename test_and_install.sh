#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${ROOT_DIR}/build"

###############################################################################
# Helper
###############################################################################

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

###############################################################################
# Components
###############################################################################

echo
echo "=================================="
echo "        mrock installer"
echo "=================================="
echo

BUILD_UTILITY=OFF
BUILD_SYMBOLIC_OPERATORS=OFF
BUILD_IEOM=OFF

if ask_yes_no "Build utility?" Y; then
    BUILD_UTILITY=ON
fi

if ask_yes_no "Build symbolic_operators?" Y; then
    BUILD_SYMBOLIC_OPERATORS=ON
fi

if ask_yes_no "Build iEoM?" Y; then
    BUILD_IEOM=ON
fi

if [[ "$BUILD_UTILITY" == OFF &&
      "$BUILD_SYMBOLIC_OPERATORS" == OFF &&
      "$BUILD_IEOM" == OFF ]]; then
    echo
    echo "Nothing selected."
    exit 1
fi

###############################################################################
# Tests
###############################################################################

RUN_TESTS=OFF

if ask_yes_no "Run tests?" Y; then
    RUN_TESTS=ON
fi

###############################################################################
# Install
###############################################################################

DO_INSTALL=OFF

if ask_yes_no "Install after build?" Y; then
    DO_INSTALL=ON
fi

INSTALL_PREFIX=""

if [[ "$DO_INSTALL" == ON ]]; then
    read -rp "Install prefix (leave empty for CMake default): " INSTALL_PREFIX
fi

###############################################################################
# Extra include directories
###############################################################################

MROCK_EXTRA_INCLUDE_DIRS=""

if ask_yes_no "Add custom include directories?" N; then

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
# Configure
###############################################################################

echo
echo "Configuring..."

CMAKE_ARGS=(
    -S "${ROOT_DIR}"
    -B "${BUILD_DIR}"
    -DMROCK_BUILD_UTILITY="${BUILD_UTILITY}"
    -DMROCK_BUILD_SYMBOLIC_OPERATORS="${BUILD_SYMBOLIC_OPERATORS}"
    -DMROCK_BUILD_IEOM="${BUILD_IEOM}"
)

if [[ -n "${MROCK_EXTRA_INCLUDE_DIRS}" ]]; then
    CMAKE_ARGS+=(
        "-DMROCK_EXTRA_INCLUDE_DIRS=${MROCK_EXTRA_INCLUDE_DIRS}"
    )
fi

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

if [[ "$RUN_TESTS" == ON ]]; then
    echo
    echo "Running tests..."

    ctest \
        --test-dir "${BUILD_DIR}" \
        --output-on-failure
fi

###############################################################################
# Install
###############################################################################

if [[ "$DO_INSTALL" == ON ]]; then

    echo
    echo "Installing..."

    if [[ -n "$INSTALL_PREFIX" ]]; then
        cmake \
            --install "${BUILD_DIR}" \
            --prefix "${INSTALL_PREFIX}"
    else
        cmake \
            --install "${BUILD_DIR}"
    fi
fi

echo
echo "Done."