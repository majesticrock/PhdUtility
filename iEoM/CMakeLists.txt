cmake_minimum_required(VERSION 3.16)

if(NOT MROCK_VERSION)
    if(TARGET mrock)
        message(WARNING "MROCK_VERSION not set!")
    endif()
    set(MROCK_VERSION 0.0.0)
endif()
project(mrock_iEoM VERSION ${MROCK_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- Library Definition ---
add_library(mrock_iEoM INTERFACE)
add_library(mrock::iEoM ALIAS mrock_iEoM)

target_include_directories(mrock_iEoM INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/mrock_install.cmake)
mrock_install(mrock_iEoM)

if(IEOM_TESTS)
    add_subdirectory(tests)
endif()