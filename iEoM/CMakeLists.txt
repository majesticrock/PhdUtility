cmake_minimum_required(VERSION 3.16)

project(
    mrock_iEoM
    VERSION ${MROCK_VERSION}
    LANGUAGES CXX
)

add_library(mrock_iEoM INTERFACE)
add_library(
    mrock::iEoM
    ALIAS
    mrock_iEoM
)
set_target_properties(mrock_iEoM PROPERTIES
    EXPORT_NAME iEoM
)

target_link_libraries(
    mrock_iEoM
    INTERFACE
        mrock::mrock
)

target_include_directories(
    mrock_iEoM
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
        $<INSTALL_INTERFACE:${MROCK_INSTALL_INCLUDEDIR}>
)

include(
    ${CMAKE_CURRENT_LIST_DIR}/../cmake/mrock_install.cmake
)

mrock_install(mrock_iEoM)

if(IEOM_TESTS)
    add_subdirectory(tests)
endif()