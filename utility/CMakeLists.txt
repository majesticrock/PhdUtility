cmake_minimum_required(VERSION 3.16)

project(mrock_utility VERSION ${MROCK_VERSION} LANGUAGES CXX)

include(CMakeFindDependencyMacro)

find_dependency(ZLIB REQUIRED)
find_dependency(Boost REQUIRED COMPONENTS iostreams)

add_library(mrock_utility INTERFACE)

add_library(
    mrock::utility
    ALIAS
    mrock_utility
)

target_link_libraries(
    mrock_utility
    INTERFACE
        mrock::mrock
        ZLIB::ZLIB
        Boost::iostreams
)

target_include_directories(
    mrock_utility
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

include(
    ${CMAKE_CURRENT_LIST_DIR}/../cmake/mrock_install.cmake
)

mrock_install(mrock_utility)

if(UTIL_TESTS)
    add_subdirectory(tests)
endif()