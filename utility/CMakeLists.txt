cmake_minimum_required(VERSION 3.16)

if(NOT MROCK_VERSION)
    if(TARGET mrock)
        message(WARNING "MROCK_VERSION not set!")
    endif()
    set(MROCK_VERSION 0.0.0)
endif()
project(mrock_utility VERSION ${MROCK_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- Dependencies ---
include(CMakeFindDependencyMacro)
find_dependency(ZLIB REQUIRED)
if(CLUSTER_BUILD) # Force boost 1.78 on the cluster - there are conflicts on that system specifically
    find_dependency(Boost 1.78 EXACT REQUIRED COMPONENTS serialization iostreams)
else()
    find_dependency(Boost REQUIRED COMPONENTS serialization iostreams)
endif()

# --- Library Definition ---
add_library(mrock_utility INTERFACE)
add_library(mrock::utility ALIAS mrock_utility)

target_include_directories(mrock_utility INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(mrock_utility INTERFACE 
    ZLIB::ZLIB 
    Boost::serialization
    Boost::iostreams
)

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/mrock_install.cmake)
mrock_install(mrock_utility)

if(UTIL_TESTS)
    add_subdirectory(tests)
endif()