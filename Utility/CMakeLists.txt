cmake_minimum_required(VERSION 3.16)

project(mrock_utility LANGUAGES CXX VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(MROCK_UTILITY_HEADERS include/mrock/utility)
# --- Dependencies ---
include(CMakeFindDependencyMacro)
find_dependency(ZLIB REQUIRED)
find_dependency(Boost REQUIRED COMPONENTS serialization iostreams)

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

include(${CMAKE_CURRENT_LIST_DIR}/cmake/mrock_utility-information.cmake)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/utility_info.h" DESTINATION ${MROCK_UTILITY_HEADERS})

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/mrock_install.cmake)
mrock_install(mrock_utility)

if(UTIL_TESTS)
    add_subdirectory(tests)
endif()