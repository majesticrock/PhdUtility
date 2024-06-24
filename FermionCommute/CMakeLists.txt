# FermionCommute/CMakeLists.txt

cmake_minimum_required(VERSION 3.16)

project(FermionCommute VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")

add_subdirectory(sources/Definitions)
include_directories(
    $ENV{HOME}/usr/local/include
)

add_executable(FermionCommute 
    sources/FermionCommute.cpp
)

# Check if using GCC and adjust C++ standard flag if necessary
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # Check GCC version
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.1)
        message(FATAL_ERROR "GCC version ${CMAKE_CXX_COMPILER_VERSION} is not supported. GCC 12.1 or newer is required for C++20 support.")
    endif()
    # Use the appropriate flag for GCC
    target_compile_options(FermionCommute PRIVATE
        -Wall
        -Wno-sign-compare
        -fopenmp
        -march=native
        -O3
    )
else()
    message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only GCC is supported.")
endif()


# Find Boost libraries
find_package(Boost REQUIRED COMPONENTS serialization)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost libraries not found")
endif()
target_include_directories(FermionCommute PRIVATE ${Boost_INCLUDE_DIRS})

target_link_libraries(FermionCommute PRIVATE Definitions ${Boost_LIBRARIES})