set(DEF_SOURCES
    Continuum.cpp
    Hubbard.cpp
    HubbardDispersions.cpp
    StandardOperators.cpp
)

add_library(Definitions STATIC ${DEF_SOURCES})

# Include directories for Definitions library
target_include_directories(Definitions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Find and include SymbolicOperators
find_package(SymbolicOperators REQUIRED CONFIG)
target_include_directories(Definitions PRIVATE ${SymbolicOperators_INCLUDE_DIRS})
target_link_libraries(Definitions PUBLIC SymbolicOperators::SymbolicOperators)

find_package(Utility REQUIRED)
target_include_directories(Definitions PRIVATE ${Utility_INCLUDE_DIRS})
target_link_libraries(Definitions PUBLIC Utility::Utility)
