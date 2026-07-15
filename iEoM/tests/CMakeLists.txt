# iEoM/tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(iEoM_tests)

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/mrock_iEoM_extra_options.cmake)

add_executable(block_matrix_test block_matrix.cpp)
mrock_set_build_options(block_matrix_test)
target_link_libraries(block_matrix_test
    PRIVATE 
        mrock::iEoM 
        mrock_iEoM_extra_options
)

add_executable(ieom_bcs_test ieom_bcs.cpp)
mrock_set_build_options(ieom_bcs_test)
target_link_libraries(ieom_bcs_test
    PRIVATE 
        mrock::iEoM 
        mrock_iEoM_extra_options
)

# Enable CTest
enable_testing()

# Add a test
add_test(NAME block_matrix_test COMMAND block_matrix_test)
add_test(NAME ieom_bcs_test COMMAND ieom_bcs_test)