cmake_minimum_required(VERSION 3.16)
project(UtilityTests)

# Add test executable
add_executable(minimization_test minimization.cpp)
add_executable(block_matrix_test block_matrix.cpp)

# Enable CTest
enable_testing()

# Add a test
add_test(NAME minimization_test COMMAND minimization_test)
add_test(NAME block_matrix_test COMMAND block_matrix_test)