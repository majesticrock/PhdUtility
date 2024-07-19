cmake_minimum_required(VERSION 3.16)
project(MyProjectTests)

# Add test executable
add_executable(minimization_test minimization.cpp)

# Enable CTest
enable_testing()

# Add a test
add_test(NAME minimization_test COMMAND minimization_test)
