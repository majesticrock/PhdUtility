cmake_minimum_required(VERSION 3.16)
project(utilityTests)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(function_time_test function_time.cpp)
SET_COMPILER_FLAGS(function_time_test)
target_link_libraries(function_time_test PRIVATE mrock::utility)
add_executable(cauchy_principal_value_test cauchy_principal_value.cpp)
SET_COMPILER_FLAGS(cauchy_principal_value_test)
target_link_libraries(cauchy_principal_value_test PRIVATE mrock::utility)
add_executable(trapezoidal_test trapezoidal_test.cpp)
SET_COMPILER_FLAGS(trapezoidal_test)
target_link_libraries(trapezoidal_test PRIVATE mrock::utility)
add_executable(hypergeometric_test hypergeometric.cpp)
SET_COMPILER_FLAGS(hypergeometric_test)
target_link_libraries(hypergeometric_test PRIVATE mrock::utility)

# Enable CTest
enable_testing()

# Add a test
add_test(NAME function_time_test COMMAND function_time_test)
add_test(NAME cauchy_principal_value_test COMMAND cauchy_principal_value_test)
add_test(NAME trapezoidal_test COMMAND trapezoidal_test)
add_test(NAME hypergeometric_test COMMAND hypergeometric_test)