# utility/tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(utilityTests)

add_executable(function_time_test function_time.cpp)
mrock_set_build_options(function_time_test)
target_link_libraries(function_time_test 
    PRIVATE 
        mrock::utility 
)

add_executable(cauchy_principal_value_test cauchy_principal_value.cpp)
mrock_set_build_options(cauchy_principal_value_test)
target_link_libraries(cauchy_principal_value_test 
    PRIVATE 
        mrock::utility 
)

add_executable(trapezoidal_test trapezoidal_test.cpp)
mrock_set_build_options(trapezoidal_test)
target_link_libraries(trapezoidal_test 
    PRIVATE 
        mrock::utility 
)

add_executable(hypergeometric_test hypergeometric.cpp)
mrock_set_build_options(hypergeometric_test)
target_link_libraries(hypergeometric_test 
    PRIVATE 
        mrock::utility 
)

# Enable CTest
enable_testing()

# Add a test
add_test(NAME function_time_test COMMAND function_time_test)
add_test(NAME cauchy_principal_value_test COMMAND cauchy_principal_value_test)
add_test(NAME trapezoidal_test COMMAND trapezoidal_test)
add_test(NAME hypergeometric_test COMMAND hypergeometric_test)