# utility/tests/Numerics/test_setup.cmake
include_guard()

add_executable(cauchy_principal_value_test ${CMAKE_CURRENT_LIST_DIR}/cauchy_principal_value.cpp)
mrock_set_build_options(cauchy_principal_value_test)
target_link_libraries(cauchy_principal_value_test 
    PRIVATE 
        mrock::utility 
)

add_executable(trapezoidal_test ${CMAKE_CURRENT_LIST_DIR}/trapezoidal.cpp)
mrock_set_build_options(trapezoidal_test)
target_link_libraries(trapezoidal_test 
    PRIVATE 
        mrock::utility 
)

add_executable(bisection_test ${CMAKE_CURRENT_LIST_DIR}/bisection.cpp)
mrock_set_build_options(bisection_test)
target_link_libraries(bisection_test 
    PRIVATE 
        mrock::utility 
)

add_executable(broydens_method_test ${CMAKE_CURRENT_LIST_DIR}/broydens_method.cpp)
mrock_set_build_options(broydens_method_test)
target_link_libraries(broydens_method_test 
    PRIVATE 
        mrock::utility 
)

add_executable(interpolation_test ${CMAKE_CURRENT_LIST_DIR}/interpolation.cpp)
mrock_set_build_options(interpolation_test)
target_link_libraries(interpolation_test 
    PRIVATE 
        mrock::utility 
)

add_executable(hypergeometric_test ${CMAKE_CURRENT_LIST_DIR}/hypergeometric.cpp)
mrock_set_build_options(hypergeometric_test)
target_link_libraries(hypergeometric_test 
    PRIVATE 
        mrock::utility 
)

# Add a test
add_test(NAME cauchy_principal_value_test COMMAND cauchy_principal_value_test)
add_test(NAME trapezoidal_test COMMAND trapezoidal_test)
add_test(NAME bisection_test COMMAND bisection_test)
add_test(NAME broydens_method_test COMMAND broydens_method_test)
add_test(NAME interpolation_test COMMAND interpolation_test)
add_test(NAME hypergeometric_test COMMAND hypergeometric_test)