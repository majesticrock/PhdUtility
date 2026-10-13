# symbolic_operators/tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.30)
project(symbolic_operatorsTests)

add_executable(symop_continuum_test continuum.cpp)
mrock_set_build_options(symop_continuum_test)
target_link_libraries(symop_continuum_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(symop_bosons_test bosons.cpp)
mrock_set_build_options(symop_bosons_test)
target_link_libraries(symop_bosons_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(symop_bcs_test bcs.cpp)
mrock_set_build_options(symop_bcs_test)
target_link_libraries(symop_bcs_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(symop_string_helper_test string_helper.cpp)
mrock_set_build_options(symop_string_helper_test)
target_link_libraries(symop_string_helper_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(symop_matrix_representation_test matrix_representation.cpp)
mrock_set_build_options(symop_matrix_representation_test)
target_link_libraries(symop_matrix_representation_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(explicit_wick_test explicit_wick.cpp)
mrock_set_build_options(explicit_wick_test)
target_link_libraries(explicit_wick_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(symop_wick_multiple_indices_test wick_multiple_indices.cpp)
mrock_set_build_options(symop_wick_multiple_indices_test)
target_link_libraries(symop_wick_multiple_indices_test
    PRIVATE
        mrock::symbolic_operators
)

add_executable(symop_wick_ordered_multiple_indices_test wick_ordered_multiple_indices.cpp)
mrock_set_build_options(symop_wick_ordered_multiple_indices_test)
target_link_libraries(symop_wick_ordered_multiple_indices_test
    PRIVATE
        mrock::symbolic_operators
)

add_executable(symop_wick_quartic_test wick_quartic.cpp)
mrock_set_build_options(symop_wick_quartic_test)
target_link_libraries(symop_wick_quartic_test
    PRIVATE
        mrock::symbolic_operators
)

add_executable(symop_wick_normal_order_test wick_normal_order.cpp)
mrock_set_build_options(symop_wick_normal_order_test)
target_link_libraries(symop_wick_normal_order_test
    PRIVATE
        mrock::symbolic_operators
)

add_executable(symop_wick_ordered_multiply_test wick_ordered_multiply.cpp)
mrock_set_build_options(symop_wick_ordered_multiply_test)
target_link_libraries(symop_wick_ordered_multiply_test
    PRIVATE
        mrock::symbolic_operators
)

add_executable(symop_wick_operator_conversion_test wick_operator_conversion.cpp)
mrock_set_build_options(symop_wick_operator_conversion_test)
target_link_libraries(symop_wick_operator_conversion_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(symop_wick_ordered_test wick_ordered.cpp)
mrock_set_build_options(symop_wick_ordered_test)
target_link_libraries(symop_wick_ordered_test
    PRIVATE
        mrock::symbolic_operators
)

add_executable(symop_wick_matrix_test wick_matrix.cpp)
mrock_set_build_options(symop_wick_matrix_test)
target_link_libraries(symop_wick_matrix_test
    PRIVATE
        mrock::symbolic_operators
)

add_executable(symop_wick_commutator_test wick_commutator.cpp)
mrock_set_build_options(symop_wick_commutator_test)
target_link_libraries(symop_wick_commutator_test
    PRIVATE
        mrock::symbolic_operators
)

enable_testing()

add_test(NAME symop_continuum_test COMMAND symop_continuum_test)
add_test(NAME symop_bosons_test COMMAND symop_bosons_test)
add_test(NAME symop_bcs_test COMMAND symop_bcs_test)
add_test(NAME symop_string_helper_test COMMAND symop_string_helper_test)
add_test(NAME symop_matrix_representation_test COMMAND symop_matrix_representation_test)
add_test(NAME explicit_wick_test COMMAND explicit_wick_test)
add_test(NAME symop_wick_multiple_indices_test COMMAND symop_wick_multiple_indices_test)
add_test(NAME symop_wick_ordered_multiple_indices_test COMMAND symop_wick_ordered_multiple_indices_test)
add_test(NAME symop_wick_quartic_test COMMAND symop_wick_quartic_test)
add_test(NAME symop_wick_normal_order_test COMMAND symop_wick_normal_order_test)
add_test(NAME symop_wick_ordered_multiply_test COMMAND symop_wick_ordered_multiply_test)
add_test(NAME symop_wick_operator_conversion_test COMMAND symop_wick_operator_conversion_test)
add_test(NAME symop_wick_ordered_test COMMAND symop_wick_ordered_test)
add_test(NAME symop_wick_matrix_test COMMAND symop_wick_matrix_test)
add_test(NAME symop_wick_commutator_test COMMAND symop_wick_commutator_test)