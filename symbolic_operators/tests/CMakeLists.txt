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


enable_testing()

add_test(NAME symop_continuum_test COMMAND symop_continuum_test)
add_test(NAME symop_bosons_test COMMAND symop_bosons_test)
add_test(NAME symop_bcs_test COMMAND symop_bcs_test)
add_test(NAME symop_string_helper_test COMMAND symop_string_helper_test)
add_test(NAME symop_matrix_representation_test COMMAND symop_matrix_representation_test)