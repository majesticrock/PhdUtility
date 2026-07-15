# symbolic_operators/tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(symbolic_operatorsTests)

add_executable(continuum_test continuum.cpp)
mrock_set_build_options(continuum_test)
target_link_libraries(continuum_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(bosons_test bosons.cpp)
mrock_set_build_options(bosons_test)
target_link_libraries(bosons_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(bcs_test bcs.cpp)
mrock_set_build_options(bcs_test)
target_link_libraries(bcs_test 
    PRIVATE 
        mrock::symbolic_operators 
)

add_executable(string_helper_test string_helper.cpp)
mrock_set_build_options(string_helper_test)
target_link_libraries(string_helper_test 
    PRIVATE 
        mrock::symbolic_operators 
)

enable_testing()

add_test(NAME continuum_test COMMAND continuum_test)
add_test(NAME bosons_test COMMAND bosons_test)
add_test(NAME bcs_test COMMAND bcs_test)
add_test(NAME string_helper_test COMMAND string_helper_test)