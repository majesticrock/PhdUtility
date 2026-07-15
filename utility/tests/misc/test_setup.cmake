# utility/tests/misc/test_setup.cmake
include_guard()

add_executable(better_to_string_test ${CMAKE_CURRENT_LIST_DIR}/better_to_string.cpp)
mrock_set_build_options(better_to_string_test)
target_link_libraries(better_to_string_test 
    PRIVATE 
        mrock::utility 
)

add_executable(BinaryIO_test ${CMAKE_CURRENT_LIST_DIR}/BinaryIO.cpp)
mrock_set_build_options(BinaryIO_test)
target_link_libraries(BinaryIO_test 
    PRIVATE 
        mrock::utility 
)

add_executable(ComplexNumberIterators_test ${CMAKE_CURRENT_LIST_DIR}/ComplexNumberIterators.cpp)
mrock_set_build_options(ComplexNumberIterators_test)
target_link_libraries(ComplexNumberIterators_test 
    PRIVATE 
        mrock::utility 
)

add_executable(constexpr_power_test ${CMAKE_CURRENT_LIST_DIR}/constexpr_power.cpp)
mrock_set_build_options(constexpr_power_test)
target_link_libraries(constexpr_power_test 
    PRIVATE 
        mrock::utility 
)

add_executable(function_time_test ${CMAKE_CURRENT_LIST_DIR}/function_time.cpp)
mrock_set_build_options(function_time_test)
target_link_libraries(function_time_test 
    PRIVATE 
        mrock::utility 
)

add_executable(UnderlyingRealType_test ${CMAKE_CURRENT_LIST_DIR}/UnderlyingRealType.cpp)
mrock_set_build_options(UnderlyingRealType_test)
target_link_libraries(UnderlyingRealType_test 
    PRIVATE 
        mrock::utility 
)

# Add a test
add_test(NAME better_to_string_test COMMAND better_to_string_test)
add_test(NAME BinaryIO_test COMMAND BinaryIO_test)
add_test(NAME ComplexNumberIterators_test COMMAND ComplexNumberIterators_test)
add_test(NAME constexpr_power_test COMMAND constexpr_power_test)
add_test(NAME function_time_test COMMAND function_time_test)
add_test(NAME UnderlyingRealType_test COMMAND UnderlyingRealType_test)
