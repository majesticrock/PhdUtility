#include "../include/mrock/Utility/FunctionTime.hpp"
#include <cmath>
#include <iostream>

double test_cos(double x, int N) {
    double y{};
    for (int i = 0; i < N; ++i) {
        y += std::cos(x + static_cast<double>(i));
    }
    return y;
}

double test_cos_2() {
    double x = 1.2;
    for (int i = 0; i < 100000; ++i) {
        x += std::cos(x + static_cast<double>(i));
    }
    return x;
}

struct TestStruct {
    double x{};

    double member_test(int N) {
        for (int i = 0; i < N; ++i) {
            x += std::cos(x + static_cast<double>(i));
        }
        return x;
    }

    double member_test_2() {
        for (int i = 0; i < 100000; ++i) {
            x += std::cos(x + static_cast<double>(i));
        }
        return x;
    }
};

using namespace mrock::Utility;

int main() {
    const double x = 0.1;
    const int N = 100000;

    // Calling f(args)
    const double test_cos_result = function_time_micro(&test_cos, x, N);
    const double base_line = test_cos(x, N);
    if(std::abs(test_cos_result - base_line) > 1e-12) return 1;
    // Calling f()
    const double test_cos_result_2 = function_time_micro(&test_cos_2);
    const double base_line_2 = test_cos_2();
    if(std::abs(test_cos_result_2 - base_line_2) > 1e-12) return 1;

    // Calling class::f(args)
    TestStruct test_struct;
    const double test_struct_result = member_function_time_micro(test_struct, &TestStruct::member_test, N);
    TestStruct base_struct;
    const double base_struct_result = base_struct.member_test(N);
    if(std::abs(test_struct_result - base_struct_result) > 1e-12) return 1;
    // Calling class::f()
    const double test_struct_result_2 = member_function_time_micro(test_struct, &TestStruct::member_test_2);
    const double base_struct_result_2 = base_struct.member_test_2();
    if(std::abs(test_struct_result_2 - base_struct_result_2) > 1e-12) return 1;

    return 0;
}