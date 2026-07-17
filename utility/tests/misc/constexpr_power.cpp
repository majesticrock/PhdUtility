#include <mrock/utility/constexpr_power.hpp>

#include <cmath>
#include <iostream>
#include <type_traits>

#define FAIL_IF_NOT(cond)                                                      \
    if (!(cond)) {                                                             \
        std::cerr << "FAILED: " #cond << " at line " << __LINE__ << std::endl; \
        return 1;                                                              \
    }

int main() {
    using mrock::utility::constexpr_power;

    constexpr double tolerance = 1e-12;

    {
        std::cout << "Test 1: positive exponent\n";

        constexpr double result = constexpr_power<3>(2.0);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 8\n\n";

        static_assert(result == 8.0);
        FAIL_IF_NOT(result == 8.0);
    }

    {
        std::cout << "Test 2: exponent zero\n";

        constexpr double result = constexpr_power<0>(123.456);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 1\n\n";

        static_assert(result == 1.0);
        FAIL_IF_NOT(result == 1.0);
    }

    {
        std::cout << "Test 3: negative exponent\n";

        constexpr double result = constexpr_power<-2>(4.0);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 0.0625\n\n";

        static_assert(result == 0.0625);
        FAIL_IF_NOT(result == 0.0625);
    }

    {
        std::cout << "Test 4: integer base with default ResultType\n";

        constexpr double result = constexpr_power<5>(2);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 32\n\n";

        static_assert(result == 32.0);
        static_assert(std::is_same_v<decltype(result), const double>);
        FAIL_IF_NOT(result == 32.0);
    }

    {
        std::cout << "Test 5: custom ResultType\n";

        constexpr int result = constexpr_power<4, int, int>(3);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 81\n\n";

        static_assert(result == 81);
        static_assert(std::is_same_v<decltype(result), const int>);
        FAIL_IF_NOT(result == 81);
    }

    {
        std::cout << "Test 6: floating-point base\n";

        constexpr double result = constexpr_power<3>(1.5);
        constexpr double expected = 3.375;

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: " << expected << "\n\n";

        static_assert(result == expected);
        FAIL_IF_NOT(std::abs(result - expected) < tolerance);
    }

    {
        std::cout << "Test 7: runtime value\n";

        double base = 2.5;
        double result = constexpr_power<3>(base);
        double expected = 15.625;

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: " << expected << "\n\n";

        FAIL_IF_NOT(std::abs(result - expected) < tolerance);
    }

    {
        std::cout << "Test 8: compare with std::pow\n";

        double base = 1.2345;
        double result = constexpr_power<6>(base);
        double expected = std::pow(base, 6);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: " << expected << "\n\n";

        FAIL_IF_NOT(std::abs(result - expected) < tolerance);
    }

    {
        std::cout << "Test 9: negative exponent compared with std::pow\n";

        double base = 2.5;
        double result = constexpr_power<-3>(base);
        double expected = std::pow(base, -3);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: " << expected << "\n\n";

        FAIL_IF_NOT(std::abs(result - expected) < tolerance);
    }

    std::cout << "All constexpr_power tests passed.\n";

    return 0;
}