#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include <mrock/utility/Numerics/Interpolation.hpp>

#define FAIL_IF_NOT(cond)                                      \
    if (!(cond)) {                                             \
        std::cerr << "FAILED: " #cond                          \
                  << " at line " << __LINE__ << std::endl;     \
        return 1;                                              \
    }

int main() {
    using namespace mrock::utility::Numerics;

    constexpr double tolerance = 1e-12;

    {
        std::cout << "Test 1: linear interpolation\n";

        double x = 0.5;
        double x0 = 0.0;
        double x1 = 1.0;
        double y0 = 0.0;
        double y1 = 10.0;

        double y = linearly_interpolate(x, x0, x1, y0, y1);

        std::cout << "Interpolated value: " << y << '\n';
        std::cout << "Expected value:     " << 5.0 << "\n\n";

        FAIL_IF_NOT(std::abs(y - 5.0) < tolerance);
    }

    {
        std::cout << "Test 2: Lagrange interpolation of f(x) = x^2\n";

        std::array<double, 3> x_data{0.0, 1.0, 2.0};
        std::array<double, 3> y_data{0.0, 1.0, 4.0};

        double x = 1.5;
        double y = interpolate_lagrange<3>(x, x_data, y_data);

        std::cout << "Interpolated value: " << y << '\n';
        std::cout << "Expected value:     " << x * x << "\n\n";

        FAIL_IF_NOT(std::abs(y - x * x) < tolerance);
    }

    {
        std::cout << "Test 3: interpolate_from_vector with f(x) = x^2\n";

        std::vector<double> x_data{0.0, 1.0, 2.0, 3.0};
        std::vector<double> y_data{0.0, 1.0, 4.0, 9.0};

        double x = 1.5;

        // Use 3 interpolation points starting at index 0:
        // x = 0, 1, 2
        double y = interpolate_from_vector<3>(x, x_data, y_data, 0);

        std::cout << "Interpolated value: " << y << '\n';
        std::cout << "Expected value:     " << x * x << "\n\n";

        FAIL_IF_NOT(std::abs(y - x * x) < tolerance);
    }

    {
        std::cout << "Test 4: interpolate_from_vector near end of vector\n";

        std::vector<double> x_data{0.0, 1.0, 2.0, 3.0};
        std::vector<double> y_data{0.0, 1.0, 4.0, 9.0};

        double x = 2.5;

        // start_index_x is intentionally too close to the end.
        // The function should adjust it internally.
        double y = interpolate_from_vector<3>(x, x_data, y_data, 3);

        std::cout << "Interpolated value: " << y << '\n';
        std::cout << "Expected value:     " << x * x << "\n\n";

        FAIL_IF_NOT(std::abs(y - x * x) < tolerance);
    }

    {
        std::cout << "Test 5: interpolate_from_vector with index_offset_y\n";

        std::vector<double> x_data{0.0, 1.0, 2.0, 3.0};

        // First entry is dummy data.
        // Actual y-values start at index 1.
        std::vector<double> y_data{-123.0, 0.0, 1.0, 4.0, 9.0};

        double x = 1.5;

        double y = interpolate_from_vector<3>(x, x_data, y_data, 0, 1);

        std::cout << "Interpolated value: " << y << '\n';
        std::cout << "Expected value:     " << x * x << "\n\n";

        FAIL_IF_NOT(std::abs(y - x * x) < tolerance);
    }

    {
        std::cout << "Test 6: constexpr linear interpolation\n";

        constexpr double y = linearly_interpolate(
            0.25,
            0.0,
            1.0,
            0.0,
            8.0
        );

        static_assert(y == 2.0, "constexpr linear interpolation failed");

        std::cout << "constexpr interpolated value: " << y << '\n';
        std::cout << "Expected value:               " << 2.0 << "\n\n";
    }

    {
        std::cout << "Test 7: constexpr Lagrange interpolation\n";

        constexpr std::array<double, 3> x_data{0.0, 1.0, 2.0};
        constexpr std::array<double, 3> y_data{1.0, 2.0, 5.0};

        // This corresponds to f(x) = x^2 + 1
        constexpr double x = 1.5;
        constexpr double y = interpolate_lagrange<3>(x, x_data, y_data);

        static_assert(y == 3.25, "constexpr Lagrange interpolation failed");

        std::cout << "constexpr interpolated value: " << y << '\n';
        std::cout << "Expected value:               " << 3.25 << "\n\n";
    }

    std::cout << "All interpolation tests passed.\n";

    return 0;
}