#include <Eigen/Dense>
#include <mrock/utility/Numerics/Roots/BroydensMethod.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>

#define FAIL_IF_NOT(cond)                                                      \
    if (!(cond)) {                                                             \
        std::cerr << "FAILED: " #cond << " at line " << __LINE__ << std::endl; \
        return 1;                                                              \
    }

int main() {
    using namespace mrock::utility::Numerics::Roots;

    {
        std::cout << "Test 1: scalar nonlinear equation x^2 - 2 = 0\n";

        BroydensMethod<double, 1> solver;

        Eigen::Matrix<double, 1, 1> x;
        x(0) = 1.0;

        auto f = [](const Eigen::Matrix<double, 1, 1>& input, Eigen::Matrix<double, 1, 1>& output) {
            output(0) = input(0) * input(0) - 2.0;
        };

        double residual = solver.compute(f, x, 100);

        std::cout << "Computed root: " << x(0) << '\n';
        std::cout << "Expected root: " << std::sqrt(2.0) << '\n';
        std::cout << "Residual: " << residual << "\n\n";

        FAIL_IF_NOT(residual < 1e-10);
        FAIL_IF_NOT(std::abs(x(0) - std::sqrt(2.0)) < 1e-8);
    }

    {
        std::cout << "Test 2: real 2D linear system\n";

        BroydensMethod<double, 2> solver;

        Eigen::Vector2d x;
        x << 0.0, 0.0;

        Eigen::Vector2d target;
        target << 1.0, -2.0;

        auto f = [&target](const Eigen::Vector2d& input, Eigen::Vector2d& output) {
            // Root is input == target
            output = input - target;
        };

        double residual = solver.compute(f, x, 50);

        std::cout << "Computed solution:\n" << x << '\n';
        std::cout << "Expected solution:\n" << target << '\n';
        std::cout << "Residual: " << residual << "\n\n";

        FAIL_IF_NOT(residual < 1e-12);
        FAIL_IF_NOT((x - target).norm() < 1e-12);
    }

    {
        std::cout << "Test 3: complex dynamic-size equation\n";

        using Complex = std::complex<double>;

        BroydensMethod<Complex, Eigen::Dynamic> solver;

        Eigen::VectorXcd z(1);
        z(0) = Complex{0.0, 0.0};

        Eigen::VectorXcd target(1);
        target(0) = Complex{1.0, -2.0};

        auto f = [&target](const Eigen::VectorXcd& input, Eigen::VectorXcd& output) {
            // Root is input == target
            output = input - target;
        };

        double residual = solver.compute(f, z, 50);

        std::cout << "Computed complex root: " << z(0) << '\n';
        std::cout << "Expected complex root: " << target(0) << '\n';
        std::cout << "Residual: " << residual << "\n\n";

        FAIL_IF_NOT(residual < 1e-12);
        FAIL_IF_NOT(std::abs(z(0) - target(0)) < 1e-12);
    }

    std::cout << "All Broyden tests passed.\n";

    return 0;
}