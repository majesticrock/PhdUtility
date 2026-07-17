#include <mrock/utility/Numerics/Roots/Bisection.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <typeinfo>

int main() {
    using mrock::utility::Numerics::Roots::bisection;
    using mrock::utility::Numerics::Roots::NoRootException;

    {
        // Test 1: Root of x^2 - 2 is sqrt(2)
        auto f = [](double x) { return x * x - 2.0; };

        double root = bisection(f, 0.0, 2.0, 1e-12, 100);

        std::cout << "Root of x^2 - 2: " << root << '\n';
        std::cout << "Expected: " << std::sqrt(2.0) << '\n';

        assert(std::abs(root - std::sqrt(2.0)) < 1e-10);
    }

    {
        // Test 2: Root of sin(x) in interval [3, 4] is pi
        auto f = [](double x) { return std::sin(x); };

        double root = bisection(f, 3.0, 4.0, 1e-12, 100);

        std::cout << "Root of sin(x): " << root << '\n';
        std::cout << "Expected: " << M_PI << '\n';

        assert(std::abs(root - M_PI) < 1e-10);
    }

    {
        // Test 3: Root at interval boundary
        auto f = [](double x) { return x - 5.0; };

        double root = bisection(f, 5.0, 10.0, 1e-12, 100);

        std::cout << "Root at boundary: " << root << '\n';

        assert(std::abs(root - 5.0) < 1e-12);
    }

    {
        // Test 4: No root in interval should throw
        auto f = [](double x) { return x * x + 1.0; };

        try {
            double root = bisection(f, -1.0, 1.0, 1e-12, 100);
            std::cerr << "Unexpected root found: " << root << '\n';
            assert(false);
        } catch (const NoRootException& e) {
            std::cout << "Caught expected NoRootException:\n";
            std::cout << e.what() << '\n';
        }
    }

    std::cout << "All tests passed.\n";

    return 0;
}