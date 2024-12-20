#define _USE_MATH_DEFINES
#include <cmath>
#include "../include/Utility/Numerics/Integration/CauchyPrincipalValue.hpp"
#include <array>
#include <iostream>

constexpr size_t N = 4;
typedef double Real;
typedef std::array<Real, N> result_set;
constexpr Real EPS = 1e-10;

Real printer(Real ana, Real num) {
    std::cout << "Analytical: " << ana << " Numerical: " << num << " Error: " << std::abs(ana - num) << std::endl;
    return std::abs(ana - num);
}

result_set get_numerics(Real a, Real b, Real c) {
    auto f1 = [](Real x) { return 1.; }; // 1/(x-c)
    auto f2 = [](Real x) { return 1. / (x*x); }; // 1 / (x^2 * (x - c))
    auto f3 = [](Real x) { return std::exp(-x); }; // e^(-x) / (x-c)
    auto f4 = [](Real x) { return std::exp(x); }; // e^(x) / (x-c)
   
    static auto integrator = Utility::Numerics::Integration::CauchyPrincipalValue<Real, 60>();

    return {
        integrator.cauchy_principal_value(f1, a, b, c),
        integrator.cauchy_principal_value(f2, a, b, c),
        integrator.cauchy_principal_value(f3, a, b, c),
        integrator.cauchy_principal_value(f4, a, b, c)
    };
}

int main() {
    Real c = 0;
    Real a = -1;
    Real b = 1;

    result_set analytical_results = {
        0., 
        0., 
        -2.114501750751457029143684709, 
        2.114501750751457029143684709
    };
    result_set numerical_results = get_numerics(a, b, c);
    for (size_t i = 0U; i < N; ++i)
    {
        const Real error = printer(analytical_results[i], numerical_results[i]);
        if (error > EPS) {
            return 1;
        }
    }
    std::cout << "[a = -1, b = 1, c = 0] was successful!" << std::endl;


    c = 2;
    a = 1;
    b = 4;
    analytical_results = {
        std::log(2.0), 
        -(3./8.) - 0.25 * std::log(2.0), 
        -0.263094270910372120986974717345145551, 
        38.22815578220014947090427021723432
    };
    numerical_results = get_numerics(a, b, c);
    for (size_t i = 0U; i < N; ++i)
    {
        const Real error = printer(analytical_results[i], numerical_results[i]);
        if (error > EPS) {
            return 1;
        }
    }
    std::cout << "[a = 1, b = 4, c = 2] was successful!" << std::endl;


    c = -sqrt(80);
    a = -M_PI*M_PI;
    b = -M_1_PI;
    analytical_results = {
        2.23237865191957994812012549140384, 
        (b * std::log(std::abs((b - c) / b)) + c) / (c * c * b) - (a * std::log(std::abs((a - c) / a)) + c) / (c * c * a), 
        -12967.05645719712729238311198086643, 
        0.098064022478553059271310184995731
    };
    numerical_results = get_numerics(a, b, c);
    for (size_t i = 0U; i < N; ++i)
    {
        const Real error = printer(analytical_results[i], numerical_results[i]);
        if (error > EPS) {
            return 1;
        }
    }
    std::cout << "[a = -pi^2, b = -1/pi, c = -sqrt(80)] was successful!" << std::endl;
    return 0;
}