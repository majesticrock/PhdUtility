#define _USE_MATH_DEFINES
#include <cmath>
#include "../include/mrock/utility/Numerics/Integration/CauchyPrincipalValue.hpp"
#include "../include/mrock/utility/Numerics/Integration/GeneralizedPrincipalValue.hpp"
#include <array>
#include <iostream>

constexpr size_t N_A = 4;
typedef double Real;
typedef std::array<Real, N_A> result_set_A;
constexpr size_t N_B = 4;
typedef std::array<Real, N_B> result_set_B;
constexpr Real EPS = 1e-10;

using base_integrator = mrock::utility::Numerics::Integration::CauchyPrincipalValue<Real, 60>;
using gen_integrator = mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<Real, 60>;

Real printer(Real ana, Real num, std::string const& algorithm) {
    std::cout << algorithm << ": \t" << "Analytical: " << ana << " Numerical: " << num << " Error: " << std::abs(ana - num) << std::endl;
    return std::abs(ana - num);
}

result_set_A get_numerics(Real a, Real b, Real c) {
    auto f1 = [](Real x) { return 1.; }; // 1/(x-c)
    auto f2 = [](Real x) { return 1. / (x*x); }; // 1 / (x^2 * (x - c))
    auto f3 = [](Real x) { return std::exp(-x); }; // e^(-x) / (x - c)
    auto f4 = [](Real x) { return std::exp(x); }; // e^(x) / (x - c)

    return {
        base_integrator::cauchy_principal_value(f1, a, b, c),
        base_integrator::cauchy_principal_value(f2, a, b, c),
        base_integrator::cauchy_principal_value(f3, a, b, c),
        base_integrator::cauchy_principal_value(f4, a, b, c)
    };
}

result_set_A get_generalized_numerics(Real a, Real b, Real c) {
    auto f1 = [&c](Real x) { return 1. / (x - c); };
    auto f2 = [&c](Real x) { return 1. / ((x - c) * x * x); }; 
    auto f3 = [&c](Real x) { return std::exp(-x) / (x - c); };
    auto f4 = [&c](Real x) { return std::exp(x) / (x - c); }; 

    const std::vector<Real> second_singularitites = (c != 0 ? std::vector<Real>{std::min(0.0, c), std::max(c, 0.0)} : std::vector<Real>{ 0.0 });
    return {
        gen_integrator::generalized_principal_value(f1, a, b, std::vector<Real>{ c }),
        gen_integrator::generalized_principal_value(f2, a, b, second_singularitites),
        gen_integrator::generalized_principal_value(f3, a, b, std::vector<Real>{ c }),
        gen_integrator::generalized_principal_value(f4, a, b, std::vector<Real>{ c })
    };
}

int main() {
    const std::string cpv_name = "CauchyPrincipalValue";
    const std::string gpv_name = "GeneralizedPrincipalValue";
    {
        Real c = 0;
        Real a = -1;
        Real b = 1;

        result_set_A analytical_results = {
            0., 
            0., 
            -2.114501750751457029143684709, 
            2.114501750751457029143684709
        };
        result_set_A numerical_results = get_numerics(a, b, c);
        result_set_A generalized_numerical_results = get_generalized_numerics(a, b, c);
        for (size_t i = 0U; i < N_A; ++i)
        {
            const Real error = printer(analytical_results[i], numerical_results[i], cpv_name);
            if (error > EPS) {
                return 1;
            }
            const Real generalized_error = printer(analytical_results[i], generalized_numerical_results[i], gpv_name);
            if (generalized_error > EPS) {
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
        generalized_numerical_results = get_generalized_numerics(a, b, c);
        for (size_t i = 0U; i < N_A; ++i)
        {
            const Real error = printer(analytical_results[i], numerical_results[i], cpv_name);
            if (error > EPS) {
                return 1;
            }
            const Real generalized_error = printer(analytical_results[i], generalized_numerical_results[i], gpv_name);
            if (generalized_error > EPS) {
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
        generalized_numerical_results = get_generalized_numerics(a, b, c);
        for (size_t i = 0U; i < N_A; ++i)
        {
            const Real error = printer(analytical_results[i], numerical_results[i], cpv_name);
            if (error > EPS) {
                return 1;
            }
            const Real generalized_error = printer(analytical_results[i], generalized_numerical_results[i], gpv_name);
            if (generalized_error > EPS) {
                return 1;
            }
        }
        std::cout << "[a = -pi^2, b = -1/pi, c = -sqrt(80)] was successful!" << std::endl;
    }
    std::cout << "\n###############################\n" << std::endl;
    {
        Real a = -2;
        Real b = 2;
        Real c = 1;
        Real d = 1.5;

        auto f1 = [&c](Real x) { return 1. / (x * x - c * c); };
        auto f2 = [&c, &d](Real x) { return 1. / ((x - c) * (x - d)); };
        auto f3 = [](Real x) { return 1. / std::sin(x); };
        auto f4 = [&c](Real x) { return 1. / std::log(std::abs(2*x - c)); };

        auto ana1 = [](Real a, Real b, Real c) {
            return ( std::log( std::abs( (a + c) / (a - c) )) 
                   - std::log( std::abs( (b + c) / (b - c) )) ) / (2 * c);
        };
        auto ana2 = [](Real a, Real b, Real c, Real d) {
            return ( std::log( std::abs( (b - d) / (b - c) )) 
                   - std::log( std::abs( (a - d) / (a - c) )) ) / (d - c);
        };
        auto ana3 = [](Real a, Real b) {
            return std::log(std::abs(std::tan(b / 2) / std::tan(a / 2)));
        };

        result_set_A numerical_results = { 
            gen_integrator::generalized_principal_value(f1, a, b, std::vector<Real>{ -c, c }),
            gen_integrator::generalized_principal_value(f2, a, b, std::vector<Real>{ c, d }),
            gen_integrator::generalized_principal_value(f3, a, b, std::vector<Real>{ -2 * M_PI, -M_PI, 0, M_PI, 2 * M_PI }),
            // Computing these kinds of logarithmic singularities is incredibly difficult, but possible
            mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<Real, 10000>::generalized_principal_value(f4, a, b, std::vector<Real>{ 0.5 * (c - 1), 0.5 * (c + 1) })
        };
        result_set_A analytical_results = { 
            ana1(a, b, c),
            ana2(a, b, c, d),
            ana3(a, b),
            2.899088452349921918634632537104266318956263881341842
        };

        for (size_t i = 0U; i < N_B; ++i)
        {
            const Real error = printer(analytical_results[i], numerical_results[i], gpv_name);
            if (error > EPS) {
                return 1;
            }
        }
        std::cout << "[a = -2, b = 2, c = 1, d=1.5] was successful!" << std::endl;

        a = -11;
        b = 3;
        c = 2.123456;
        d = -1.123456;
        numerical_results = { 
            gen_integrator::generalized_principal_value(f1, a, b, std::vector<Real>{ -c, c }),
            gen_integrator::generalized_principal_value(f2, a, b, std::vector<Real>{ d, c }),
            gen_integrator::generalized_principal_value(f3, a, b, std::vector<Real>{ -4 * M_PI, -3 * M_PI, -2 * M_PI, -M_PI, 0, M_PI, 2 * M_PI, 3 * M_PI, 4 * M_PI }),
            // Computing these kinds of logarithmic singularities is incredibly difficult, but possible
            mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<Real, 10000>::generalized_principal_value(f4, a, b, std::vector<Real>{ 0.5 * (c - 1), 0.5 * (c + 1) })
        };
        analytical_results = { 
            ana1(a, b, c),
            ana2(a, b, c, d),
            ana3(a, b),
            7.058322126020065070227595524877212606453
        };

        for (size_t i = 0U; i < N_B; ++i)
        {
            const Real error = printer(analytical_results[i], numerical_results[i], gpv_name);
            if (error > EPS) {
                return 1;
            }
        }
        std::cout << "[a = -11, b = 3, c = 2.123456, d=-1.123456] was successful!" << std::endl;
    }

    return 0;
}