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

// Returns the relative error if ana is not 0 otherwise it returns the absolute error
Real printer(Real ana, Real num, std::string const& algorithm) {
    std::cout << algorithm << ": \t" << "Analytical: " << ana << " Numerical: " << num << " Error: " << std::abs(ana - num) << std::endl;
    if (std::abs(ana) > 1e-12)
        return std::abs(ana - num) / std::abs(ana);
    return std::abs(ana - num);
}

template<int N=60>
result_set_A get_numerics(Real a, Real b, Real c) {
    auto f1 = [](Real x) { return 1.; }; // 1/(x-c)
    auto f2 = [](Real x) { return 1. / (x*x); }; // 1 / (x^2 * (x - c))
    auto f3 = [](Real x) { return std::exp(-x); }; // e^(-x) / (x - c)
    auto f4 = [](Real x) { return std::exp(x); }; // e^(x) / (x - c)

    return {
        mrock::utility::Numerics::Integration::CauchyPrincipalValue<Real, N>::cauchy_principal_value(f1, a, b, c),
        mrock::utility::Numerics::Integration::CauchyPrincipalValue<Real, N>::cauchy_principal_value(f2, a, b, c),
        mrock::utility::Numerics::Integration::CauchyPrincipalValue<Real, N>::cauchy_principal_value(f3, a, b, c),
        mrock::utility::Numerics::Integration::CauchyPrincipalValue<Real, N>::cauchy_principal_value(f4, a, b, c)
    };
}

template<int N=60>
result_set_A get_generalized_numerics(Real a, Real b, Real c) {
    auto f1 = [&c](Real x) { return 1. / (x - c); };
    auto f2 = [&c](Real x) { return 1. / ((x - c) * x * x); }; 
    auto f3 = [&c](Real x) { return std::exp(-x) / (x - c); };
    auto f4 = [&c](Real x) { return std::exp(x) / (x - c); }; 

    const std::vector<Real> second_singularitites = (c != 0 ? std::vector<Real>{std::min(0.0, c), std::max(c, 0.0)} : std::vector<Real>{ 0.0 });
    return {
        mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<Real, N>::generalized_principal_value(f1, a, b, std::vector<Real>{ c }),
        mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<Real, N>::generalized_principal_value(f2, a, b, second_singularitites),
        mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<Real, N>::generalized_principal_value(f3, a, b, std::vector<Real>{ c }),
        mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<Real, N>::generalized_principal_value(f4, a, b, std::vector<Real>{ c })
    };
}

result_set_A get_analytics(Real a, Real b, Real c) {
    return {
           std::log(std::abs(c - b) / std::abs(c - a)),
            (
                c != Real{0} ? 
                    std::log(std::abs(c - b) / std::abs(c - a)) / (c*c)
                    - std::log(std::abs(b / a)) / (c*c)
                    + (1. / b - 1. / a) / c
                : 0.5 * (1. / (a*a) - 1. / (b*b))
            ),
            (std::expint(c - b) - std::expint(c - a)) * std::exp(-c),
            (std::expint(b - c) - std::expint(a - c)) * std::exp(c)
        };
}

int main() {
    const std::string cpv_name = "CauchyPrincipalValue";
    const std::string gpv_name = "GeneralizedPrincipalValue";
    {
        Real c = 0;
        Real a = -1;
        Real b = 1;

        result_set_A analytical_results = get_analytics(a, b, c);
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
        analytical_results = get_analytics(a, b, c);
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
        analytical_results = get_analytics(a, b, c);
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
    std::cout << "\n###############################\n" << std::endl;
    {
        std::cout << "Stress test: Singularity is very close one of the limits of integration" << std::endl;

        Real a = -1;
        Real b = 0.1;
        Real c = 0;
        
        result_set_A analytical_results = get_analytics(a, b, c);
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
        std::cout << "[a = -1, b = 0.1, c = 0] was successful!" << std::endl;

        a = -1e-4;
        b = 1e-4;
        c = 0;
        
        analytical_results = get_analytics(a, b, c);
        numerical_results = get_numerics<100>(a, b, c);
        generalized_numerical_results = get_generalized_numerics<100>(a, b, c);
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
        std::cout << "[a = -1e-4, b = 1e-4, c = 0] was successful!" << std::endl;

        /* 
        * This is a truly difficult case, but it seems that the issue is actually the normal integration and not the principal value
        * because setting the integration limits to similar magnitude (i.e. computing mainly the PV) yields very good results
        * even for small polynomial degrees N
        */
        a = -1e-8;
        b = 1e-2;
        c = 0;
        
        analytical_results = get_analytics(a, b, c);
        numerical_results = get_numerics<10000>(a, b, c);
        generalized_numerical_results = get_generalized_numerics<10000>(a, b, c);
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
        std::cout << "[a = -1e9, b = 1e-2, c = 0] was successful!" << std::endl;
    }

    return 0;
}