#pragma once
#include "CauchyPrincipalValue.hpp"
#include <iostream>
#include <algorithm>

namespace Utility::Numerics::Integration {
    template <class Real, int PolynomialDegree>
    class GeneralizedPrincipalValue {
    private:
        static_assert(PolynomialDegree % 2 == 0, "Polynomial degree must be even");
        using gauss = boost::math::quadrature::gauss<Real, PolynomialDegree>;
        using cauchy = Utility::Numerics::Integration::CauchyPrincipalValue<Real, PolynomialDegree>;
        template<class F, class G> using result_type = std::common_type_t<decltype(std::declval<F>()(Real{})), decltype(std::declval<G>()(Real{}))>;

        template <class F>
        static decltype(std::declval<F>()(Real{})) __pv(F const& f, Real const& symmetrized_bound, Real const& singularity) {
            auto proper_integrand = [&f, &symmetrized_bound, &singularity](const Real& x) {
                assert(x >= std::numeric_limits<Real>::epsilon());
                assert(x < Real{1});
                return (f(symmetrized_bound * x + singularity) + f(singularity - symmetrized_bound * x ));
            };
            decltype(std::declval<F>()(Real{})) result{};
            for (int n = 0; n < PolynomialDegree / 2; ++n) {
                // We only need to consider the positive abscissas and boosts only saves the positive ones
                result += gauss::weights()[n] * proper_integrand(gauss::abscissa()[n]);
            }
            result *= symmetrized_bound;
            return result;
        }
    public:
        // Computes the generalized principal value of f(x) / g(x) over the interval [lower, upper]
        // Requires the zeros of g(x) and the derivative of g(x) at the zeros
        // Assumes that the zeros are ordered and that the first order Taylor expansion of g(x) is a good approximation near the zeros
        // taylor_tolerance is the distance from the zero where the Taylor expansion is used
        // RegularQuadrature is used for the regular part of the integration. It requires the syntax RegularQuadrature::integrate(f, lower, upper).
        template<class RegularQuadrature, class F, class G, class RandomAccessContainer>
        static result_type<F, G> principal_value_f_over_g(F const& f, G const& g, Real lower, Real upper, 
            RandomAccessContainer const& zeros, RandomAccessContainer const& g_prime_at_zeros, Real taylor_tolerance = 1e-3) 
        {
            if (lower > upper) {
                return -principal_value_f_over_g(f, g, upper, lower, zeros, g_prime_at_zeros, taylor_tolerance);
            }
            result_type<F, G> result{};

            auto integrand = [&f, &g](const Real& x) {
                return f(x) / g(x);
            };

            for(int i = 0; i < zeros.size(); ++i) {
                result += RegularQuadrature::integrate(integrand, lower, zeros[i] - taylor_tolerance);
                result += cauchy::cauchy_principal_value(f, zeros[i] - taylor_tolerance, zeros[i] + taylor_tolerance, zeros[i]) / g_prime_at_zeros[i];
                lower = zeros[i] + taylor_tolerance;
            }
            result += RegularQuadrature::integrate(integrand, lower, upper);

            return result;
        }

        // Computes the principal value of any f(x) with a singularity at singularity over the interval [lower, upper]
        template<class F>
        static decltype(std::declval<F>()(Real{})) generalized_principal_value(F const& f, Real lower, Real upper, Real singularity) {
            if (lower > upper) {
                return -generalized_principal_value(f, upper, lower, singularity);
            }
            if (lower > singularity || upper < singularity) {
                return gauss::integrate(f, lower, upper);
            }

            // Shift singularity to 0
            // x -> x + singularity
            lower -= singularity;
            upper -= singularity;
            if (std::abs(lower) < upper) {
                const Real middle_bound = std::abs(lower);
                const Real cpv_integral = __pv(f, middle_bound, singularity);
                const Real regular_integral = gauss::integrate(f, middle_bound + singularity, upper + singularity);
                return cpv_integral + regular_integral;
            }
            else {
                const Real middle_bound = upper;
                const Real cpv_integral = __pv(f, middle_bound, singularity);
                const Real regular_integral = gauss::integrate(f, lower + singularity, singularity - middle_bound);
                return cpv_integral + regular_integral;
            }
        }

        template<class F, class RandomAccessContainer>
        static decltype(std::declval<F>()(Real{})) generalized_principal_value(F const& f, Real lower, Real upper, RandomAccessContainer const& singularities) {
            int n = 0;
            while(singularities[n] < lower) {
                ++n;
            }
            Real halfway;
            decltype(std::declval<F>()(Real{})) result{};
            while(++n < singularities.size() && singularities[n] < upper) {
                assert((singularities[n] > singularities[n - 1]) && "Singularities must be ordered");
                halfway = 0.5 * (singularities[n - 1] + singularities[n]);
                result += generalized_principal_value(f, lower, halfway, singularities[n - 1]);
                lower = halfway;
            }
            result += generalized_principal_value(f, lower, upper, singularities[n - 1]);
            return result;
        }
    };
}