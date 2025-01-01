#pragma once
#include <array>
#include <limits>
#include <boost/math/quadrature/gauss.hpp>
#include <type_traits>

namespace mrock::utility::Numerics::Integration {
    /* Computes the cauchy principal value of f(x) / (x - singularity) over the interval [lower, upper].
    *  Based on the algorithm proposed in https://link.springer.com/article/10.1007/BF01935567
    *  Uses the Gauss-Legendre implementation from boost as a basis
    */
    template <class Real, int PolynomialDegree>
    class CauchyPrincipalValue {
    protected:
        static_assert(PolynomialDegree % 2 == 0, "Polynomial degree must be even");
        using gauss = boost::math::quadrature::gauss<Real, PolynomialDegree>;
        // Computes the cauchy principal value for int_-1^1 f(x) / x dx
        template <class F>
        static decltype(std::declval<F>()(Real{})) __cpv(F const& f, Real const& symmetrized_bound, Real const& singularity) {
            auto proper_integrand = [&f, &symmetrized_bound, &singularity](const Real& x) {
                assert(x >= std::numeric_limits<Real>::epsilon());
                assert(x < Real{1});
                return (f(symmetrized_bound * x + singularity) - f(singularity - symmetrized_bound * x ));
            };
            Real result{};
            for(int n = 0; n < PolynomialDegree / 2; ++n) {
                // We only need to consider the positive abscissas and boosts only saves the positive ones
                result += (gauss::weights()[n] / gauss::abscissa()[n]) * proper_integrand(gauss::abscissa()[n]);
            }
            // No additional contribution from symmetrized_bound required as we get one 
            // contribution from the Jacobi determinant and one contribution from the 1/x
            return result;
        }
        // Computes an integral of f(x + c) / x for and interval that does not include 0 
        template <class F>
        static decltype(std::declval<F>()(Real{})) __regular(F const& f, Real lower, Real upper, Real singularity) {
            assert(lower * upper > 0);
            auto proper_integrand = [&f, &singularity](const Real& x) {
                return f(x + singularity) / x;
            };
            return gauss::integrate(proper_integrand, lower, upper);
        }
    public:
        template <class F>
        static decltype(std::declval<F>()(Real{})) cauchy_principal_value(F const& f, Real lower, Real upper, Real singularity) {
            // assert that the lower edge is _actually_ the lower edge
            if (lower > upper) {
                return -cauchy_principal_value(f, upper, lower, singularity);
            }
            // If there is no singulartiy, we can just use the regular integration
            if (singularity < lower || singularity > upper) {
                return __regular(f, lower - singularity, upper - singularity, singularity);
            }

            // Shift singularity to 0
            // x -> x + singularity
            lower -= singularity;
            upper -= singularity;
            if(std::abs(lower) < upper) {
                const Real middle_bound = std::abs(lower);
                const Real cpv_integral = __cpv(f, middle_bound, singularity);
                const Real regular_integral = __regular(f, middle_bound, upper, singularity);
                return cpv_integral + regular_integral;
            }
            else {
                const Real middle_bound = upper;
                const Real cpv_integral = __cpv(f, middle_bound, singularity);
                const Real regular_integral = __regular(f, lower, -middle_bound, singularity);
                return cpv_integral + regular_integral;
            }
        }
    };
}