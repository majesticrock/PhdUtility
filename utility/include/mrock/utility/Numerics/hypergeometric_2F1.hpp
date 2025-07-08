#pragma once
#include <type_traits>
#include <cmath>
#include <stdexcept>
#include <limits>

#include "../ThrowException.hpp"
#include "is_integer.hpp"
#include "../IsComplex.hpp"

#include <iostream>

namespace mrock::utility::Numerics {
    namespace detail {
        template<class T>
        struct KahanSum {
            T sum = 0.0;
            T c = 0.0;  // Compensation for lost low-order bits

            KahanSum& operator+=(T value) {
                T y = value - c;
                T t = sum + y;
                c = (t - sum) - y;
                sum = t;
                return *this;
            }

            T result() const { return sum; }
        };

        /**
         * Analytical continuation as a power series from https://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/02/02/
         * 
         * The maximum reachable precision is currently around 1e-8.
         */
        template<Real real, Complex z_t>
        auto hg_2F1_continuation(real a, real b, real c, z_t z, real const& tol)
        {
            using result_t = std::common_type_t<real, z_t>;
            using std::tgamma;
            using std::pow;
            using std::abs;

            // If a-b is an integer, we need to work in a limit
            if (is_integer(a - b, tol)) {
                b += real{2} * tol;
            }
            // If |z| = 1, we need to work in a limit as well
            if (abs(z) - real{1} < tol) {
                z += tol;
            }

            // The rising Pochhammer symbol is defined as (x)_n = Gamma(x + n) / (Gamma(x)) ~> (x+n)! / x! ~> x(x+1)(x+2)...(x+n-1)(x+n)
            // It can thus be computed as (x)_n+1 = (x + n)(x)_n

            result_t summand = a * (a - c + real{1}) / ((a - b + real{1}) * z);
            KahanSum<result_t> first_term{real{1} + summand, result_t{}};
            real n{1};
            while( abs(summand) > tol ) {
                ++n;
                summand *= (a + n - real{1}) * (a - c + n) / (n * (a - b + n) * z);
                first_term += summand;
            }

            n = real{1};
            summand = b * (b - c + real{1}) / ((b - a + real{1}) * z);
            KahanSum<result_t> second_term{real{1} + summand, result_t{}};
            while( abs(summand) > tol) {
                ++n;
                summand *= (b + n - real{1}) * (b - c + n) / (n * (b - a + n) * z);
                second_term += summand;
            }

            // To correctly respect the branch cut of the complex square root, we need to avoid -0 (which exists in c++)
            const z_t negative_z = {-std::real(z), std::imag(z) == typename z_t::value_type{} ? typename z_t::value_type{} : -std::imag(z)};
            const result_t prefactor_1 = tgamma(b - a) * tgamma(c) * pow(negative_z, -a) / (tgamma(b) * tgamma(c - a));
            const result_t prefactor_2 = tgamma(a - b) * tgamma(c) * pow(negative_z, -b) / (tgamma(a) * tgamma(c - b));
            return prefactor_1 * first_term.result() + prefactor_2 * second_term.result();
        }

        template<Real real, class z_t>
        auto hg_2F1(real a, real b, real c, z_t z, real const& tol) 
        {
            using result_t = std::common_type_t<real, z_t>;
            using std::abs;

            result_t summand = a * b * z / c;
            result_t result = real{1} + summand; // (x)_n = 1
            real n{1};
            while ( abs(summand) > tol ) {
                ++a, ++b, ++c, ++n;
                summand *= a * b * z / (n * c);
                result += summand;
            }
            return result;
        }
    }

    /**
     * The maximum reachable precision for the analytic continuation (i.e. z) is currently around 1e-8.
     * Asking for more precise results actually *decreases* the precision due to floating point arithmetic.
     * If higher precisions are needed use long double. If that doesn't help, I need to split the summation as to avoid adding 1 + 1e-10 a million times and getting bad results.
     */
    template<Real real, Complex z_t>
    auto hypergeometric_2F1(real a, real b, real c, z_t z, real const& tol = real{1e-10})
    {
        using std::abs;

        throw_exception<std::invalid_argument>((is_integer(c) && c < real{0}) || abs(c) < tol, "c must neither be a negative integer nor 0!");

        if (abs(z) < real{1}) {
            return detail::hg_2F1(a, b, c, z, tol);
        }
        return detail::hg_2F1_continuation(a, b, c, z, tol);
    }

    /**
     * The maximum reachable precision for the analytic continuation (i.e. z) is currently around 1e-8.
     * Asking for more precise results actually *decreases* the precision due to floating point arithmetic.
     * If higher precisions are needed use long double. If that doesn't help, I need to split the summation as to avoid adding 1 + 1e-10 a million times and getting bad results.
     */
    template<Real real, Real z_t>
    auto hypergeometric_2F1(real a, real b, real c, z_t z, real const& tol = real{1e-10})
    {
        return hypergeometric_2F1(a, b, c, std::complex<z_t>{z, z_t{}}, tol);
    }
}