/**
 * @file AdaptiveTrapezoidalRule.hpp
 * @brief Provides numerical integration using the adaptive trapezoidal rule.
 */

#pragma once

#include "../../ThrowException.hpp"
#include "../ErrorFunctors.hpp"

#include <cmath>
#include <type_traits>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <functional>

namespace mrock::utility::Numerics::Integration {
	struct adapative_trapezoidal_rule_print_policy {
		bool print_values;
		bool print_error;
	};
	constexpr adapative_trapezoidal_rule_print_policy adapative_trapezoidal_rule_print_nothing = {false, false};

	/**
	 * @struct adapative_trapezoidal_rule
	 * @brief Performs numerical integration using the adaptive trapezoidal rule.
	 * 
	 * @tparam RealType The type of the real numbers (default is double).
	 * @tparam print_steps Whether to print the steps of the integration process (default is false).
	 */
	template<class RealType = double, adapative_trapezoidal_rule_print_policy print_steps = adapative_trapezoidal_rule_print_nothing>
	struct adapative_trapezoidal_rule {
	private:
		/**
		 * @brief Type alias for the result of the unary function with const, volatile, and reference removed.
		 * 
		 * @tparam UnaryFunction The type of the unary function.
		 */
		template<class UnaryFunction>
		using decayed_result = std::remove_cvref_t< std::invoke_result_t< UnaryFunction, RealType > >;

		/**
		 * @brief Type alias for the result of the error function.
		 * 
		 * @tparam UnaryFunction The type of the unary function.
		 * @tparam ErrorFunction The type of the error function.
		 */
		template<class UnaryFunction, class ErrorFunction>
		using error_result = std::invoke_result_t< ErrorFunction, decayed_result<UnaryFunction>, decayed_result<UnaryFunction> >;
	public:
		/**
		 * @brief Integrates a function over a given interval using the adaptive trapezoidal rule.
		 * 
		 * @tparam UnaryFunction The type of the function to integrate.
		 * @tparam ErrorFunction The type of the error function (default is scalar_error).
		 * @param function The function to integrate.
		 * @param begin The beginning of the interval.
		 * @param end The end of the interval.
		 * @param num_steps The initial number of steps.
		 * @param max_error The maximum allowable error.
		 * @param error_func The error function (default is ErrorFunction()).
		 * @param zero The zero value for the result type (default is the default-constructed value).
		 * @return The integral of the function over the interval.
		 */
		template <class UnaryFunction, class ErrorFunction = scalar_error<decayed_result<UnaryFunction>>>
		decayed_result<UnaryFunction> integrate(const UnaryFunction& function, const RealType begin, const RealType end, unsigned int num_steps,  
			const error_result<UnaryFunction, ErrorFunction> max_error, const ErrorFunction& error_func = ErrorFunction(),
			const decayed_result<UnaryFunction>& zero = decayed_result<UnaryFunction>{}) 
		{
			throw_exception<std::domain_error>(MROCK_NDEBUG_CONDITION(!(std::isfinite(begin) && std::isfinite(end))), "The integration domain is not sensible!");
			using std::abs;
			if (abs(begin - end) < std::numeric_limits<RealType>::epsilon()) return zero;
			if (begin < end) {
				return -integrate(function, end, begin, num_steps, max_error, error_func, zero);
			}
			using result_type = decayed_result<UnaryFunction>;
			RealType step = (end - begin) / num_steps;

			result_type new_value{};
			result_type old_value{0.5 * (function(begin) + function(end))};
			for (unsigned int n = 1U; n < num_steps; ++n) {
				old_value += function(begin + n * step);
			}
			old_value *= step;
			error_result<UnaryFunction, ErrorFunction> current_error{};

			do {
				step *= 0.5;
				num_steps *= 2;

				new_value = zero;
				for (unsigned int n = 1U; n < num_steps; n+=2) {
					new_value += function(begin + n * step);
				}
				new_value = 0.5 * old_value + step * new_value;

				current_error = error_func(new_value, old_value);

				if constexpr (print_steps.print_values) {
					std::cout << "I_n = " << old_value << "\nI_(n+1) = " << new_value;
				}
				if constexpr (print_steps.print_error) {
					std::cout << "\n\t||\terror = " << current_error << "\tCurrent step = " << step << "\tCurrent number of steps = " << num_steps << std::endl;
				}

				old_value = new_value;
			} while (current_error > max_error);
			return new_value;
		}
	};
}