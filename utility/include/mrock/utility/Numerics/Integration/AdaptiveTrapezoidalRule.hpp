/**
 * @file AdaptiveTrapezoidalRule.hpp
 * @brief Provides numerical integration using the adaptive trapezoidal rule.
 */

#pragma once
#include <utility>
#include <cmath>
#include <type_traits>
#include <iostream>
#include <numeric>
#include <functional>

namespace mrock::utility::Numerics::Integration {
	/**
	 * @struct relative_error
	 * @brief Computes the relative error between two values.
	 * 
	 * @tparam T Type of the values.
	 */
	template<class T>
	struct relative_error {
		/**
		 * @brief Computes the relative error between two values.
		 * 
		 * @param _new The new value.
		 * @param _old The old value.
		 * @return double The relative error.
		 */
		double operator()(const T& _new, const T& _old) const {
			using std::abs;
			const double abs_new = abs(_new);
			return (abs_new > 1e-8 ?  abs(_new - _old) / abs_new : abs(_new - _old));
		};
	};

	/**
	 * @struct vector_norm_error
	 * @brief Computes the error based on the norm of vectors.
	 * 
	 * @tparam vector_type Type of the vector.
	 * @tparam value_type Type of the values in the vector.
	 */
	template<class vector_type, class value_type = typename vector_type::value_type>
	struct vector_norm_error {
	private:
		/**
		 * @brief Computes the norm of a vector.
		 * 
		 * @param vec The vector.
		 * @return value_type The norm of the vector.
		 */
		static value_type norm(const vector_type& vec) {
			return sqrt(std::transform_reduce(vec.begin(), vec.end(), vec.begin(), value_type{}));
		}
	public:
		/**
		 * @brief Computes the error based on the norm of two vectors.
		 * 
		 * @param _new The new vector.
		 * @param _old The old vector.
		 * @return value_type The error based on the norm.
		 */
		value_type operator()(const vector_type& _new, const vector_type& _old) const {
			const value_type abs_new = norm(_new);
			return (abs_new > value_type{1e-8} ?  norm(_new - _old) / abs_new : norm(_new - _old));
		}
	};

	/**
	 * @struct vector_elementwise_error
	 * @brief Computes the element-wise error between two vectors.
	 * 
	 * @tparam vector_type Type of the vector.
	 * @tparam value_type Type of the values in the vector.
	 */
	template<class vector_type, class value_type = typename vector_type::value_type>
	struct vector_elementwise_error {
		/**
		 * @brief Computes the element-wise error between two vectors.
		 * 
		 * @param _new The new vector.
		 * @param _old The old vector.
		 * @return value_type The element-wise error.
		 */
		value_type operator()(const vector_type& _new, const vector_type& _old) const {
			const value_type abs_new = std::transform_reduce(_new.begin(), _new.end(), value_type{}, 
										std::plus<value_type>(), [](const value_type& val){ using std::abs; return abs(val); });

			const value_type abs_diff = std::transform_reduce(_new.begin(), _new.end(), _old.begin(), value_type{}, 
										std::plus<value_type>(), [](const value_type& left, const value_type& right){ using std::abs; return abs(left - right); });
			return (abs_new > value_type{1e-8} ?  abs_diff / abs_new : abs_diff);
		}
	};

	/**
	 * @struct adapative_trapezoidal_rule
	 * @brief Performs numerical integration using the adaptive trapezoidal rule.
	 * 
	 * @tparam RealType The type of the real numbers (default is double).
	 * @tparam print_steps Whether to print the steps of the integration process (default is false).
	 */
	template<class RealType = double, bool print_steps = false>
	struct adapative_trapezoidal_rule {
	private:
		/**
		 * @brief Type alias for the result of the error function.
		 * 
		 * @tparam UnaryFunction The type of the unary function.
		 * @tparam ErrorFunction The type of the error function.
		 */
		template<class UnaryFunction, class ErrorFunction>
		using error_result = std::invoke_result_t<ErrorFunction, std::invoke_result_t<UnaryFunction, RealType>, std::invoke_result_t<UnaryFunction, RealType>>;

	public:
		/**
		 * @brief Integrates a function over a given interval using the adaptive trapezoidal rule.
		 * 
		 * @tparam UnaryFunction The type of the function to integrate.
		 * @tparam ErrorFunction The type of the error function (default is relative_error).
		 * @param function The function to integrate.
		 * @param begin The beginning of the interval.
		 * @param end The end of the interval.
		 * @param num_steps The initial number of steps.
		 * @param max_error The maximum allowable error.
		 * @param error_func The error function (default is ErrorFunction()).
		 * @param zero The zero value for the result type (default is the default-constructed value).
		 * @return The integral of the function over the interval.
		 */
		template <class UnaryFunction, class ErrorFunction = relative_error<std::invoke_result_t<UnaryFunction, RealType>>>
		auto integrate(const UnaryFunction& function, const RealType begin, const RealType end, unsigned int num_steps,  
			const error_result<UnaryFunction, ErrorFunction> max_error, const ErrorFunction& error_func = ErrorFunction(),
			const std::invoke_result_t<UnaryFunction, RealType>& zero = std::invoke_result_t<UnaryFunction, RealType>{}) 
		{
			using result_type = decltype(function(begin));
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

				if constexpr (print_steps) {
					std::cout << "I_n = " << old_value << "\nI_(n+1) = " << new_value << "\n\t||\terror = " << current_error 
						<< "\tCurrent step = " << step << std::endl;
				}

				old_value = new_value;
			} while (current_error > max_error);
			return new_value;
		}
	};
}