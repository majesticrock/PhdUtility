/**
 * @file ErrorFunctors.hpp
 * @brief Provides some functors for computing errors
 */

#pragma once
#include <numeric>
#include <limits>
#include <functional>
#include <cmath>

#include "../UnderlyingFloatingPoint.hpp"

namespace mrock::utility::Numerics {
    /**
	 * @struct scalar_error
	 * @brief Computes the relative error between two values.
	 * 
	 * @tparam T Type of the values.
	 */
	template<class T, class RealType = double, bool compute_relative_error = true>
	struct scalar_error {
		/**
		 * @brief Computes the relative error between two values.
		 * 
		 * @param _new The new value.
		 * @param _old The old value.
		 * @return double The relative error.
		 */
		double operator()(const T& _new, const T& _old) const {
			using std::abs;
            if constexpr (compute_relative_error) {
                const double abs_new = abs(_new);
			    return (abs_new > sqrt(std::numeric_limits<RealType>::epsilon()) ? abs(_new - _old) / abs_new : abs(_new - _old));
            }
			else {
                return abs(_new - _old);
            }
		};
	};

	/**
	 * @struct vector_norm_error
	 * @brief Computes the error based on the norm of vectors. Not designed for complex numbers.
	 * 
	 * @tparam vector_type Type of the vector.
	 * @tparam value_type Type of the values in the vector.
	 */
	template<class vector_type, class value_type = typename vector_type::value_type, bool compute_relative_error = true>
	struct vector_norm_error {
	public:
		typedef UnderlyingFloatingPoint_t<value_type> error_type;
	private:
		/**
		 * @brief Computes the norm of a vector.
		 * 
		 * @param vec The vector.
		 * @return error_type The norm of the vector.
		 */
		static error_type norm(const vector_type& vec) {
			return sqrt(std::transform_reduce(vec.begin(), vec.end(), vec.begin(), error_type{}));
		}
	public:
		/**
		 * @brief Computes the error based on the norm of two vectors.
		 * 
		 * @param _new The new vector.
		 * @param _old The old vector.
		 * @return error_type The error based on the norm.
		 */
		error_type operator()(const vector_type& _new, const vector_type& _old) const {
            if constexpr (compute_relative_error) {
			    const error_type abs_new = norm(_new);
			    return (abs_new > sqrt(std::numeric_limits<error_type>::epsilon()) ?  norm(_new - _old) / abs_new : norm(_new - _old) / _new.size());
            }
			else {
                return norm(_new - _old) / static_cast<error_type>(_new.size());
            }
		}
	};

	/**
	 * @struct vector_elementwise_error
	 * @brief Computes the element-wise error between two vectors. Not designed for complex numbers.
	 * 
	 * @tparam vector_type Type of the vector.
	 * @tparam value_type Type of the values in the vector.
	 */
	template<class vector_type, class value_type = typename vector_type::value_type, bool compute_relative_error = true>
	struct vector_elementwise_error {
		/**
		 * @brief Computes the element-wise error between two vectors.
		 * 
		 * @param _new The new vector.
		 * @param _old The old vector.
		 * @return value_type The element-wise error.
		 */
		typedef UnderlyingFloatingPoint_t<value_type> error_type;

		error_type operator()(const vector_type& _new, const vector_type& _old) const {
            const error_type abs_diff = std::transform_reduce(_new.begin(), _new.end(), _old.begin(), error_type{}, 
										std::plus<error_type>(), [](const value_type& left, const value_type& right){ using std::abs; return abs(left - right); });
			
            if constexpr (compute_relative_error) {
                const error_type abs_new = std::transform_reduce(_new.begin(), _new.end(), error_type{}, 
										std::plus<error_type>(), [](const value_type& val){ using std::abs; return abs(val); });
			    return (abs_new > sqrt(std::numeric_limits<error_type>::epsilon()) ?  abs_diff / abs_new : abs_diff / _new.size());
            }
            else {
                return abs_diff / static_cast<error_type>(_new.size());
            }
		}
	};
}