/**
 * @file ErrorFunctors.hpp
 * @brief Provides some functors for computing errors
 */

#pragma once
#include <numeric>
#include <limits>
#include <functional>
#include <cmath>

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
            if constexpr (compute_relative_error) {
			    const value_type abs_new = norm(_new);
			    return (abs_new > sqrt(std::numeric_limits<value_type>::epsilon()) ?  norm(_new - _old) / abs_new : norm(_new - _old) / _new.size());
            }
			else {
                return norm(_new - _old) / _new.size();
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
		value_type operator()(const vector_type& _new, const vector_type& _old) const {
            const value_type abs_diff = std::transform_reduce(_new.begin(), _new.end(), _old.begin(), value_type{}, 
										std::plus<value_type>(), [](const value_type& left, const value_type& right){ using std::abs; return abs(left - right); });
			
            if constexpr (compute_relative_error) {
                const value_type abs_new = std::transform_reduce(_new.begin(), _new.end(), value_type{}, 
										std::plus<value_type>(), [](const value_type& val){ using std::abs; return abs(val); });
			    return (abs_new > sqrt(std::numeric_limits<value_type>::epsilon()) ?  abs_diff / abs_new : abs_diff / _new.size());
            }
            else {
                return abs_diff / _new.size();
            }
		}
	};
}