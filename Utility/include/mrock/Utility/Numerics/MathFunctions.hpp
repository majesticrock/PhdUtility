#pragma once
#include <type_traits>

namespace mrock::Utility::Numerics {
	// Computes n!!!... (m times)
	template <unsigned int factorial_m, typename IntegerType, typename RealType = size_t>
	constexpr RealType m_fold_factorial(const IntegerType n) {
		RealType result{ 1 };
		for (std::make_signed_t<IntegerType> i = n; i > 0; i -= factorial_m) {
			result *= i;
		}
		return result;
	}

	// Computes n!
	template <typename IntegerType, typename RealType = size_t>
	constexpr RealType factorial(const IntegerType n) {
		return m_fold_factorial<1>(n);
	}

	// Computes n!!
	template <typename IntegerType, typename RealType = size_t>
	constexpr RealType double_factorial(const IntegerType n) {
		return m_fold_factorial<2>(n);
	}
}