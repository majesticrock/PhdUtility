#pragma once
#include <type_traits>
#include <complex>

namespace mrock::utility {
	template <class T>
	struct IsComplex_t : public std::false_type {};

	template <class T>
	struct IsComplex_t<std::complex<T>> : public std::true_type {};

	template <class T>
	constexpr bool is_complex() {
		return IsComplex_t<T>::value;
	};

	template <class T>
	concept Complex = IsComplex_t<T>::value;

	template <class T>
	concept Real = !IsComplex_t<T>::value;
};