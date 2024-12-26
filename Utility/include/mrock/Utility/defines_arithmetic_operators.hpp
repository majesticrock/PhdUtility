#pragma once
#include <utility>
#include <type_traits>

namespace mrock::Utility {
	namespace _internal_utility {
		template<class _internal>
		using plus_type = std::decay_t<decltype((std::declval<_internal>().operator+=(std::declval<_internal>())))>;
		template<class _internal>
		using minus_type = std::decay_t<decltype((std::declval<_internal>().operator-=(std::declval<_internal>())))>;
		template<class _internal>
		using product_type = std::decay_t<decltype((std::declval<_internal>().operator*=(std::declval<_internal>())))>;
		template<class _internal>
		using division_type = std::decay_t<decltype((std::declval<_internal>().operator/=(std::declval<_internal>())))>;

		template < template< class T > class operation_type, class T >
		class _defines_operation {
			template <class _internal>
			static constexpr std::true_type test(operation_type<_internal>*) {
				return std::true_type{};
			};
			template <class _internal>
			static constexpr std::false_type test(...) {
				return std::false_type{};
			};
		public:
			static constexpr bool value = std::is_arithmetic<T>::value || decltype(test<T>(nullptr))::value;
		};
	}

	// value is true if T defines operator+=(T) and false otherwise
	template <class T>
	using defines_plus = _internal_utility::_defines_operation<_internal_utility::plus_type, T>;
	// retruns true if T defines operator+=(T) and false otherwise
	template <class T>
	constexpr bool defines_plus_v() { return defines_plus<T>::value; };

	// value is true if T defines operator-=(T) and false otherwise
	template <class T>
	using defines_minus = _internal_utility::_defines_operation<_internal_utility::minus_type, T>;
	// retruns true if T defines operator-=(T) and false otherwise
	template <class T>
	constexpr bool defines_minus_v() { return defines_minus<T>::value; };

	// value is true if T defines operator*=(T) and false otherwise
	template <class T>
	using defines_product = _internal_utility::_defines_operation<_internal_utility::product_type, T>;
	// retruns true if T defines operator*=(T) and false otherwise
	template <class T>
	constexpr bool defines_product_v() { return defines_product<T>::value; };

	// value is true if T defines operator/=(T) and false otherwise
	template <class T>
	using defines_division = _internal_utility::_defines_operation<_internal_utility::division_type, T>;
	// retruns true if T defines operator/=(T) and false otherwise
	template <class T>
	constexpr bool defines_division_v() { return defines_division<T>::value; };

	// value is true if T defines both += and -= operations and false otherwise
	template <class T>
	struct is_linearly_combinable {
		static constexpr bool value = defines_minus<T>::value && defines_plus<T>::value;
	};
	// retruns true if T defines both += and -= operations and false otherwise
	template <class T>
	constexpr bool is_linearly_combinable_v() { return is_linearly_combinable<T>::value; };
}

//Debug and test area
/*
#include <complex>
#include <iostream>
struct has_plus {
	int i;

	has_plus& operator+=(const has_plus& rhs) {
		return *this;
	};
	has_plus() = delete;
	has_plus(int _i) : i(_i) {};
};

struct has_plus_minus : public has_plus {
	has_plus_minus& operator-=(const has_plus_minus& rhs) {
		return *this;
	}
	has_plus_minus(int _i) : has_plus(_i) {};
};

struct no_plus {};
enum enum_test { one, two, three };

void test() {
	using namespace mrock::Utility;

	std::cout << "std::complex<double>: " << is_linearly_combinable<std::complex<double>>::value << std::endl;
	std::cout << "has_plus: " << is_linearly_combinable<has_plus>::value << std::endl;
	std::cout << "has_plus_minus: " << is_linearly_combinable<has_plus_minus>::value << std::endl;
	std::cout << "no_plus: " << is_linearly_combinable<no_plus>::value << std::endl;
	std::cout << "double: " << is_linearly_combinable<double>::value << std::endl;
	std::cout << "enum: " << is_linearly_combinable<enum_test>::value << std::endl;
}
*/