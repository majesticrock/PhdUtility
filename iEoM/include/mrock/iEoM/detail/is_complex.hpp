#pragma once
#include <type_traits>
#include <complex>

namespace mrock::iEoM::detail {

template <class T>
struct is_complex : std::false_type {};

template <class T>
struct is_complex<std::complex<T>> : std::true_type {};

template <class T>
inline constexpr bool is_complex_v = is_complex<std::remove_cvref_t<T>>::value;

template <class T>
concept Complex = is_complex_v<T>;
}