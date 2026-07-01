#pragma once
#include <complex>
#include <type_traits>

namespace mrock::utility {
template <class T>
struct UnderlyingRealType : std::type_identity<T> {};

template <class T>
struct UnderlyingRealType<std::complex<T>> : std::type_identity<T> {};

template <class T>
using UnderlyingRealType_t = typename UnderlyingRealType<std::remove_cvref_t<T>>::type;
}