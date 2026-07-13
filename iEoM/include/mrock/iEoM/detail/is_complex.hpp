#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_IS_COMPLEX_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_IS_COMPLEX_HPP
#include <type_traits>
#include <complex>

namespace mrock::iEoM::detail {

/**
 * @brief Type trait that determines whether a type is std::complex.
 *
 * Defaults to false for all types, and is specialized to true for std::complex.
 *
 * @tparam T Type to inspect.
 */
template <class T>
struct is_complex : std::false_type {};

/**
 * @brief Specialization of is_complex for std::complex types.
 *
 * @tparam T Underlying real type inside the complex wrapper.
 */
template <class T>
struct is_complex<std::complex<T>> : std::true_type {};

/**
 * @brief Boolean helper indicating whether a type is complex.
 *
 * Removes cv/ref qualifiers before evaluation.
 *
 * @tparam T Type to inspect.
 */
template <class T>
inline constexpr bool is_complex_v = is_complex<std::remove_cvref_t<T>>::value;

/**
 * @brief Concept satisfied by complex scalar types.
 *
 * @tparam T Type to inspect.
 */
template <class T>
concept Complex = is_complex_v<T>;
}
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_IS_COMPLEX_HPP
