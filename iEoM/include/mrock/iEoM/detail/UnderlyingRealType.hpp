#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_UNDERLYINGREALTYPE_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_UNDERLYINGREALTYPE_HPP
#include <complex>
#include <type_traits>

namespace mrock::iEoM::detail {

/**
 * @brief Determine the underlying real scalar type for a numeric type.
 *
 * This trait returns the same type for real scalar types and strips the
 * complex wrapper for std::complex types.
 *
 * @tparam T Input scalar type.
 */
template <class T>
struct UnderlyingRealType : std::type_identity<T> {};

/**
 * @brief Specialization for std::complex types.
 *
 * Extracts the underlying real scalar type from a complex type.
 *
 * @tparam T Underlying real component type.
 */
template <class T>
struct UnderlyingRealType<std::complex<T>> : std::type_identity<T> {};

/**
 * @brief Helper alias for obtaining the underlying real type.
 *
 * @tparam T Input scalar type, with cv/ref qualifiers removed.
 */
template <class T>
using UnderlyingRealType_t = typename UnderlyingRealType<std::remove_cvref_t<T>>::type;
}  // namespace mrock::iEoM::detail
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_UNDERLYINGREALTYPE_HPP
