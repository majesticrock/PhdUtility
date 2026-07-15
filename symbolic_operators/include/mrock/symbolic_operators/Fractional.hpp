#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_FRACTIONAL_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_FRACTIONAL_HPP
/**
 * @file Fractional.hpp
 * @brief A simple compile-time rational number wrapper for integer types.
 *
 * This header provides a templated `Fractional` type for representing and
 * manipulating rational numbers using integer numerators and denominators.
 * It supports arithmetic operations, comparisons, conversion to built-in
 * numeric types, and ostream formatting.
 */

#include <iostream>
#include <numeric>
#include <type_traits>

namespace mrock::symbolic_operators {

/**
 * @brief Represents a rational number as a numerator and denominator.
 *
 * @tparam _int The integer type used for numerator and denominator.
 */
template <class _int>
struct Fractional {
    _int numerator{};     ///< Numerator of the fraction.
    _int denominator{1};  ///< Denominator of the fraction. Defaults to 1.

    /**
     * @brief Default constructs the fraction to 0/1.
     */
    constexpr Fractional(){};

    /**
     * @brief Constructs a fraction from an integer value.
     *
     * @param integer The integer numerator; denominator is set to 1.
     */
    constexpr Fractional(_int integer) : numerator(integer){};

    /**
     * @brief Constructs a fraction from a numerator and denominator.
     *
     * @param _numerator The numerator.
     * @param _denominator The denominator.
     */
    constexpr Fractional(_int _numerator, _int _denominator) : numerator(_numerator), denominator(_denominator){};

    /**
     * @brief Serializes the fraction using Boost.Serialization-like archives.
     *
     * @tparam Archive The archive type.
     * @param ar The archive to serialize to or deserialize from.
     * @param version Archive version number.
     */
    template <class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & numerator;
        ar & denominator;
    }

    /**
     * @brief Returns the fraction as a floating-point value.
     *
     * @tparam _float The floating-point type to convert to.
     * @return The value of numerator / denominator.
     */
    template <class _float = double>
    constexpr _float value() const {
        return static_cast<_float>(numerator) / static_cast<_float>(denominator);
    }

    /**
     * @brief Checks whether the fraction represents an integer.
     *
     * @return True if numerator is an exact multiple of denominator.
     */
    constexpr bool is_integer() const { return (numerator % denominator == 0); }

    /**
     * @brief Reduces the fraction to lowest terms.
     *
     * If the underlying integer type is signed, the denominator is kept
     * positive by moving the sign into the numerator.
     */
    inline void reduce_fraction() {
        if constexpr (std::is_signed_v<_int>) {
            if (denominator < 0) {
                numerator *= -1;
                denominator *= -1;
            }
        }
        const _int gcd = std::gcd(numerator, denominator);
        numerator /= gcd;
        denominator /= gcd;
    }

    /**
     * @brief Converts the fraction to an integral type by truncating.
     *
     * @tparam Integer The target integral type.
     * @return Truncated integer division result.
     */
    template <typename Integer, std::enable_if_t<std::is_integral<Integer>::value, bool> = true>
    constexpr operator Integer() const {
        return numerator / denominator;
    }

    /**
     * @brief Converts the fraction to a floating-point type.
     *
     * @tparam Floating The target floating-point type.
     * @return Floating-point division result.
     */
    template <typename Floating, std::enable_if_t<std::is_floating_point<Floating>::value, bool> = true>
    constexpr operator Floating() const {
        return value<Floating>();
    }

    /**
     * @brief Three-way comparison with another Fractional.
     *
     * @param other The fraction to compare with.
     * @return Comparison result based on cross-multiplication.
     */
    constexpr auto operator<=>(const Fractional<_int>& other) const {
        return (this->numerator * other.denominator) <=> (other.numerator * this->denominator);
    }

    /**
     * @brief Checks equality with another Fractional.
     *
     * @param other The fraction to compare with.
     * @return True if both fractions represent the same rational number.
     */
    constexpr bool operator==(const Fractional<_int>& other) const {
        return (this->numerator * other.denominator) == (other.numerator * this->denominator);
    };

    /**
     * @brief Checks inequality with another Fractional.
     *
     * @param other The fraction to compare with.
     * @return True if the fractions represent different values.
     */
    constexpr bool operator!=(const Fractional<_int>& other) const { return !(*this == other); };

    /**
     * @brief Three-way comparison with an integer.
     *
     * @param other The integer to compare with.
     * @return Comparison result using the fraction's denominator.
     */
    constexpr auto operator<=>(_int other) const { return this->numerator <=> (other * this->denominator); }

    /**
     * @brief Checks equality with an integer.
     *
     * @param other The integer to compare with.
     * @return True if the fraction equals the integer.
     */
    constexpr bool operator==(_int other) const { return this->numerator == (other * this->denominator); };

    /**
     * @brief Checks inequality with an integer.
     *
     * @param other The integer to compare with.
     * @return True if the fraction does not equal the integer.
     */
    constexpr bool operator!=(_int other) const { return !(*this == other); };

    /**
     * @brief Adds another Fractional to this one.
     *
     * @param other The fraction to add.
     * @return Reference to this fraction.
     */
    inline Fractional& operator+=(Fractional const& other) {
        if (other.denominator == this->denominator) {
            this->numerator += other.numerator;
            this->reduce_fraction();
            return *this;
        }
        this->numerator = this->numerator * other.denominator + other.numerator * this->denominator;
        this->denominator *= other.denominator;
        this->reduce_fraction();
        return *this;
    }

    /**
     * @brief Subtracts another Fractional from this one.
     *
     * @param other The fraction to subtract.
     * @return Reference to this fraction.
     */
    inline Fractional& operator-=(Fractional const& other) {
        if (other.denominator == this->denominator) {
            this->numerator -= other.numerator;
            this->reduce_fraction();
            return *this;
        }
        this->numerator = this->numerator * other.denominator - other.numerator * this->denominator;
        this->denominator *= other.denominator;
        this->reduce_fraction();
        return *this;
    }

    /**
     * @brief Multiplies this fraction by another Fractional.
     *
     * @param other The fraction multiplier.
     * @return Reference to this fraction.
     */
    inline Fractional& operator*=(Fractional const& other) {
        this->denominator *= other.denominator;
        this->numerator *= other.numerator;
        this->reduce_fraction();
        return *this;
    }

    /**
     * @brief Divides this fraction by another Fractional.
     *
     * @param other The fraction divisor.
     * @return Reference to this fraction.
     */
    inline Fractional& operator/=(Fractional const& other) {
        this->denominator *= other.numerator;
        this->numerator *= other.denominator;
        this->reduce_fraction();
        return *this;
    }

    /**
     * @brief Adds an integer to the fraction.
     *
     * @param other The integer to add.
     * @return Reference to this fraction.
     */
    constexpr Fractional& operator+=(_int other) {
        this->numerator += other * this->denominator;
        return *this;
    }

    /**
     * @brief Subtracts an integer from the fraction.
     *
     * @param other The integer to subtract.
     * @return Reference to this fraction.
     */
    constexpr Fractional& operator-=(_int other) {
        this->numerator -= other * this->denominator;
        return *this;
    }

    /**
     * @brief Multiplies the fraction by an integer.
     *
     * @param other The integer multiplier.
     * @return Reference to this fraction.
     */
    inline Fractional& operator*=(_int other) {
        this->numerator *= other;
        this->reduce_fraction();
        return *this;
    }

    /**
     * @brief Divides the fraction by an integer.
     *
     * @param other The integer divisor.
     * @return Reference to this fraction.
     */
    inline Fractional& operator/=(_int other) {
        this->denominator *= other;
        this->reduce_fraction();
        return *this;
    }
};

/**
 * @brief Computes an integer power of a Fractional value.
 *
 * @tparam _int The integer type used by the fraction.
 * @param base The base fraction.
 * @param exponent The exponent to raise to.
 * @return The result of raising the base to the exponent.
 */
template <class _int>
constexpr Fractional<_int> pow(Fractional<_int> base, int exponent) {
    if (exponent == 0)
        return Fractional<_int>{_int(1), _int(1)};
    if (exponent < 0)
        return pow(Fractional<_int>{base.denominator, base.numerator}, exponent);
    if (exponent == 1)
        return base;
    if (exponent & 1) {
        return base * pow(base) * pow(base);
    } else {
        return pow(base) * pow(base);
    }
}

/**
 * @brief Formats a fraction for output.
 *
 * @param os Output stream.
 * @param frac Fraction to format.
 * @return The output stream.
 */
template <class _int>
std::ostream& operator<<(std::ostream& os, const Fractional<_int>& frac) {
    if (frac.is_integer()) {
        os << frac.numerator / frac.denominator;
    } else {
        if (frac.numerator < 0) {
            os << "-";
        }
        os << "(" << std::abs(frac.numerator) << "/" << frac.denominator << ")";
    }
    return os;
}

/**
 * @name Fractional arithmetic operators
 * @{
 */
template <class _int>
inline Fractional<_int> operator+(Fractional<_int> lhs, const Fractional<_int>& rhs) {
    return lhs += rhs;
}
template <class _int>
inline Fractional<_int> operator-(Fractional<_int> lhs, const Fractional<_int>& rhs) {
    return lhs -= rhs;
}
template <class _int>
inline Fractional<_int> operator*(Fractional<_int> lhs, const Fractional<_int>& rhs) {
    return lhs *= rhs;
}
template <class _int>
inline Fractional<_int> operator/(Fractional<_int> lhs, const Fractional<_int>& rhs) {
    return lhs /= rhs;
}

template <class _int>
inline Fractional<_int> operator-(Fractional<_int> rhs) {
    return rhs *= _int{-1};
}

template <class _int>
constexpr Fractional<_int> operator+(Fractional<_int> lhs, _int rhs) {
    return lhs += rhs;
}
template <class _int>
constexpr Fractional<_int> operator-(Fractional<_int> lhs, _int rhs) {
    return lhs -= rhs;
}
template <class _int>
inline Fractional<_int> operator*(Fractional<_int> lhs, _int rhs) {
    return lhs *= rhs;
}
template <class _int>
inline Fractional<_int> operator/(Fractional<_int> lhs, _int rhs) {
    return lhs /= rhs;
}

template <class _int>
constexpr Fractional<_int> operator+(_int lhs, Fractional<_int> rhs) {
    return rhs += lhs;
}
template <class _int>
constexpr Fractional<_int> operator-(_int lhs, Fractional<_int> rhs) {
    return rhs -= lhs;
}
template <class _int>
inline Fractional<_int> operator*(_int lhs, Fractional<_int> rhs) {
    return rhs *= lhs;
}
template <class _int>
inline Fractional<_int> operator/(_int lhs, Fractional<_int> rhs) {
    return rhs /= lhs;
}

template <class _Number, class _int>
inline _Number operator+(const Fractional<_int>& lhs, const _Number& rhs) {
    return rhs + static_cast<_Number>(lhs);
}
template <class _Number, class _int>
inline _Number operator-(const Fractional<_int>& lhs, const _Number& rhs) {
    return rhs - static_cast<_Number>(lhs);
}
template <class _Number, class _int>
inline _Number operator*(const Fractional<_int>& lhs, const _Number& rhs) {
    return rhs * static_cast<_Number>(lhs);
}
template <class _Number, class _int>
inline _Number operator/(const Fractional<_int>& lhs, const _Number& rhs) {
    return rhs / static_cast<_Number>(lhs);
}
template <class _Number, class _int>
inline _Number operator+(const _Number& lhs, const Fractional<_int>& rhs) {
    return lhs + static_cast<_Number>(rhs);
}
template <class _Number, class _int>
inline _Number operator-(const _Number& lhs, const Fractional<_int>& rhs) {
    return lhs - static_cast<_Number>(rhs);
}
template <class _Number, class _int>
inline _Number operator*(const _Number& lhs, const Fractional<_int>& rhs) {
    return lhs * static_cast<_Number>(rhs);
}
template <class _Number, class _int>
inline _Number operator/(const _Number& lhs, const Fractional<_int>& rhs) {
    return lhs / static_cast<_Number>(rhs);
}
/** @} */

/**
 * @typedef IntFractional
 * @brief A Fractional type instantiated with `int`.
 */
typedef Fractional<int> IntFractional;
}  // namespace mrock::symbolic_operators
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_FRACTIONAL_HPP
