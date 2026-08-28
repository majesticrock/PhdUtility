#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_TERMCOLLECTOR_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_TERMCOLLECTOR_HPP

#include "AbstractCollector.hpp"
#include "Term.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace mrock::symbolic_operators {

/**
 * @class TermCollector
 * @brief A wrapper for a vector of Term objects.
 */
struct TermCollector : public AbstractCollector<Term> {
    /**
     * @brief Serializes the TermCollector object.
     *
     * @tparam Archive The type of the archive.
     * @param ar The archive object.
     * @param version The version of the serialization.
     */
    template <class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & terms;
    };

    MROCK_FORWARD_CONSTRUCTORS(TermCollector, AbstractCollector<Term>);

    MROCK_FORWARD_ASSIGNMENT(TermCollector, terms);

    /**
     * @brief Normal orders the terms by using the canoncical (anti-)commutation relations
     * The result is stored in the input vector.
     * A simple example is
     * \f$ b b^\dagger = 1 \pm b^\dagger b \f$, where the + applies to bosons and the minus to fermions.
     * @param terms The terms to normal order.
     */
    void normal_order();

    /**
     * @brief Sorts the terms, adds identical ones together and removes those that are equal to 0.
     */
    void clean_up();

    /**
     * @brief Applies the Hermitian conjugate to \c *this
     */
    void hermitian_conjugate();

    /**
     * @brief Compute other * this
     * IMPOARTANT: The result will not be normal ordered! If you require
     * a normal ordered expression, please call normal_order!
     *
     * @param other the other TermCollector object.
     * @return Reference to \c *this containing the result.
     */
    TermCollector& multiply_from_the_left(const TermCollector& other);

    /**
     * @brief Compute this * other
     * IMPOARTANT: The result will not be normal ordered! If you require
     * a normal ordered expression, please call normal_order!
     *
     * @param other the other TermCollector object.
     * @return Reference to \c *this containing the result.
     */
    TermCollector& multiply_from_the_right(const TermCollector& other);

    /**
     * @brief Converts \c *this to a string without the prefactor.
     *
     * @return The string representation.
     */
    std::string to_string_without_prefactor() const;

    /**
     * @brief Adds the right-hand side terms to the left-hand side terms.
     *
     * @param other The right-hand side terms to add.
     * @return A reference to the modified left-hand side terms in \c *this.
     */
    TermCollector& operator+=(const TermCollector& other);

    /**
     * @brief Subtracts the right-hand side terms from the left-hand side terms.
     *
     * @param other The right-hand side terms to add.
     * @return A reference to the modified left-hand side terms in \c *this.
     */
    TermCollector& operator-=(const TermCollector& other);

    /**
     * @brief Multiplies the left-hand side term vector by the right-hand side term vector.
     *
     * @param other The right-hand side terms to add.
     * @return A reference to the modified left-hand side terms in \c *this.
     */
    TermCollector& operator*=(const TermCollector& other);
};

/**
 * @brief Computes the commutator of two terms: \f$ [A, B] = AB - BA \f$.
 * @param left The left term.
 * @param right The right term.
 * @return The commutation result.
 */
TermCollector commutator(const Term& left, const Term& right);

/**
 * @brief Computes the commutator of two sets of terms: \f$ [A, B] = AB - BA \f$.
 * @param left The left-hand side terms.
 * @param right The right-hand side terms.
 * @return The result of [left, right]
 */
TermCollector commutator(const TermCollector& left, const TermCollector& right);

/**
 * @brief Computes the commutator of a term and a set of terms: \f$ [A, B] = AB - BA \f$.
 * @param left The left-hand side term.
 * @param right The right-hand side terms.
 * @return The result of [left, right]
 */
TermCollector commutator(const Term& left, const TermCollector& right);

/**
 * @brief Computes the commutator of a set of terms and a term: \f$ [A, B] = AB - BA \f$.
 * @param left The left-hand side terms.
 * @param right The right-hand side term.
 * @return The result of [left, right]
 */
TermCollector commutator(const TermCollector& left, const Term& right);

/**
 * @brief Negates each term in the vector.
 *
 * @param terms The terms to negate.
 * @return A vector containing the negated terms.
 */
TermCollector operator-(TermCollector terms);

/**
 * @brief Adds two vectors of terms and returns the result.
 *
 * @param lhs The left-hand side terms.
 * @param rhs The right-hand side terms.
 * @return A vector containing the sum of the operands.
 */
inline TermCollector operator+(TermCollector lhs, const TermCollector& rhs) {
    return (lhs += rhs);
}

/**
 * @brief Subtracts the right-hand side term vector from the left-hand side term vector.
 *
 * @param lhs The left-hand side terms.
 * @param rhs The right-hand side terms.
 * @return A vector containing the difference of the operands.
 */
inline TermCollector operator-(TermCollector lhs, const TermCollector& rhs) {
    return (lhs -= rhs);
}

/**
 * @brief Multiplies two vectors of terms and returns the result.
 *
 * @param lhs The left-hand side terms.
 * @param rhs The right-hand side terms.
 * @return A vector containing the product of the operands.
 */
inline TermCollector operator*(TermCollector lhs, const TermCollector& rhs) {
    return (lhs *= rhs);
}

/**
 * @brief Outputs a vector of terms to a stream.
 * @param os The output stream.
 * @param terms The terms.
 * @return The output stream.
 */
std::ostream& operator<<(std::ostream& os, const TermCollector& terms);

}  // namespace mrock::symbolic_operators

#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_TERMCOLLECTOR_HPP