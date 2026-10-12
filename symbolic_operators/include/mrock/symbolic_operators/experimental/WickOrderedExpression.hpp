#pragma once

#include <mrock/symbolic_operators/detail/vector_macro.hpp>
#include <mrock/symbolic_operators/Operator.hpp>

#include <vector>
#include <iosfwd>

namespace mrock::symbolic_operators::experimental {

// Represents a string of Wick-ordered (or normal-ordered) operators, e.g., :c^dagger c^dagger c c:
struct WickOrderedExpression {
    std::vector<Operator> operators;  ///< The operators in the Wick ordered expression.

    WickOrderedExpression& hermitian_conjugate_inplace() noexcept;
    WickOrderedExpression hermitian_conjugate() const noexcept;

    MROCK_VECTOR_WRAPPER_FILL_MEMBERS(Operator, operators);
};

/**
 * @brief Outputs a Wick ordered expression to a stream.
 * @param os The output stream.
 * @param expression The Wick ordered expression.
 * @return The output stream.
 */
std::ostream& operator<<(std::ostream& os, const WickOrderedExpression& expression);

/**
 * @brief Compares two Wick ordered expressions for equality.
 * @param lhs The left-hand side expression.
 * @param rhs The right-hand side expression.
 * @return true if the expressions are equal, false otherwise.
 */
bool operator==(const WickOrderedExpression& lhs, const WickOrderedExpression& rhs);

/**
 * @brief Compares two Wick ordered expressions for inequality.
 * @param lhs The left-hand side expression.
 * @param rhs The right-hand side expression.
 * @return true if the expressions are not equal, false otherwise.
 */
bool operator!=(const WickOrderedExpression& lhs, const WickOrderedExpression& rhs);

} // namespace mrock::symbolic_operators::experimental