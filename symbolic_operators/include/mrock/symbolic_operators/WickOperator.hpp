#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICKOPERATOR_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICKOPERATOR_HPP
/**
 * @file WickOperator.hpp
 * @brief Defines the WickOperator structure used in symbolic operators.
 */

#include "IndexWrapper.hpp"
#include "Momentum.hpp"
#include "MomentumSymbol.hpp"
#include "Operator.hpp"
#include "OperatorType.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace mrock::symbolic_operators {

/**
 * @class WickOperator
 * @brief A structure representing a Wick operator.
 */
struct WickOperator {
    OperatorType type{OperatorType::Undefined};  ///< The type of the operator.
    bool is_daggered{};                          ///< Indicates if the operator is daggered.
    Momentum momentum;                           ///< The momentum associated with the operator.
    IndexWrapper indices;                        ///< The indices associated with the operator.

    /**
     * @brief Serializes the WickOperator object.
     *
     * @tparam Archive The type of the archive.
     * @param ar The archive object.
     * @param version The version of the serialization.
     */
    template <class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & type;
        ar & is_daggered;
        ar & momentum;
        ar & indices;
    }

    /**
     * @brief Constructs a WickOperator object.
     *
     * @param _type The type of the operator.
     * @param _is_daggered Whether the operator is daggered.
     * @param _momentum The momentum of the operator.
     * @param _indices The indices of the operator.
     */
    WickOperator(const OperatorType& _type,
                 const bool _is_daggered,
                 const Momentum& _momentum,
                 const IndexWrapper& _indices = IndexWrapper());

    /**
     * @brief Constructs a WickOperator object.
     *
     * @param _type The type of the operator.
     * @param _is_daggered Whether the operator is daggered.
     * @param _momentum The momentum of the operator.
     * @param _index The index of the operator.
     */
    WickOperator(const OperatorType& _type, const bool _is_daggered, const Momentum& _momentum, const Index _index);

    /**
     * @brief Default constructor for WickOperator.
     */
    WickOperator() = default;

    /**
     * @brief Constructs a WickOperator object from a string expression.
     *
     * @param expression The string expression.
     */
    WickOperator(const std::string& expression);

    /**
     * @brief Checks if the operator uses a specific index.
     *
     * @param index The index to check.
     * @return true if the operator uses the index.
     * @return false otherwise.
     */
    inline bool uses_index(const Index index) const noexcept;

    /**
     * @brief Checks if the operator depends on a specific momentum.
     *
     * @param momentum The momentum to check.
     * @return true if the operator depends on the momentum.
     * @return false otherwise.
     */
    inline bool depends_on(const MomentumSymbol::name_type momentum) const noexcept;

    /**
     * @brief Removes a momentum contribution from the operator.
     *
     * @param value The momentum value to remove.
     */
    inline void remove_momentum_contribution(const MomentumSymbol::name_type value);

    /**
     * @brief Transforms \c this to the equivalent operator expression.
     * Example: < \c OperatorType::SC , \c Momentum k, \c is_dagger \c false > becomes < c_{-k down} c_{k up} >
     *
     * @return The transformed expression
     */
    std::vector<Operator> to_operator_expression() const;
};

/**
 * @brief Equality operator for WickOperator.
 *
 * @param lhs The left-hand side WickOperator.
 * @param rhs The right-hand side WickOperator.
 * @return true if the two WickOperator objects are equal.
 * @return false otherwise.
 */
bool operator==(const WickOperator& lhs, const WickOperator& rhs);

/**
 * @brief Inequality operator for WickOperator.
 *
 * @param lhs The left-hand side WickOperator.
 * @param rhs The right-hand side WickOperator.
 * @return true if the two WickOperator objects are not equal.
 * @return false otherwise.
 */
bool operator!=(const WickOperator& lhs, const WickOperator& rhs);

/**
 * @brief Stream insertion operator for WickOperator.
 *
 * @param os The output stream.
 * @param op The WickOperator object.
 * @return std::ostream& The updated output stream.
 */
std::ostream& operator<<(std::ostream& os, const WickOperator& op);

/**
 * @brief Stream insertion operator for a vector of WickOperator objects.
 *
 * @param os The output stream.
 * @param ops The vector of WickOperator objects.
 * @return std::ostream& The updated output stream.
 */
std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops);

/**
 * @brief Compares two WickOperator objects for greater-than ordering.
 * @param lhs The left-hand side WickOperator.
 * @param rhs The right-hand side WickOperator.
 * @return True if lhs is greater than rhs, false otherwise.
 */
bool operator>(const WickOperator& lhs, const WickOperator& rhs);

/**
 * @brief Compares two WickOperator objects for less-than ordering.
 * @param lhs The left-hand side WickOperator.
 * @param rhs The right-hand side WickOperator.
 * @return True if lhs is less than rhs, false otherwise.
 */
bool operator<(const WickOperator& lhs, const WickOperator& rhs);

/**
 * @brief Compares two WickOperator objects for greater-or-equal ordering.
 * @param lhs The left-hand side WickOperator.
 * @param rhs The right-hand side WickOperator.
 * @return True if lhs is greater than or equal to rhs, false otherwise.
 */
bool operator>=(const WickOperator& lhs, const WickOperator& rhs);

/**
 * @brief Compares two WickOperator objects for less-or-equal ordering.
 * @param lhs The left-hand side WickOperator.
 * @param rhs The right-hand side WickOperator.
 * @return True if lhs is less than or equal to rhs, false otherwise.
 */
bool operator<=(const WickOperator& lhs, const WickOperator& rhs);

// Inline definitions
bool WickOperator::uses_index(const Index index) const noexcept {
    for (const auto& idx : this->indices) {
        if (idx == index)
            return true;
    }
    return false;
}
inline bool WickOperator::depends_on(const MomentumSymbol::name_type momentum) const noexcept {
    return this->momentum.is_used_at(momentum) != -1;
}
inline void WickOperator::remove_momentum_contribution(const MomentumSymbol::name_type value) {
    momentum.remove_contribution(value);
}
}  // namespace mrock::symbolic_operators
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICKOPERATOR_HPP
