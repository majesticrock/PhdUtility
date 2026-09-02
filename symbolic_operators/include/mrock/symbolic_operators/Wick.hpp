#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICK_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICK_HPP
/**
 * @file Wick.hpp
 * @brief Functions for applying Wick's theorem and manipulating Wick terms.
 */

#include "TermCollector.hpp"
#include "WickOperatorTemplate.hpp"
#include "WickSymmetry.hpp"
#include "WickTermCollector.hpp"
#include <vector>

namespace mrock::symbolic_operators {

/**
 * @brief Identifies Wick operators in a given Wick term.
 *
 * @param source The source Wick term.
 * @param operator_templates The vector of Wick operator templates.
 * @return WickTermCollector The collected Wick terms.
 */
WickTermCollector identify_wick_operators(const WickTerm& source,
                                          const std::vector<WickOperatorTemplate>& operator_templates);

/**
 * @brief Applies Wick's theorem to a set of terms.
 *
 * @param terms The vector of terms.
 * @param operator_templates The vector of Wick operator templates.
 * @param reciever The WickTermCollector to receive the results.
 */
void wicks_theorem(const TermCollector& terms,
                   const std::vector<WickOperatorTemplate>& operator_templates,
                   WickTermCollector& reciever);
}  // namespace mrock::symbolic_operators
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICK_HPP
