#pragma once
#include "WickOrderedExpression.hpp"

#include <mrock/symbolic_operators/AbstractTerm.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>

#include <compare>
#include <ostream>
#include <vector>

namespace mrock::symbolic_operators::experimental {

// Represents a Wick-ordered (or normal-ordered) Term, including the appropriate expectation values if
// the ordering procedure is not done with respect to the vacuum state, e.g.,
// <n_k> :n_l:
class WickOrderedTerm : public AbstractTerm<WickOperator> {
public:
    void for_each_momentum_except_deltas(const std::function<void(Momentum&)>& f) override;

    void for_each_index_except_deltas(const std::function<void(IndexWrapper&)>& f) override;

    WickOrderedExpression wick_expression;  ///< The Wick ordered expression associated with the term.
    
    /**
     * @brief Default constructor for WickOrderedTerm.
     */
    WickOrderedTerm() = default;

    /**
     * @brief Constructs a WickOrderedTerm object using the coefficients, sums, etc provided from a base Term
     * 
     * @param base The base Term object.
     */
    WickOrderedTerm(const Term& base);

    /**
     * @brief Constructs a WickOrderedTerm object using the coefficients, sums, etc provided from a base WickTerm
     * 
     * @param base The base WickTerm object.
     */
    WickOrderedTerm(const WickTerm& base);

    /**
     * @brief Constructs a Term with a summation over momenta and spins and a coefficient
     *
     * @param _multiplicity The _multiplicity of the term
     * @param _coefficient The coefficient
     * @param _sums The sums
     * @param _wick_operator_string The operators of the term, i.e., :cccc... :
     * @param _expectation_values The expectation values, i.e., <cc> ...
     */
    WickOrderedTerm(const IntFractional& _multiplicity,
        const Coefficient& _coefficient,
        const SumContainer& _sums,
        const std::vector<Operator>& _wick_operator_string,
        const std::vector<WickOperator>& _expectation_values);

    /**
     * @brief Constructs a WickOrderedTerm object using the coefficients, sums, etc provided from a base WickOrderedTerm, 
     * while ignoring the Wick expression of the base term.
     * 
     * @param base The base WickTerm object.
     */
    static WickOrderedTerm from_wick_term_skip_wick_expression(const WickOrderedTerm& base);

    /**
     * @brief Transforms the expectation values into Coefficient objects and fills the operators of the term
     * with the WickOrderedExpression.
     * 
     * @return The resulting Term object
     */
    Term to_term(const std::string& base) const;

    /**
     * @brief Permutes the operator string in wick_expression such that all creation operators are on the left,
     * and all annihilation operators on the right. The function keeps track of sign changes due to fermionic exchange.
     * No new terms are generated because the operators anticommute under Wick-ordering.
     * 
     * @return Returns true if everything ran successfully.
     * Returns false if two identical fermionic creation or annihilation operators exist in the term.
     */
    bool creation_operators_to_the_left();

    /**
     * @brief Includes a template result in the term.
     *
     * @param result The TemplateResult::SingleResult to include.
     */
    void include_template_result(const TemplateResult::SingleResult& result);

    /**
     * @brief Resolves the deltas in the term.
     */
    bool resolve_deltas();

    /**
     * @brief Discards momenta that are zero from the term.
     */
    void discard_zero_momenta();

    /**
     * @brief Sorts the operators in the term.
     */
    void sort();

    /**
     * @brief Returns true if the term is bilinear
     * 
     * @return true if the term is bilinear
     */
    bool is_identity() const noexcept;

    /**
     * @brief Returns true if the term is bilinear
     * 
     * @return true if the term is bilinear
     */
    bool is_bilinear() const noexcept;

    /**
     * @brief Returns true if the term is quartic
     * 
     * @return true if the term is quartic
     */
    bool is_quartic() const noexcept;

    /**
     * @brief Takes the Hermitian conjugate of the current expression in place
     * 
     * @return A reference to \c *this after the Hermitian conjugation 
     */
    WickOrderedTerm& hermitian_conjugate_inplace() noexcept;

    /**
     * @brief Takes the Hermitian conjugate of the current expression in a copy
     * 
     * @return A copy of \c *this after the Hermitian conjugation 
     */
    WickOrderedTerm hermitian_conjugate() const noexcept;

    std::strong_ordering operator<=>(const WickOrderedTerm& other) const;

    bool operator==(const WickOrderedTerm& other) const;
};

/**
 * @brief Overloads the output stream operator for WickOrderedTerm.
 * 
 * @param os The output stream.
 * @param term The WickOrderedTerm to output.
 * @return The output stream.
 */
std::ostream& operator<<(std::ostream& os, const WickOrderedTerm& term);

} // namespace mrock::symbolic_operators::experimental