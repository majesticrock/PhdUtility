#pragma once

#include "WickOrderedTerm.hpp"

#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/TermCollector.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>
#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/AbstractCollector.hpp>
#include <mrock/symbolic_operators/detail/vector_macro.hpp>

#include <memory>
#include <ostream>
#include <vector>


namespace mrock::symbolic_operators::experimental {

// Represents multiple Wick-ordered (or normal-ordered) terms; for instance
// :n_k n_l: + <n_k> :n_l: + <n_l> :n_k: + <n_k><n_l>
struct WickOrderedCollector : public mrock::symbolic_operators::AbstractCollector<WickOrderedTerm> {
    MROCK_FORWARD_CONSTRUCTORS(WickOrderedCollector, mrock::symbolic_operators::AbstractCollector<WickOrderedTerm>)

    MROCK_FORWARD_ASSIGNMENT(WickOrderedCollector, terms)

    /**
     * @brief Takes the Hermitian conjugate of the current expression in place
     * 
     * @return A reference to \c *this after the Hermitian conjugation 
     */
    WickOrderedCollector& hermitian_conjugate_inplace() noexcept;

    /**
     * @brief Takes the Hermitian conjugate of the current expression in a copy
     * 
     * @return A copy of \c *this after the Hermitian conjugation 
     */
    WickOrderedCollector hermitian_conjugate() const noexcept;

    // Sorts the terms vector according to a predefined predicate
    void sort();

    WickOrderedCollector& operator+=(const WickOrderedCollector& other);

    WickOrderedCollector& operator-=(const WickOrderedCollector& other);
};

/**
 * @brief Decomposes a Term into WickOrderedTerm objects using the provided WickOperatorTemplates.
 * 
 * @param term The Term to decompose.
 * @param templates The vector of WickOperatorTemplate objects to use for decomposition.
 * @return WickOrderedCollector The collector containing the resulting WickOrderedTerm objects.
 */
WickOrderedCollector wick_decompose(const Term& term, const std::vector<WickOperatorTemplate>& templates, const bool clean_up = false);

WickOrderedCollector wick_decompose(const TermCollector& terms, const std::vector<WickOperatorTemplate>& templates, const bool clean_up = false);

WickOrderedCollector multiply(const WickOrderedTerm& lhs, const WickOrderedTerm& rhs, std::vector<WickOperatorTemplate> templates);

WickOrderedCollector commutator(const WickOrderedTerm& lhs, const WickOrderedTerm& rhs, std::vector<WickOperatorTemplate> templates);

WickOrderedCollector commutator(const WickOrderedCollector& lhs,
                                const WickOrderedCollector& rhs,
                                std::vector<WickOperatorTemplate> templates);

void clean_wick_ordered_terms(WickOrderedCollector& terms,
    const std::vector<std::unique_ptr<WickSymmetry>>& symmetries = std::vector<std::unique_ptr<WickSymmetry>>{});

std::ostream& operator<<(std::ostream& os, const WickOrderedCollector& terms);

WickOrderedCollector operator+(WickOrderedCollector lhs, const WickOrderedCollector& rhs);

WickOrderedCollector operator-(WickOrderedCollector lhs, const WickOrderedCollector& rhs);

} // namespace mrock::symbolic_operators::experimental