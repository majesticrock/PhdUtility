#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_ABSTRACTCOLLECTOR_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_ABSTRACTCOLLECTOR_HPP

#include "detail/vector_macro.hpp"
#include "MomentumSymbol.hpp"

#include <vector>

namespace mrock::symbolic_operators {

template <class Collected>
struct AbstractCollector {
    std::vector<Collected> terms;  ///< The collected objects

    MROCK_VECTOR_WRAPPER_FILL_MEMBERS(Collected, terms);

    MROCK_FORWARD_CONSTRUCTORS(AbstractCollector, terms);

    MROCK_FORWARD_ASSIGNMENT(AbstractCollector, terms);

    virtual ~AbstractCollector() = default;

    /**
     * @brief Combines duplicate terms in \c *this.
     */
    void combine_duplicates();

    /**
     * @brief Renames momenta in \c *this.
     * 
     * @param what The momentum to rename.
     * @param to The new momentum.
     */
    void rename_momenta(const MomentumSymbol::name_type what, const MomentumSymbol::name_type to);
};

// Inline definitions
template <class Collected>
void AbstractCollector<Collected>::combine_duplicates() {
    // remove duplicates
    for (std::size_t i = 0U; i < terms.size(); i++) {
        for (std::size_t j = i + 1U; j < terms.size(); j++) {
            if (terms[i] == terms[j]) {
                terms[i].multiplicity += terms[j].multiplicity;
                terms.erase(terms.begin() + j);
                --i;
                break;
            }
        }
    }
    // removes any terms that have a 0 prefactor
    for (auto it = terms.begin(); it != terms.end();) {
        if (it->multiplicity == 0) {
            it = terms.erase(it);
        } else {
            ++it;
        }
    }
}

template <class Collected>
void AbstractCollector<Collected>::rename_momenta(const MomentumSymbol::name_type what, const MomentumSymbol::name_type to) {
    for (auto& t : terms) {
        t.rename_momenta(what, to);
    }
}

}  // namespace mrock::symbolic_operators
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_ABSTRACTCOLLECTOR_HPP