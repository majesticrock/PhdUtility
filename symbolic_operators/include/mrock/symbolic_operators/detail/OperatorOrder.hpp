#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_DETAIL_OPERATOR_ORDER_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_DETAIL_OPERATOR_ORDER_HPP

#include "../AbstractTerm.hpp"
#include "../Operator.hpp"

#include <vector>

namespace mrock::symbolic_operators {

// Forward declarations
class Term;

namespace detail {
/**
 * @brief Tries to bring the momentum dependencies of the operators to a fixed notation
 * 
 * The passed container of operators should be the operators of the term that is being structured.
 * Allowing this additional parameter allows this function to work with more implementations
 * than just the standard \c Term class.
 * 
 * @tparam TermType The type of the term containing the operators.
 * @tparam OperatorContainer The type of the container holding the operators.
 * @param term The term containing the operators.
 * @param operators The container of operators to be structured.
 */
template <DerivedFromAbstractTerm TermType, class OperatorContainer>
void structure_momentum_dependencies_impl(TermType& term, OperatorContainer& operators)
{
// This function does not do anything to an identity term
    if (operators.empty())
        return;

    auto most_similar_daggered_operator = [&operators](const IndexWrapper& indices) {
        int best_match_count = 0;
        std::vector<Operator>::const_iterator best_it = operators.end();

        for (std::vector<Operator>::const_iterator it = operators.begin(); it != operators.end() && it->is_daggered;
             ++it) {
            // Penalize size mismatch
            int current_match_count = static_cast<int>(indices.size()) - static_cast<int>(it->indices.size());
            for (std::size_t i = 0U; i < it->indices.size() && i < indices.size(); ++i) {
                if (it->indices[i] == indices[i]) {
                    ++current_match_count;
                }
            }
            if (current_match_count > best_match_count) {
                best_match_count = current_match_count;
                best_it = it;
            }
        }
        return best_it;
    };

    MomentumSum::const_iterator sum_it = term.sums.momenta.begin();
    std::vector<MomentumSymbol::name_type> do_not_touch;
    do_not_touch.reserve(term.sums.momenta.size());

    for (auto op_it = operators.begin(); op_it < operators.end() && sum_it != term.sums.momenta.end(); ++op_it) {
        if (!op_it->is_daggered)
            break;

        try {
            const MomentumSymbol::name_type target = *sum_it;
            term.redistribute_momenta(op_it->momentum, target, do_not_touch);
            do_not_touch.push_back(target);
            ++sum_it;
        } catch (redistribution_error& e) {
        }
    }

    // Iterate backwards
    for (auto op_it = operators.rbegin(); op_it < operators.rend() && sum_it != term.sums.momenta.end(); ++op_it) {
        if (op_it->is_daggered)
            break;

        auto best_it = most_similar_daggered_operator(op_it->indices);
        if (best_it == operators.end())
            continue;
        try {
            const MomentumSymbol::name_type target = *sum_it;
            term.redistribute_momenta(op_it->momentum, target, do_not_touch);
            do_not_touch.push_back(target);
            ++sum_it;
            term.transform_momentum_sum(target, best_it->momentum + Momentum(PLACEHOLDER_SYMBOL), PLACEHOLDER_SYMBOL);
            term.rename_momenta(PLACEHOLDER_SYMBOL, target);
        } catch (redistribution_error& e) {
        }
    }

    for (auto coeff_it = term.coefficients.rbegin(); coeff_it != term.coefficients.rend(); ++coeff_it) {
        for (auto mom_it = coeff_it->momenta.rbegin(); mom_it != coeff_it->momenta.rend(); ++mom_it) {
            if (sum_it == term.sums.momenta.end())
                return;
            try {
                const MomentumSymbol::name_type target = *sum_it;
                term.redistribute_momenta(*mom_it, target, do_not_touch);
                do_not_touch.push_back(target);
                ++sum_it;
            } catch (redistribution_error& e) {
            }
        }
    }
}
} // namespace detail

/**
 * @brief Sorts a vector of operators by their indices while preserving the normal order.
 * 
 * @param operators The vector of operators to be sorted.
 * @return returns \c true if the sorting induced a fermionic sign change, \c false otherwise.
 */
bool sort_operators_by_indices(std::vector<Operator>& operators);

/**
 * @brief Tries to bring the momentum dependencies of the operators in the \c Term object to a fixed notation
 * 
 * @param term The term containing the operators.
 */
void structure_momentum_dependencies(Term& term);

} // namespace mrock::symbolic_operators
#endif // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_DETAIL_OPERATOR_ORDER_HPP