#include <mrock/symbolic_operators/experimental/WickOrderedCollector.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace mrock::symbolic_operators;
using namespace mrock::symbolic_operators::experimental;

namespace {
bool contains_term(const WickOrderedCollector& terms, const WickOrderedTerm& expected) {
    return std::find(terms.begin(), terms.end(), expected) != terms.end();
}

WickOrderedCollector decompose_and_clean(const TermCollector& terms,
                                         const std::vector<WickOperatorTemplate>& templates) {
    WickOrderedCollector result = wick_decompose(terms, templates);
    clean_wick_ordered_terms(result);
    return result;
}
}  // namespace

int main() {
    const Momentum k('k');
    const IndexWrapper spin_and_orbital_a{ Index::SpinUp, Index::TypeA };
    const IndexWrapper spin_and_orbital_b{ Index::SpinUp, Index::TypeB };
    const Operator c_k_a_dagger(k, spin_and_orbital_a, true);
    const Operator c_k_a(k, spin_and_orbital_a, false);
    const Operator c_k_b(k, spin_and_orbital_b, false);

    const TermCollector different_orbitals({ Term(1, std::vector<Operator>{ c_k_a_dagger, c_k_b }) });
    const TermCollector matching_orbitals({ Term(1, std::vector<Operator>{ c_k_a_dagger, c_k_a }) });
    const TermCollector reversed_matching_orbitals({ Term(1, std::vector<Operator>{ c_k_a, c_k_a_dagger }) });

    const std::vector<WickOperatorTemplate> orbital_independent_template({
        WickOperatorTemplate{{ Num_Comparison }, Momentum(), OperatorType::Number}
    });
    const std::vector<WickOperatorTemplate> orbital_matching_template({
        WickOperatorTemplate{{ Num_Comparison, Num_Comparison }, Momentum(), OperatorType::Number}
    });

    const WickOrderedCollector independent_result =
        decompose_and_clean(different_orbitals, orbital_independent_template);
    const WickOrderedCollector mismatch_result =
        decompose_and_clean(different_orbitals, orbital_matching_template);
    const WickOrderedCollector matching_result =
        decompose_and_clean(matching_orbitals, orbital_matching_template);
    const WickOrderedCollector reversed_result =
        decompose_and_clean(reversed_matching_orbitals, orbital_matching_template);

    std::cout << "% Experimental Wick decomposition, orbital ignored\n\\begin{align*}\n\t"
              << different_orbitals << " = " << independent_result
              << "\\end{align*}" << std::endl;
    std::cout << "% Experimental Wick decomposition, orbital required\n\\begin{align*}\n\t"
              << different_orbitals << " = " << mismatch_result << "\\\\\n"
              << matching_orbitals << " = " << matching_result
              << "\\end{align*}" << std::endl;
    std::cout << "% Reversed order uses the hole expectation, with h=1-n\n\\begin{align*}\n\t"
              << reversed_matching_orbitals << " = " << reversed_result
              << "\\end{align*}" << std::endl;

    WickOrderedTerm uncontracted_different;
    uncontracted_different.multiplicity = 1;
    uncontracted_different.wick_expression.operators = { c_k_a_dagger, c_k_b };
    WickOrderedTerm spin_only_contraction;
    spin_only_contraction.multiplicity = 1;
    spin_only_contraction.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);
    WickOrderedTerm uncontracted_matching;
    uncontracted_matching.multiplicity = 1;
    uncontracted_matching.wick_expression.operators = { c_k_a_dagger, c_k_a };
    WickOrderedTerm orbital_contraction;
    orbital_contraction.multiplicity = 1;
    orbital_contraction.operators.emplace_back(OperatorType::Number, false, k, spin_and_orbital_a);
    WickOrderedTerm reversed_residual;
    reversed_residual.multiplicity = -1;
    reversed_residual.wick_expression.operators = { c_k_a_dagger, c_k_a };
    WickOrderedTerm hole_contraction;
    hole_contraction.multiplicity = 1;
    hole_contraction.operators.emplace_back(OperatorType::HoleNumber, false, k, spin_and_orbital_a);

    if (independent_result.size() != 2U || !contains_term(independent_result, uncontracted_different) ||
        !contains_term(independent_result, spin_only_contraction) || mismatch_result.size() != 1U ||
        !contains_term(mismatch_result, uncontracted_different) || matching_result.size() != 2U ||
        !contains_term(matching_result, uncontracted_matching) || !contains_term(matching_result, orbital_contraction) ||
        reversed_result.size() != 2U || !contains_term(reversed_result, reversed_residual) ||
        !contains_term(reversed_result, hole_contraction)) {
        return 1;
    }
    return 0;
}