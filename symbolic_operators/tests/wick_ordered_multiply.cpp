#include <mrock/symbolic_operators/experimental/WickOrderedCollector.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace mrock::symbolic_operators;
using namespace mrock::symbolic_operators::experimental;

namespace {
bool contains_term(const std::vector<WickOrderedTerm>& terms, const WickOrderedTerm& expected) {
    return std::find(terms.begin(), terms.end(), expected) != terms.end();
}
}  // namespace

int main() {
    const Momentum k('k');
    const Momentum p('p');
    const Operator c_k_up_dagger(k, Index::SpinUp, true);
    const Operator c_k_down(k, Index::SpinDown, false);
    const Operator c_p_up_dagger(p, Index::SpinUp, true);
    const Operator c_p_up(p, Index::SpinUp, false);

    WickOrderedTerm lhs;
    lhs.multiplicity = 1;
    lhs.wick_expression.operators = {c_k_up_dagger, c_k_down};

    WickOrderedTerm rhs;
    rhs.multiplicity = 1;
    rhs.wick_expression.operators = {c_p_up_dagger, c_p_up};

    const std::vector<WickOperatorTemplate> templates(
        {WickOperatorTemplate{{Num_Comparison}, Momentum(), OperatorType::Number}});
    const WickOrderedCollector result = multiply(lhs, rhs, templates);

    WickOrderedTerm expected_uncontracted;
    expected_uncontracted.multiplicity = -1;
    expected_uncontracted.wick_expression.operators = {c_k_up_dagger, c_p_up_dagger, c_k_down, c_p_up};

    WickOrderedTerm expected_contraction;
    expected_contraction.multiplicity = -1;
    expected_contraction.delta_momenta.push_back({p, k});
    expected_contraction.operators.emplace_back(OperatorType::Number, false, p, Index::SpinUp);
    expected_contraction.wick_expression.operators = {c_p_up_dagger, c_k_down.with_momentum(p)};
    const WickOrderedCollector expected{expected_uncontracted, expected_contraction};

    std::cout << "% Product of two bilinear Wick-ordered strings\n\\begin{align*}\n\t" << lhs.wick_expression << " "
              << rhs.wick_expression << " = " << result << "\\end{align*}" << std::endl;
    std::cout << "% Analytical result, including the cross-contraction delta\n\\begin{align*}\n\t" << expected
              << "\\end{align*}" << std::endl;

    const bool has_expected_terms = result.size() == 2U && contains_term(result.terms, expected_uncontracted) &&
                                    contains_term(result.terms, expected_contraction);
    return has_expected_terms ? 0 : 1;
}