#include <mrock/symbolic_operators/ExpectationValues>

#include <iostream>
#include <vector>

using namespace mrock::symbolic_operators;

int main() {
    const Momentum k('k');
    const IndexWrapper spin_and_orbital_a{Index::SpinUp, Index::TypeA};
    const IndexWrapper spin_and_orbital_b{Index::SpinUp, Index::TypeB};
    const Operator c_k_a_dagger(k, spin_and_orbital_a, true);
    const Operator c_k_a(k, spin_and_orbital_a, false);
    const Operator c_k_b(k, spin_and_orbital_b, false);

    const TermCollector different_orbitals({Term(1, std::vector<Operator>{c_k_a_dagger, c_k_b})});
    const TermCollector matching_orbitals({Term(1, std::vector<Operator>{c_k_a_dagger, c_k_a})});

    const std::vector<WickOperatorTemplate> orbital_independent_template(
        {WickOperatorTemplate{{Num_Comparison}, Momentum(), OperatorType::Number}});
    const std::vector<WickOperatorTemplate> orbital_matching_template(
        {WickOperatorTemplate{{Num_Comparison, Num_Comparison}, Momentum(), OperatorType::Number}});

    const WickTermCollector independent_result = wicks_theorem(different_orbitals, orbital_independent_template, true);
    const WickTermCollector mismatch_result = wicks_theorem(different_orbitals, orbital_matching_template, true);
    const WickTermCollector matching_result = wicks_theorem(matching_orbitals, orbital_matching_template, true);

    std::cout << "% Orbital-independent expectation value\n\\begin{align*}\n\t" << "\\langle" << different_orbitals
              << "\\rangle = " << independent_result << "\\end{align*}" << std::endl;
    std::cout << "% Orbital indices must match\n\\begin{align*}\n\t" << "\\langle" << different_orbitals
              << "\\rangle = " << mismatch_result << "\\\\\n"
              << "\\langle" << matching_orbitals << "\\rangle = " << matching_result << "\\end{align*}" << std::endl;

    WickTerm expected_independent;
    expected_independent.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);
    WickTerm expected_matching;
    expected_matching.operators.emplace_back(OperatorType::Number, false, k, spin_and_orbital_a);
    const WickTermCollector expected_independent_result({expected_independent});
    const WickTermCollector expected_matching_result({expected_matching});

    if (independent_result.terms != expected_independent_result.terms || !mismatch_result.empty() ||
        matching_result.terms != expected_matching_result.terms) {
        return 1;
    }
    return 0;
}