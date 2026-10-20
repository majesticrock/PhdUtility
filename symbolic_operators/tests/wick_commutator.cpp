#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/experimental/WickOrderedCollector.hpp>
#include <mrock/symbolic_operators/experimental/verify/Verifier.hpp>

#include <array>
#include <iostream>
#include <vector>

using namespace mrock::symbolic_operators;
using namespace mrock::symbolic_operators::experimental;

using mrock::symbolic_operators::experimental::verify::Verifier;

template <class LhsType, class RhsType>
void print_commutator(const LhsType& lhs, const RhsType& rhs, const WickOrderedCollector& result) {
    std::cout << "\\begin{align*}\n\t[" << lhs << ", " << rhs << "] = " << result << "\\end{align*}\n";
}

WickOrderedTerm make_bilinear(const Momentum& creation_momentum, const Momentum& annihilation_momentum) {
    WickOrderedTerm result;
    result.multiplicity = 1;
    result.wick_expression.operators = {Operator(creation_momentum, Index::SpinUp, true),
                                        Operator(annihilation_momentum, Index::SpinUp, false)};
    return result;
}

WickOrderedTerm make_expectation_bilinear(const Momentum& expectation_momentum,
                                          const Momentum& creation_momentum,
                                          const Momentum& annihilation_momentum) {
    WickOrderedTerm result = make_bilinear(creation_momentum, annihilation_momentum);
    result.operators.emplace_back(OperatorType::Number, false, expectation_momentum, Index::SpinUp);
    return result;
}

WickOrderedTerm make_summed_quartic(const Momentum& creation_momentum, const Momentum& annihilation_momentum) {
    WickOrderedTerm result;
    result.multiplicity = 1;
    const Momentum summed_momentum('K');
    result.wick_expression.operators = {
        Operator(creation_momentum, Index::SpinUp, true), Operator(summed_momentum, Index::SpinDown, true),
        Operator(summed_momentum, Index::SpinDown, false), Operator(annihilation_momentum, Index::SpinUp, false)};
    result.sums = SumContainer{MomentumSum{'K'}, IndexSum{}};
    return result;
}

bool verify_commutator(const WickOrderedTerm& lhs,
                       const WickOrderedTerm& rhs,
                       const char* name,
                       const Verifier& verifier) {
    const WickOrderedCollector lhs_collector{lhs};
    const WickOrderedCollector rhs_collector{rhs};
    const Verifier::SparseMatrix expected =
        verifier.symbolic_to_matrix(lhs_collector) * verifier.symbolic_to_matrix(rhs_collector) -
        verifier.symbolic_to_matrix(rhs_collector) * verifier.symbolic_to_matrix(lhs_collector);

    WickOrderedCollector result = commutator(lhs, rhs, verifier.templates);
    clean_wick_ordered_terms(result, verifier.symmetries);
    print_commutator(lhs, rhs, result);
    return verifier(result, name, expected);
}

bool verify_commutator(const WickOrderedCollector& lhs,
                       const WickOrderedCollector& rhs,
                       const char* name,
                       const Verifier& verifier) {
    const Verifier::SparseMatrix lhs_matrix = verifier.symbolic_to_matrix(lhs);
    const Verifier::SparseMatrix rhs_matrix = verifier.symbolic_to_matrix(rhs);
    const Verifier::SparseMatrix expected = lhs_matrix * rhs_matrix - rhs_matrix * lhs_matrix;

    WickOrderedCollector result = commutator(lhs, rhs, verifier.templates);
    clean_wick_ordered_terms(result, verifier.symmetries);
    print_commutator(lhs, rhs, result);
    return verifier(result, name, expected);
}

int main() {
    const Verifier verifier;
    const Momentum k('k');
    Momentum p('p');
    p.add_PI = true;

    const WickOrderedTerm number_k = make_bilinear(k, k);
    const WickOrderedTerm number_p = make_bilinear(p, p);
    const WickOrderedTerm hop_kp = make_bilinear(k, p);
    const WickOrderedTerm hop_pk = make_bilinear(p, k);
    const WickOrderedTerm occupation_hop_kp = make_expectation_bilinear(k, k, p);
    const WickOrderedTerm occupation_hop_pk = make_expectation_bilinear(p, p, k);
    const WickOrderedTerm summed_quartic_kp = make_summed_quartic(k, p);
    const WickOrderedTerm summed_quartic_pk = make_summed_quartic(p, k);

    const Verifier fractional_occupation_verifier(std::array<double, 4>{0.2, 0.4, 0.6, 0.8});
    WickOrderedTerm number_expectation;
    number_expectation.multiplicity = 1;
    number_expectation.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);
    WickOrderedTerm hole_number_expectation;
    hole_number_expectation.multiplicity = 1;
    hole_number_expectation.operators.emplace_back(OperatorType::HoleNumber, false, k, Index::SpinUp);
    const auto number_matrix =
        fractional_occupation_verifier.symbolic_to_matrix(WickOrderedCollector{number_expectation});
    const auto hole_number_matrix =
        fractional_occupation_verifier.symbolic_to_matrix(WickOrderedCollector{hole_number_expectation});
    if (!fractional_occupation_verifier.matrices_equal(hole_number_matrix,
                                                       fractional_occupation_verifier.identity - number_matrix)) {
        std::cerr << "HoleNumber expectation did not evaluate to 1 - Number expectation: difference norm = "
                  << (hole_number_matrix - (fractional_occupation_verifier.identity - number_matrix)).norm()
                  << ", number = " << number_matrix.coeff(0, 0) << ", hole number = " << hole_number_matrix.coeff(0, 0)
                  << '\n';
        return 1;
    }

    if (!verify_commutator(number_k, number_p, "[n_k, n_p]", verifier)) {
        return 1;
    }
    if (!verify_commutator(number_k, hop_kp, "[n_k, c_k^dagger c_p]", verifier)) {
        return 1;
    }
    if (!verify_commutator(hop_kp, hop_pk, "[c_k^dagger c_p, c_p^dagger c_k]", verifier)) {
        return 1;
    }
    if (!verify_commutator(occupation_hop_kp, hop_pk, "[<n_k> :c_k^dagger c_p:, c_p^dagger c_k]", verifier)) {
        return 1;
    }
    if (!verify_commutator(occupation_hop_kp, occupation_hop_pk, "[<n_k> :c_k^dagger c_p:, <n_p> :c_p^dagger c_k:]",
                           verifier)) {
        return 1;
    }
    if (!verify_commutator(number_k, summed_quartic_kp, "[n_k, sum_p c_k^dagger c_p^dagger_down c_p_down c_q]",
                           verifier)) {
        return 1;
    }
    if (!verify_commutator(summed_quartic_kp, summed_quartic_pk,
                           "[sum_p c_k^dagger c_p^dagger_down c_p_down c_q, "
                           "sum_p c_q^dagger c_p^dagger_down c_p_down c_k]",
                           verifier)) {
        return 1;
    }

    const WickOrderedCollector lhs{number_k, hop_kp};
    const WickOrderedCollector rhs{number_p, hop_pk};
    if (!verify_commutator(lhs, rhs, "[n_k + c_k^dagger c_p, n_p + c_p^dagger c_k]", verifier)) {
        return 1;
    }
    return 0;
}
