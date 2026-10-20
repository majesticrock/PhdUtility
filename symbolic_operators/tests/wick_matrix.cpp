#include <Eigen/Sparse>
#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/experimental/WickOrderedCollector.hpp>
#include <mrock/symbolic_operators/experimental/verify/Verifier.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numbers>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace mrock::symbolic_operators;
using namespace mrock::symbolic_operators::experimental;
using SparseMatrix = Eigen::SparseMatrix<double>;

using mrock::symbolic_operators::experimental::verify::L;
using mrock::symbolic_operators::experimental::verify::matrix_size;

int wrap_momentum(int momentum) {
    momentum %= L;
    return momentum < 0 ? momentum + L : momentum;
}

double cosine(int momentum) {
    return std::cos(std::numbers::pi * (2. * momentum / L - 1.));
}

TermCollector get_symbolic_H() {
    const Operator c_K_sigma(Momentum('K'), Index::Sigma, false);
    const Operator c_K_sigma_prime(Momentum('K'), Index::SigmaPrime, false);
    return {Term(1, Coefficient::RealInversionSymmetric("\\tilde{\\varepsilon}", MomentumList(Momentum('K'))),
                 SumContainer{MomentumSum{'K'}, IndexSum{Index::Sigma}}, {c_K_sigma.hermitian_conjugate(), c_K_sigma}),
            Term(1,
                 Coefficient::RealInteraction("U", MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}),
                                              IndexWrapper{Index::Sigma, Index::SigmaPrime},
                                              [](Coefficient& coeff) {
                                                  if (coeff.indices[0] > coeff.indices[1]) {
                                                      std::swap(coeff.indices[0], coeff.indices[1]);
                                                  }
                                              }),
                 SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
                 {c_K_sigma.hermitian_conjugate(), c_K_sigma_prime.with_momentum('P').hermitian_conjugate(),
                  c_K_sigma_prime.with_momentum(Momentum("P-Q")), c_K_sigma.with_momentum(Momentum("K+Q"))})};
}

TermCollector get_symbolic_eta() {
    const Operator c_K_sigma(Momentum('K'), Index::Sigma, false);
    const Operator c_K_sigma_prime(Momentum('K'), Index::SigmaPrime, false);
    return {Term(2,
                 Coefficient::RealInteraction("\\alpha", MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}),
                                              IndexWrapper{Index::Sigma, Index::SigmaPrime},
                                              [](Coefficient& coeff) {
                                                  if (coeff.indices[0] > coeff.indices[1]) {
                                                      std::swap(coeff.indices[0], coeff.indices[1]);
                                                  }
                                              }),
                 SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
                 {c_K_sigma.hermitian_conjugate(), c_K_sigma_prime.with_momentum('P').hermitian_conjugate(),
                  c_K_sigma_prime.with_momentum(Momentum("P-Q")), c_K_sigma.with_momentum(Momentum("K+Q"))})};
}

int main() {
    mrock::symbolic_operators::experimental::verify::Verifier verifier;

    auto explicit_total_number_operator = [&verifier]() {
        SparseMatrix result(matrix_size, matrix_size);
        for (int momentum = 0; momentum < L; ++momentum) {
            for (int spin = 0; spin < 2; ++spin) {
                const int mode = 2 * momentum + spin;
                result += verifier.operator_string({mode, mode}, {true, false});
            }
        }
        return result;
    };

    auto explicit_total_energy = [&verifier]() {
        constexpr double U = 1.;
        constexpr double V = 0.5;

        SparseMatrix result(matrix_size, matrix_size);
        for (int momentum = 0; momentum < L; ++momentum) {
            for (int spin = 0; spin < 2; ++spin) {
                const int mode = 2 * momentum + spin;
                result += (cosine(momentum) - 0.5 * U - 2. * V) * verifier.operator_string({mode, mode}, {true, false});
            }
        }
        return result;
    };

    // Verify the algebra
    {
        const SparseMatrix creation = verifier.operator_string({L / 2}, {true});
        const SparseMatrix annihilation = verifier.operator_string({L / 2}, {false});
        const SparseMatrix annihilation_transpose = annihilation.transpose();

        double diff = (creation - annihilation_transpose).norm();
        if (diff != 0.) {
            std::cerr << "c^+ != c     diff=" << diff << std::endl;
            return 1;
        }
        diff = (creation * annihilation - verifier.operator_string({L / 2, L / 2}, {true, false})).norm();
        if (diff != 0.) {
            std::cerr << "n != c^+ c     diff=" << diff << std::endl;
            return 1;
        }

        diff = (creation * creation).norm();
        if (diff != 0.) {
            std::cerr << "c c != 0     diff=" << diff << std::endl;
            return 1;
        }
        diff = (creation * annihilation + annihilation * creation).norm() - (matrix_size >> L);
        if (diff != 0.) {
            std::cerr << "{c, c^+} != 1     diff=" << diff << std::endl;
            return 1;
        }

        const SparseMatrix annihilation2 = verifier.operator_string({0}, {false});
        diff = (creation * annihilation2 + annihilation2 * creation).norm();
        if (diff != 0.) {
            std::cerr << "{c_i, c_j^+} != 0     diff=" << diff << std::endl;
            return 1;
        }
        diff = (annihilation * annihilation2 + annihilation2 * annihilation).norm();
        if (diff != 0.) {
            std::cerr << "{c_i^+, c_j^+} != 0     diff=" << diff << std::endl;
            return 1;
        }

        diff = verifier.operator_string({0, 0, 0, 0}, {true, true, false, false}).norm();
        if (diff != 0.) {
            std::cerr << "c^+ c^+ c c != 0     diff=" << diff << std::endl;
            return 1;
        }

        for (int i = 1; i < 2 * L; ++i) {
            const SparseMatrix some_interaction = verifier.operator_string({0, i, 2, 1}, {true, true, false, false});
            const SparseMatrix some_interaction_transpose =
                verifier.operator_string({1, 2, i, 0}, {true, true, false, false}).transpose();
            diff = (some_interaction - some_interaction_transpose).norm();
            if (diff != 0.) {
                std::cerr << "c^+ c^+ c c != (c^+ c^+ c c)^+     diff=" << diff << std::endl;
                return 1;
            }
        }
    }

    // Bilinear tests 1
    {
        for (int K = 0; K < L; ++K) {
            SparseMatrix n = verifier.operator_string({2 * K, 2 * K}, {true, false});
            SparseMatrix no = verifier.normal_ordered_operator_string({2 * K, 2 * K}, {true, false});
            SparseMatrix expected = n;
            expected -= verifier.occupations[K] * verifier.identity;

            if ((no - expected).norm() != 0.) {
                std::cerr << "bilinear normal_ordered_operator_string does not yield the expected result!" << "K=" << K
                          << ": " << (no - expected).norm() << std::endl;
                return 1;
            }
        }
    }

    // Bilinear tests 2
    {
        // Not using the check() lambda here, because we hard-code the expected results
        // This provides yet another independent check
        const Operator c_K_sigma(Momentum('K'), Index::Sigma, false);
        const Term total_number(1, SumContainer{MomentumSum{'K'}, IndexSum{Index::Sigma}},
                                {c_K_sigma.hermitian_conjugate(), c_K_sigma});
        WickOrderedCollector total_number_wick = wick_decompose(total_number, verifier.templates);
        clean_wick_ordered_terms(total_number_wick);
        if (!verifier.matrices_equal(explicit_total_number_operator(),
                                     verifier.symbolic_to_matrix(total_number_wick))) {
            std::cerr << "sum_(k,sigma) c_(k,sigma)^dagger c_(k,sigma) changed under Wick decomposition\n";
            std::cerr << "\\begin{align*}\n"
                      << total_number << "\\end{align*}\n\\begin{align*}" << total_number_wick << "\\end{align*}"
                      << std::endl;
            return 1;
        }

        const Term total_energy(1, Coefficient("\\tilde{\\varepsilon}", Momentum('K')),
                                SumContainer{MomentumSum{'K'}, IndexSum{Index::Sigma}},
                                {c_K_sigma.hermitian_conjugate(), c_K_sigma});
        WickOrderedCollector total_energy_wick = wick_decompose(total_energy, verifier.templates);
        clean_wick_ordered_terms(total_energy_wick);
        if (!verifier.matrices_equal(explicit_total_energy(), verifier.symbolic_to_matrix(total_energy_wick))) {
            std::cerr << "sum_(k,sigma) eps(k) c_(k,sigma)^dagger c_(k,sigma) changed under Wick decomposition\n";
            std::cerr << "\\begin{align*}\n"
                      << total_energy << "\\end{align*}\n\\begin{align*}" << total_energy_wick << "\\end{align*}"
                      << std::endl;
            return 1;
        }

        WickOrderedCollector total_number_backwards = {WickOrderedTerm()};
        total_number_backwards[0].multiplicity = 1;
        total_number_backwards[0].wick_expression.operators = {c_K_sigma.hermitian_conjugate(), c_K_sigma};
        total_number_backwards[0].sums = SumContainer{MomentumSum{'K'}, IndexSum{Index::Sigma}};
        SparseMatrix backwards_should_be = explicit_total_number_operator();
        double constant_to_substract = 0.;
        for (int K = 0; K < L; ++K) {
            constant_to_substract += 2 * verifier.occupations[K];
        }
        backwards_should_be -= constant_to_substract * verifier.identity;
        if (!verifier.matrices_equal(backwards_should_be, verifier.symbolic_to_matrix(total_number_backwards))) {
            std::cerr << "\\begin{align*}\n" << total_number_backwards << "\\end{align*}" << std::endl;
            std::cerr << "Matrices for total_number_backwards do not match!" << std::endl;
            return 1;
        }
    }

    // Quartic tests
    {
        {
            SparseMatrix no = verifier.normal_ordered_operator_string({0, 0, 0, 0}, {true, true, false, false});
            if ((no).norm() != 0.) {
                std::cerr << "quartic normal_ordered_operator_string does not yield the expected result!" << "i=" << 0
                          << ": " << (no).norm() << std::endl;
                return 1;
            }
        }
        for (int i = 1; i < 2 * L; ++i) {
            SparseMatrix raw = verifier.operator_string({0, i, i, 0}, {true, true, false, false}) -
                               verifier.occupations[0] * verifier.operator_string({i, i}, {true, false}) -
                               verifier.occupations[i / 2] * verifier.operator_string({0, 0}, {true, false}) +
                               verifier.occupations[0] * verifier.occupations[i / 2] * verifier.identity;
            SparseMatrix no = verifier.normal_ordered_operator_string({0, i, i, 0}, {true, true, false, false});

            if ((no - raw).norm() != 0.) {
                std::cerr << "quartic normal_ordered_operator_string does not yield the expected result!" << "i=" << i
                          << ": " << (no - raw).norm() << std::endl;
                return 1;
            }
        }
        for (int i = 0; i < 2 * L; ++i) {
            if (i == L)
                continue;
            SparseMatrix raw = verifier.operator_string({L, i, i, L}, {true, true, false, false}) -
                               verifier.occupations[L / 2] * verifier.operator_string({i, i}, {true, false}) -
                               verifier.occupations[i / 2] * verifier.operator_string({L, L}, {true, false}) +
                               verifier.occupations[L / 2] * verifier.occupations[i / 2] * verifier.identity;
            SparseMatrix no = verifier.normal_ordered_operator_string({L, i, i, L}, {true, true, false, false});

            if ((no - raw).norm() != 0.) {
                std::cerr << "quartic normal_ordered_operator_string does not yield the expected result!"
                          << "base=L | i=" << i << ": " << (no - raw).norm() << std::endl;
                return 1;
            }
        }
    }

    const TermCollector symbolic_H = get_symbolic_H();
    const TermCollector symbolic_eta = get_symbolic_eta();

    if (!verifier(symbolic_H, "H", true)) {
        return 1;
    }
    if (!verifier(symbolic_eta, "eta", false)) {
        return 1;
    }

    // sextic tests
    {
        // Degenerate case: all operators on the same mode.
        {
            SparseMatrix no =
                verifier.normal_ordered_operator_string({0, 0, 0, 0, 0, 0}, {true, true, true, false, false, false});

            if (no.norm() != 0.) {
                std::cerr << "sextic normal_ordered_operator_string does not yield the expected result!"
                          << " i=j=0: " << no.norm() << std::endl;
                return 1;
            }
        }

        // Generic case: three distinct modes 0,i,j.
        for (int i = 1; i < 2 * L; ++i) {
            for (int j = i + 1; j < 2 * L; ++j) {
                SparseMatrix raw =
                    verifier.operator_string({0, i, j, j, i, 0}, {true, true, true, false, false, false})

                    - verifier.occupations[0] * verifier.operator_string({i, j, j, i}, {true, true, false, false}) -
                    verifier.occupations[i / 2] * verifier.operator_string({0, j, j, 0}, {true, true, false, false}) -
                    verifier.occupations[j / 2] * verifier.operator_string({0, i, i, 0}, {true, true, false, false})

                    + verifier.occupations[0] * verifier.occupations[i / 2] *
                          verifier.operator_string({j, j}, {true, false}) +
                    verifier.occupations[0] * verifier.occupations[j / 2] *
                        verifier.operator_string({i, i}, {true, false}) +
                    verifier.occupations[i / 2] * verifier.occupations[j / 2] *
                        verifier.operator_string({0, 0}, {true, false})

                    - verifier.occupations[0] * verifier.occupations[i / 2] * verifier.occupations[j / 2] *
                          verifier.identity;

                SparseMatrix no = verifier.normal_ordered_operator_string({0, i, j, j, i, 0},
                                                                          {true, true, true, false, false, false});

                if ((no - raw).norm() != 0.) {
                    std::cerr << "sextic normal_ordered_operator_string does not yield the expected result!"
                              << " i=" << i << " j=" << j << ": " << (no - raw).norm() << std::endl;
                    return 1;
                }
            }
        }
    }

    TermCollector symbolic_commutator = commutator(symbolic_eta, symbolic_H);
    symbolic_commutator.clean_up();

    if (!verifier(symbolic_commutator, "[eta, H]", true)) {
        return 1;
    }

    return 0;
}
