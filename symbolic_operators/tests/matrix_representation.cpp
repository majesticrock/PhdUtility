#include "IntMomentum.hpp"

#include <Eigen/Sparse>
#include <mrock/symbolic_operators/Commutation>

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace mrock::symbolic_operators;
using SparseMatrix = Eigen::SparseMatrix<double>;

// Consider a chain of length L
constexpr int L = 4;
// 2 spins per lattice site
constexpr int matrix_size = 1 << (2 * L);

const Operator c_K_sigma = Operator(Momentum('K'), Index::Sigma, false);
const Operator c_K_sigma_prime = Operator(Momentum('K'), Index::SigmaPrime, false);

// Filled in main
std::array<double, L> cosines;

Eigen::SparseMatrix<double> build_creation_operator(int site, int spin) {
    const int mode = 2 * site + spin;

    Eigen::SparseMatrix<double> C(matrix_size, matrix_size);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(matrix_size / 2);

    const uint32_t modeMask = 1u << mode;
    const uint32_t lowerMask = mode == 0 ? 0u : ((1u << mode) - 1u);

    for (uint32_t state = 0; state < matrix_size; ++state) {
        // Creation annihilates already-occupied states.
        if (state & modeMask)
            continue;

        uint32_t newState = state | modeMask;

        // Fermionic sign: (-1)^(number of occupied modes before 'mode').
        int parity = std::popcount(state & lowerMask) & 1;
        double sign = parity ? -1.0 : 1.0;

        triplets.emplace_back(static_cast<int>(newState), static_cast<int>(state), sign);
    }

    C.setFromTriplets(triplets.begin(), triplets.end());
    return C;
}

const Eigen::SparseMatrix<double>& creation_operator(int site, int spin) {
    static const std::array<Eigen::SparseMatrix<double>, 2 * L> cache = [] {
        std::array<Eigen::SparseMatrix<double>, 2 * L> ops;
        for (int s = 0; s < L; ++s)
            for (int sp = 0; sp < 2; ++sp)
                ops[2 * s + sp] = build_creation_operator(s, sp);
        return ops;
    }();

    return cache[2 * site + spin];
}

auto annihilation_operator(int lattice_site, int spin) {
    return creation_operator(lattice_site, spin).transpose();
}

SparseMatrix operator_string(const std::vector<int>& modes, const std::vector<bool>& creation) {
    if (modes.size() != creation.size()) {
        std::cerr << "operator_string: mode <-> operator mismatch!" << std::endl;
        abort();
    }

    SparseMatrix ret(matrix_size, matrix_size);
    if (modes.empty()) {
        ret.setIdentity();
        return ret;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(matrix_size);

    for (uint32_t input_state = 0; input_state < static_cast<uint32_t>(matrix_size); ++input_state) {
        uint32_t output_state = input_state;
        double sign = 1.;
        bool valid = true;

        // Matrix products act on a state from right to left.
        for (std::size_t i = modes.size(); i-- > 0;) {
            const uint32_t mode_mask = 1u << modes[i];
            const uint32_t lower_mask = mode_mask - 1u;
            const bool occupied = (output_state & mode_mask) != 0u;
            if (occupied == creation[i]) {
                valid = false;
                break;
            }
            if (creation[i]) {
                output_state |= mode_mask;
            } else {
                output_state &= ~mode_mask;
            }
            if ((std::popcount(output_state & lower_mask) & 1) != 0) {
                sign = -sign;
            }
        }

        if (valid) {
            triplets.emplace_back(static_cast<int>(output_state), static_cast<int>(input_state), sign);
        }
    }
    ret.setFromTriplets(triplets.begin(), triplets.end());
    return ret;
}

// Returns the product of creation and annihilation operators.
// The first half will be creation and the second half annihilation operators.
// For example c^dagger c^dagger c c
SparseMatrix product_string(const std::vector<int>& lattice_sites, const std::vector<int>& spins) {
    if (lattice_sites.empty() || lattice_sites.size() != spins.size() || lattice_sites.size() % 2 != 0) {
        std::cerr << "product_string: lattice_site <-> spin mismatch!" << std::endl;
        abort();
    }

    const std::size_t creation_count = lattice_sites.size() / 2;
    std::vector<int> modes;
    modes.reserve(lattice_sites.size());
    std::vector<bool> creation(lattice_sites.size());
    for (std::size_t i = 0; i < lattice_sites.size(); ++i) {
        modes.push_back(2 * lattice_sites[i] + spins[i]);
        creation[i] = i < creation_count;
    }
    return operator_string(modes, creation);
}

TermCollector get_symbolic_H() {
    const Term H_kin_symbolic(1,
                              Coefficient::RealInversionSymmetric("\\tilde{\\varepsilon}", MomentumList(Momentum('K'))),
                              SumContainer{MomentumSum{'K'}, IndexSum{Index::Sigma}},
                              std::vector<Operator>({c_K_sigma.hermitian_conjugate(), c_K_sigma}));

    const Term H_int_symbolic(
        1,
        Coefficient::RealInteraction("U", MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}),
                                     IndexWrapper{Index::Sigma, Index::SigmaPrime}),
        SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
        std::vector<Operator>(
            {c_K_sigma.hermitian_conjugate(), c_K_sigma_prime.with_momentum('P').hermitian_conjugate(),
             c_K_sigma_prime.with_momentum(Momentum("P-Q")), c_K_sigma.with_momentum(Momentum("K+Q"))}));

    return {H_kin_symbolic, H_int_symbolic};
}

TermCollector get_symbolic_eta() {
    const Term eta(2,
                   Coefficient::RealInteraction("\\alpha", MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}),
                                                IndexWrapper{Index::Sigma, Index::SigmaPrime}),
                   SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
                   std::vector<Operator>(
                       {c_K_sigma.hermitian_conjugate(), c_K_sigma_prime.with_momentum('P').hermitian_conjugate(),
                        c_K_sigma_prime.with_momentum(Momentum("P-Q")), c_K_sigma.with_momentum(Momentum("K+Q"))}));

    return {eta};
}

SparseMatrix get_matrix_H() {
    SparseMatrix H(matrix_size, matrix_size);
    for (int k = 0; k < L; ++k) {
        H += cosines[k] * (product_string({k, k}, {0, 0}) + product_string({k, k}, {1, 1}));
    }
    for (int k = 0; k < L; ++k) {
        for (int p = 0; p < L; ++p) {
            for (int q = 0; q < L; ++q) {
                const int k_q = IntMomentum<L>(k) + IntMomentum<L>(q);
                const int p_q = IntMomentum<L>(p) - IntMomentum<L>(q);
                H += (cosines[k_q] * cosines[p_q] + cosines[k] * cosines[q]) *
                     (product_string({k, p, p_q, k_q}, {0, 1, 1, 0}) + product_string({k, p, p_q, k_q}, {1, 0, 0, 1}));
                H += cosines[q] * (cosines[k_q] * cosines[p_q] + cosines[k] * cosines[q]) *
                     (product_string({k, p, p_q, k_q}, {0, 0, 0, 0}) + product_string({k, p, p_q, k_q}, {1, 1, 1, 1}));
            }
        }
    }
    return H;
}

SparseMatrix get_matrix_eta() {
    SparseMatrix eta(matrix_size, matrix_size);
    for (int k = 0; k < L; ++k) {
        for (int p = 0; p < L; ++p) {
            for (int q = 0; q < L; ++q) {
                const int k_q = IntMomentum<L>(k) + IntMomentum<L>(q);
                const int p_q = IntMomentum<L>(p) - IntMomentum<L>(q);

                eta += (cosines[k_q] * cosines[p_q] + cosines[k] * cosines[q]) 
                    * (cosines[k] + cosines[p] - cosines[k_q] - cosines[p_q]) *
                    (product_string({k, p, p_q, k_q}, {0, 1, 1, 0}) + product_string({k, p, p_q, k_q}, {1, 0, 0, 1}));
                eta += cosines[q] * (cosines[k_q] * cosines[p_q] + cosines[k] * cosines[q])
                * cosines[k_q] * cosines[p_q] * (cosines[k] + cosines[p] - cosines[k_q] - cosines[p_q]) *
                    (product_string({k, p, p_q, k_q}, {0, 0, 0, 0}) + product_string({k, p, p_q, k_q}, {1, 1, 1, 1}));
            }
        }
    }
    return eta;
}

SparseMatrix symbolic_to_matrix(const TermCollector& terms) {
    SparseMatrix result(matrix_size, matrix_size);

    for (const auto& term : terms) {
        std::array<IntMomentum<L>, 256> momentum_values{};
        std::array<int, 256> index_values{};

        const auto evaluate_momentum = [&momentum_values](const Momentum& momentum) {
            IntMomentum<L> value(momentum.add_PI ? 0 : L/2);
            for (const auto& symbol : momentum.momentum_list) {
                value += symbol.factor * momentum_values[static_cast<unsigned char>(symbol.name)];
            }
            return value;
        };
        const auto evaluate_index = [&index_values](const Index& index) {
            return index_values[static_cast<unsigned char>(index)];
        };

        const auto evaluate_coefficient = [&evaluate_momentum, &evaluate_index](const Coefficient& coefficient) {
            if (coefficient.name == "\\tilde{\\varepsilon}") {
                return cosines[evaluate_momentum(coefficient.momenta.front())];
            } 

            double value = cosines[evaluate_momentum(coefficient.momenta[0] + coefficient.momenta[2])] *
                           cosines[evaluate_momentum(coefficient.momenta[1] - coefficient.momenta[2])]
                        + cosines[evaluate_momentum(coefficient.momenta[0])] *
                          cosines[evaluate_momentum(coefficient.momenta[1])];
            if (evaluate_index(coefficient.indices[0]) != evaluate_index(coefficient.indices[1])) {
                value *= cosines[evaluate_momentum(coefficient.momenta[2])];
            }

            if (coefficient.name == "U") {
                return value;
                
            } 
            else if (coefficient.name == "\\alpha") {
                value *= (cosines[evaluate_momentum(coefficient.momenta[0])] +
                          cosines[evaluate_momentum(coefficient.momenta[1])] -
                          cosines[evaluate_momentum(coefficient.momenta[0] + coefficient.momenta[2])] -
                          cosines[evaluate_momentum(coefficient.momenta[1] - coefficient.momenta[2])]);
            } 
            else {
                throw std::runtime_error("No matrix representation for coefficient " + coefficient.name);
            }
            return value;
        };

        const auto evaluate_term = [&]() {
            double coefficient = static_cast<double>(term.multiplicity);
            for (const auto& term_coefficient : term.coefficients) {
                coefficient *= evaluate_coefficient(term_coefficient);
            }

            for (const auto& delta : term.delta_momenta) {
                if (evaluate_momentum(delta.first) != evaluate_momentum(delta.second)) {
                    return;
                }
            }
            for (const auto& delta : term.delta_indices) {
                const auto index_value = [&index_values](const Index index) {
                    return is_mutable(index) ? index_values[static_cast<unsigned char>(index)]
                                             : (index == Index::SpinDown ? 1 : 0);
                };
                if (index_value(delta.first) != index_value(delta.second)) {
                    return;
                }
            }

            std::vector<int> modes;
            std::vector<bool> creation;
            modes.reserve(term.operators.size());
            creation.reserve(term.operators.size());
            for (const auto& op : term.operators) {
                const int momentum = evaluate_momentum(op.momentum);
                const int spin = is_mutable(op.first_index())
                                     ? index_values[static_cast<unsigned char>(op.first_index())]
                                     : (op.first_index() == Index::SpinDown ? 1 : 0);
                modes.push_back(2 * momentum + spin);
                creation.push_back(op.is_daggered);
            }
            result += coefficient * operator_string(modes, creation);
        };

        // std function instead of lambda so that it can be called recursively
        std::function<void(std::size_t)> assign_spins;
        assign_spins = [&](const std::size_t sum_index) {
            if (sum_index == term.sums.spins.size()) {
                evaluate_term();
                return;
            }
            index_values[static_cast<unsigned char>(term.sums.spins[sum_index])] = 0;
            assign_spins(sum_index + 1U);
            index_values[static_cast<unsigned char>(term.sums.spins[sum_index])] = 1;
            assign_spins(sum_index + 1U);
        };

        std::function<void(std::size_t)> assign_momenta;
        assign_momenta = [&](const std::size_t sum_index) {
            if (sum_index == term.sums.momenta.size()) {
                assign_spins(0U);
                return;
            }
            const auto name = static_cast<unsigned char>(term.sums.momenta[sum_index]);
            for (int value = 0; value < L; ++value) {
                momentum_values[name] = IntMomentum<L>(value);
                assign_momenta(sum_index + 1U);
            }
        };
        assign_momenta(0U);
    }
    return result;
}

int main() {
    {
        // Verify the algebra
        const SparseMatrix creation = creation_operator(L / 2, 0);
        const SparseMatrix annihilation = annihilation_operator(L / 2, 0);
        const SparseMatrix annihilation_transpose = annihilation_operator(L / 2, 0).transpose();

        double diff = (creation - annihilation_transpose).norm();
        if (diff != 0.) {
            std::cerr << "c^+ != c     diff=" << diff << std::endl;
            return 1;
        }
        diff = (creation_operator(L / 2, 0) * annihilation_operator(L / 2, 0) - product_string({L / 2, L / 2}, {0, 0}))
                   .norm();
        if (diff != 0.) {
            std::cerr << "n != c^+ c     diff=" << diff << std::endl;
            return 1;
        }

        diff = (creation * creation + creation * creation).norm();
        if (diff != 0.) {
            std::cerr << "{c, c} != 0     diff=" << diff << std::endl;
            return 1;
        }
        diff = (creation * annihilation + annihilation * creation).norm() - (matrix_size >> L);
        if (diff != 0.) {
            std::cerr << "{c, c^+} != 1     diff=" << diff << std::endl;
            return 1;
        }

        const SparseMatrix annihilation2 = annihilation_operator(0, 1);
        diff = (creation * annihilation2 + annihilation2 * creation).norm();
        if (diff != 0.) {
            std::cerr << "{c_i, c_j^+} != 0     diff=" << diff << std::endl;
            return 1;
        }
    }

    for (int i = 0; i < L; ++i) {
        cosines[i] = std::cos(std::numbers::pi * (2 * i - L));
    }

    const TermCollector symbolic_H = get_symbolic_H();
    const SparseMatrix matrix_H = get_matrix_H();
    if ((symbolic_to_matrix(symbolic_H) - matrix_H).norm() > 1.e-12) {
        std::cerr << "symbolic_H != matrix_H" << std::endl;
        return 1;
    }

    const TermCollector symbolic_eta = get_symbolic_eta();
    const SparseMatrix matrix_eta = get_matrix_eta();
    if ((symbolic_to_matrix(symbolic_eta) - matrix_eta).norm() > 1.e-12) {
        std::cerr << "symbolic_eta != matrix_eta" << std::endl;
        return 1;
    }

    TermCollector symbolic_commutator = commutator(symbolic_eta, symbolic_H);
    symbolic_commutator.clean_up();

    const SparseMatrix matrix_commutator = matrix_eta * matrix_H - matrix_H * matrix_eta;
    if ((symbolic_to_matrix(symbolic_commutator) - matrix_commutator).norm() > 1.e-12) {
        std::cerr << "symbolic_commutator != matrix_commutator" << std::endl;
        return 1;
    }

    return 0;
}