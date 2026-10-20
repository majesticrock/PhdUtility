#pragma once

#include <Eigen/Sparse>
#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/IntMomentum.hpp>
#include <mrock/symbolic_operators/experimental/WickOrderedCollector.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace mrock::symbolic_operators::experimental::verify {
constexpr int L = 4;
constexpr int matrix_size = 1 << (2 * L);

class Verifier {
public:
    using SparseMatrix = Eigen::SparseMatrix<double>;

    const std::array<double, L> occupations;
    const std::array<double, L> cosines;

    const std::vector<mrock::symbolic_operators::WickOperatorTemplate> templates;
    const std::vector<std::unique_ptr<mrock::symbolic_operators::WickSymmetry>> symmetries;

    const SparseMatrix identity;

    Verifier();
    Verifier(const std::array<double, L>& _occupations);

    // We assume that creation is ordered: all true then all false (that means strings look like c^+ c^+ c^+ c c c c c)
    // The modes vector contains the numbered indizes
    // Example: c_1^+ c_4^+ c_3 c_1 would be
    // modes = {1,4,3,1}        creation = {true,true,false,false}
    SparseMatrix operator_string(const std::vector<int>& modes, const std::vector<bool>& creation) const;

    // We assume that creation is ordered: all true then all false (that means strings look like c^+ c^+ c^+ c c c c c)
    // The modes vector contains the numbered indizes
    // Example: c_1^+ c_4^+ c_3 c_1 would be
    // modes = {1,4,3,1}        creation = {true,true,false,false}
    SparseMatrix normal_ordered_operator_string(const std::vector<int>& modes, const std::vector<bool>& creation) const;

    SparseMatrix symbolic_to_matrix(const mrock::symbolic_operators::TermCollector& terms) const;
    SparseMatrix symbolic_to_matrix(const mrock::symbolic_operators::experimental::WickOrderedCollector& terms) const;

    bool is_particle_hole_invariant(const mrock::symbolic_operators::TermCollector& terms, const char* name) const;
    bool is_particle_hole_invariant(const mrock::symbolic_operators::experimental::WickOrderedCollector& terms,
                                    const char* name) const;

    // Computes (lhs - rhs).norm() and checks whether it is smaller than 1e-12
    bool matrices_equal(const SparseMatrix& lhs, const SparseMatrix& rhs) const;

    // Applies Wick's theorem to expression and checks whether the matrix representations before
    // and after are identical. If so, this function returns true.
    // If not, it prints debug information and returns false.
    bool operator()(const mrock::symbolic_operators::TermCollector& expression, const char* name, bool hermitian) const;

    // Applies Wick's theorem to expression and checks whether the matrix representations before
    // and after are identical. If so, this function returns true.
    // If not, it prints debug information and returns false.
    // Saves the resulting matrix representation to save_matrix
    bool operator()(const mrock::symbolic_operators::TermCollector& expression,
                    const char* name,
                    bool hermitian,
                    SparseMatrix& save_matrix) const;

    // Compares the matrix representation against compare
    bool operator()(const mrock::symbolic_operators::experimental::WickOrderedCollector& expression,
                    const char* name,
                    const SparseMatrix compare) const;

private:
    mutable std::array<IntMomentum<L>, 256> momenta{};
    mutable std::array<int, 256> indices{};

    // Private functions
    static std::array<double, 4> make_cosines();
    static std::array<double, 4> make_default_occupations();
    static SparseMatrix make_identity();

    IntMomentum<L> momentum_lookup(const mrock::symbolic_operators::Momentum& momentum) const noexcept;
    int index_lookup(const mrock::symbolic_operators::Index& index) const noexcept;

    double coefficient_value(const mrock::symbolic_operators::Coefficient& coefficient) const;

    template <class CollectorType>
    SparseMatrix symbolic_to_matrix_impl(const CollectorType& terms) const;
};

// Template defintion
template <class CollectorType>
Verifier::SparseMatrix Verifier::symbolic_to_matrix_impl(const CollectorType& terms) const {
    constexpr bool is_wick =
        std::is_same_v<CollectorType, mrock::symbolic_operators::experimental::WickOrderedCollector>;
    assert(occupations.size() == (std::size_t)L);

    SparseMatrix result(matrix_size, matrix_size);
    for (const auto& term : terms) {
        momenta.fill(IntMomentum<L>(0));
        indices.fill(0);
        indices[static_cast<unsigned char>(Index::SpinUp)] = 0;
        indices[static_cast<unsigned char>(Index::SpinDown)] = 1;

        std::function<void(std::size_t)> assign_spins;

        const auto evaluate = [&]() {
            for (const auto& delta : term.delta_momenta) {
                if (momentum_lookup(delta.first) != momentum_lookup(delta.second)) {
                    return;
                }
            }
            for (const auto& delta : term.delta_indices) {
                if (index_lookup(delta.first) != index_lookup(delta.second)) {
                    return;
                }
            }

            double prefactor = static_cast<double>(term.multiplicity);
            for (const auto& coefficient : term.coefficients) {
                prefactor *= coefficient_value(coefficient);
            }

            if constexpr (is_wick) {
                for (const auto& expec : term.operators) {
                    if (expec.type == mrock::symbolic_operators::OperatorType::Number) {
                        prefactor *= occupations[momentum_lookup(expec.momentum)];
                    } else if (expec.type == mrock::symbolic_operators::OperatorType::HoleNumber) {
                        prefactor *= (1. - occupations[momentum_lookup(expec.momentum)]);
                    } else {
                        throw std::runtime_error("Expected only Number- or HoleNumber-type expectation values!");
                    }
                }
            }

            std::vector<int> modes;
            std::vector<bool> creation;

            const std::vector<mrock::symbolic_operators::Operator>& symbolic_operator_string =
                ([&]() -> const std::vector<mrock::symbolic_operators::Operator>& {
                    if constexpr (is_wick) {
                        return term.wick_expression.operators;
                    } else {
                        return term.operators;
                    }
                })();
            for (const auto& op : symbolic_operator_string) {
                const int spin = indices[index_lookup(op.first_index())];

                modes.push_back(2 * momentum_lookup(op.momentum).i + spin);
                creation.push_back(op.is_daggered);
            }
            if (prefactor != 0.) {
                if constexpr (is_wick) {
                    result += prefactor * normal_ordered_operator_string(modes, creation);
                } else {
                    result += prefactor * operator_string(modes, creation);
                }
            }
        };

        assign_spins = [&](std::size_t n) {
            if (n == term.sums.spins.size()) {
                evaluate();
                return;
            }
            indices[static_cast<unsigned char>(term.sums.spins[n])] = 0;
            assign_spins(n + 1);
            indices[static_cast<unsigned char>(term.sums.spins[n])] = 1;
            assign_spins(n + 1);
        };

        std::function<void(std::size_t)> assign_momenta;
        assign_momenta = [&](std::size_t n) {
            if (n == term.sums.momenta.size()) {
                assign_spins(0);
                return;
            }
            const auto name = static_cast<unsigned char>(term.sums.momenta[n]);
            for (int value = 0; value < L; ++value) {
                momenta[name].i = value;
                assign_momenta(n + 1);
            }
        };

        assign_momenta(0);
    }
    return result;
}

}  // namespace mrock::symbolic_operators::experimental::verify
