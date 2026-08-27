#include <mrock/symbolic_operators/Commutation>

#include <Eigen/Sparse>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace mrock::symbolic_operators;
using SparseMatrix = Eigen::SparseMatrix<double>;

// Consider a chain of length L
constexpr int L = 10;
// Factor 2 for spins + 2 for occupied or unoccupied + one per lattice site
const int matrix_size = 2 * 2 * L;

const Operator c_K_sigma        = Operator(Momentum('K'), Index::Sigma,      false);
const Operator c_K_sigma_prime  = Operator(Momentum('K'), Index::SigmaPrime, false);

// Filled in main
std::array<double, L> cosines;

double wrap_momentum(int k) {
    while (k <  0) k += L;
    while (k >= L) k -= L; 
    return k;
}

SparseMatrix creation_operator(int lattice_site, int spin) {
    SparseMatrix ret(matrix_size, matrix_size);
    ret.insert(4 * lattice_site + 2 * spin + 1, 4 * lattice_site + 2 * spin) = 1.;
    return ret;
}
SparseMatrix annihilation_operator(int lattice_site, int spin) {
    SparseMatrix ret(matrix_size, matrix_size);
    ret.insert(4 * lattice_site + 2 * spin, 4 * lattice_site + 2 * spin + 1) = 1.;
    return ret;
}
SparseMatrix number_operator(int lattice_site, int spin) {
    SparseMatrix ret(matrix_size, matrix_size);
    ret.insert(4 * lattice_site + 2 * spin + 1, 4 * lattice_site + 2 * spin + 1) = 1.;
    return ret;
}

TermCollector get_symbolic_H() {
    const Term H_kin_symbolic(1, 
        Coefficient::RealInversionSymmetric("\\tilde{\\varepsilon}", MomentumList(Momentum('K'))),
        SumContainer{MomentumSum{'K'}, IndexSum{ Index::Sigma }},
        std::vector<Operator>({
            c_K_sigma.hermitian_conjugate(),
            c_K_sigma
        })
    );

    const Term H_int_symbolic(1, 
        Coefficient::RealInteraction("U", 
            MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}), 
            IndexWrapper{}
        ),
        SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
        std::vector<Operator>({
            c_K_sigma.hermitian_conjugate(),
            c_K_sigma_prime.with_momentum('P').hermitian_conjugate(),
            c_K_sigma_prime.with_momentum(Momentum("P-Q")),
            c_K_sigma.with_momentum(Momentum("K+Q"))
        })
    );

    return {H_kin_symbolic, H_int_symbolic};
}

TermCollector get_symbolic_eta() {
    const Term eta(2, 
        Coefficient::RealInteraction("\\alpha", 
            MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}), 
            IndexWrapper{}
        ),
        SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
        std::vector<Operator>({
            c_K_sigma.hermitian_conjugate(),
            c_K_sigma_prime.with_momentum('P').hermitian_conjugate(),
            c_K_sigma_prime.with_momentum(Momentum("P-Q")),
            c_K_sigma.with_momentum(Momentum("K+Q"))
        })
    );

    return {eta};
}

SparseMatrix get_matrix_H() {
    SparseMatrix H(matrix_size, matrix_size);
    for (int k=0; k<L; ++k) {
        H += cosines[k] * (number_operator(k, 0) + number_operator(k, 1));
    }
    for (int k=0; k<L; ++k) {
        for (int p=0; p<L; ++p) {
            for (int q=0; q<L; ++q) {
                H += cosines[wrap_momentum(k+q)] * cosines[wrap_momentum(p-q)] * (
                    creation_operator(k, 0) * creation_operator(p, 0) * annihilation_operator(wrap_momentum(p-q), 0) * annihilation_operator(wrap_momentum(k+q), 0)
                    + creation_operator(k, 0) * creation_operator(p, 1) * annihilation_operator(wrap_momentum(p-q), 1) * annihilation_operator(wrap_momentum(k+q), 0)
                    + creation_operator(k, 1) * creation_operator(p, 0) * annihilation_operator(wrap_momentum(p-q), 0) * annihilation_operator(wrap_momentum(k+q), 1)
                    + creation_operator(k, 1) * creation_operator(p, 1) * annihilation_operator(wrap_momentum(p-q), 1) * annihilation_operator(wrap_momentum(k+q), 1)
                );
            }
        }
    }
    return H;
}

SparseMatrix get_matrix_eta() {
    SparseMatrix eta(matrix_size, matrix_size);
    for (int k=0; k<L; ++k) {
        for (int p=0; p<L; ++p) {
            for (int q=0; q<L; ++q) {
                eta += cosines[wrap_momentum(k+q)] * cosines[wrap_momentum(p-q)]
                    * (cosines[wrap_momentum(k)] + cosines[wrap_momentum(p)] - cosines[wrap_momentum(k+q)] + cosines[wrap_momentum(p-q)])
                    * (
                        creation_operator(k, 0) * creation_operator(p, 0) * annihilation_operator(wrap_momentum(p-q), 0) * annihilation_operator(wrap_momentum(k+q), 0)
                        + creation_operator(k, 0) * creation_operator(p, 1) * annihilation_operator(wrap_momentum(p-q), 1) * annihilation_operator(wrap_momentum(k+q), 0)
                        + creation_operator(k, 1) * creation_operator(p, 0) * annihilation_operator(wrap_momentum(p-q), 0) * annihilation_operator(wrap_momentum(k+q), 1)
                        + creation_operator(k, 1) * creation_operator(p, 1) * annihilation_operator(wrap_momentum(p-q), 1) * annihilation_operator(wrap_momentum(k+q), 1)
                    );
            }
        }
    }
    return eta;
}

SparseMatrix symbolic_to_matrix(const TermCollector& terms) {
    SparseMatrix result(matrix_size, matrix_size);

    for (const auto& term : terms) {
        std::array<int, 256> momentum_values{};
        std::array<int, 256> index_values{};

        const auto evaluate_momentum = [&momentum_values](const Momentum& momentum) {
            int value = 0;
            for (const auto& symbol : momentum.momentum_list) {
                value += symbol.factor * momentum_values[static_cast<unsigned char>(symbol.name)];
            }
            if (momentum.add_PI) {
                value += L / 2;
            }
            return wrap_momentum(value);
        };

        const auto evaluate_coefficient = [&evaluate_momentum](const Coefficient& coefficient) {
            double value = 1.;
            if (coefficient.name == "\\tilde{\\varepsilon}") {
                value *= cosines[evaluate_momentum(coefficient.momenta.front())];
            } 
            else if (coefficient.name == "U") {
                value *= cosines[evaluate_momentum(coefficient.momenta[0] + coefficient.momenta[2])]
                        * cosines[evaluate_momentum(coefficient.momenta[1] - coefficient.momenta[2])];
            } 
            else if (coefficient.name == "\\alpha") {
                value *= cosines[evaluate_momentum(coefficient.momenta[0] + coefficient.momenta[2])]
                        * cosines[evaluate_momentum(coefficient.momenta[1] - coefficient.momenta[2])]
                        * ( 
                            cosines[evaluate_momentum(coefficient.momenta[0])] + cosines[evaluate_momentum(coefficient.momenta[1])]
                            - cosines[evaluate_momentum(coefficient.momenta[0] + coefficient.momenta[2])]
                            - cosines[evaluate_momentum(coefficient.momenta[1] - coefficient.momenta[2])]
                        );
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
                if (wrap_momentum(evaluate_momentum(delta.first)) !=
                    wrap_momentum(evaluate_momentum(delta.second))) {
                    return;
                }
            }
            for (const auto& delta : term.delta_indizes) {
                const auto index_value = [&index_values](const Index index) {
                    return is_mutable(index) ? index_values[static_cast<unsigned char>(index)]
                                             : (index == Index::SpinDown ? 1 : 0);
                };
                if (index_value(delta.first) != index_value(delta.second)) {
                    return;
                }
            }

            SparseMatrix term_matrix(matrix_size, matrix_size);
            term_matrix.setIdentity();
            for (const auto& op : term.operators) {
                const int momentum = wrap_momentum(evaluate_momentum(op.momentum));
                const int spin = is_mutable(op.first_index())
                                     ? index_values[static_cast<unsigned char>(op.first_index())]
                                     : (op.first_index() == Index::SpinDown ? 1 : 0);
                const SparseMatrix operator_matrix = op.is_daggered ? creation_operator(momentum, spin)
                                                                       : annihilation_operator(momentum, spin);
                term_matrix = term_matrix * operator_matrix;
            }
            result += coefficient * term_matrix;
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
                momentum_values[name] = value;
                assign_momenta(sum_index + 1U);
            }
        };
        assign_momenta(0U);
    }
    return result;
}

int main() {
    const SparseMatrix creation = creation_operator(L/2, 0);
    const SparseMatrix annihilation_transpose = annihilation_operator(L/2, 0).transpose();
    if ((creation - annihilation_transpose).norm() != 0.) {
        std::cerr << "c^+ != c" << std::endl;
        return 1;
    }
    if ((creation_operator(L/2, 0) * annihilation_operator(L/2,0) - number_operator(L/2, 0)).norm() != 0.) {
        std::cerr << "n != c^+ c" << std::endl;
        return 1;
    }

    for (int i=0; i < L; ++i) {
        cosines[i] = std::cos(std::numbers::pi * (2*i - L));
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