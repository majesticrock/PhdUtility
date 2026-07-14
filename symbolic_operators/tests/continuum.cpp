/**
 * @file continuum.cpp
 * @brief Example code for defining and using the continuum-model operators from doi:10.21468/SciPostPhys.19.3.067.
 *
 * This file demonstrates how to construct a Hamiltonian for an interacting electron gas
 * in the symbolic-operator framework. It defines the kinetic, phonon-mediated pairing,
 * Coulomb-interaction, and background-density terms, then evaluates the commutator of
 * the Hamiltonian with a pair-annihilation operator and applies Wick's theorem.
 * The resulting expressions are compared against serialized reference data stored in the
 * comparison directory.
 *
 * These reference data should be generated before changes are made to the library.
 * Thereby, one can validate that any changes do not break existing results.
 */

#include "compare_test.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/ExpectationValues>

#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace mrock::symbolic_operators;

int main(int argc, char** argv) {
    const std::string COMPARE_DIR = "../../../symbolic_operators/tests/correct_continuum/";
    // Setup a momentum object representing k. This is later reused to define the electron operators.
    const Momentum base_k = Momentum('k');

    // Setup a few operators to be used later.
    // The constructor takes the momentum, the spin index, and a boolean flag indicating whether
    // the operator is a creation (true) or annihilation (false) operator.
    const Operator c_k = Operator{base_k, Index::SpinUp, false};
    const Operator c_minus_k = Operator{-base_k, Index::SpinDown, false};
    const Operator c_k_dagger = Operator{base_k, Index::SpinUp, true};
    const Operator c_minus_k_dagger = Operator{-base_k, Index::SpinDown, true};

    // Setup the kinetic term: sum_{q,sigma} epsilon_0(q) c^\dagger_{q,sigma} c_{q,sigma}.
    // A Term is constructed from a prefactor, a coefficient, the sums over momenta and indices,
    // and the ordered list of operators appearing in the term.
    const Term H_Kin(1, Coefficient("\\epsilon_0", Momentum('q')), SumContainer{MomentumSum({'q'}), Index::Sigma},
                     std::vector<Operator>(
                         {Operator('q', 1, false, Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)}));

    // Setup the phonon-mediated pairing term: -sum_{q,p} g(q,p) c^\dagger_{q,up} c^\dagger_{-q,down} c_{-p,down}
    // c_{p,up}. The operators are created from the previously defined base operators and then assigned new momenta.
    const Term H_Ph(-1, Coefficient("g", MomentumList({'q', 'p'})), SumContainer{MomentumSum({'q', 'p'})},
                    std::vector<Operator>({c_k_dagger.with_momentum('q'), c_minus_k_dagger.with_momentum('q'),
                                           c_minus_k.with_momentum('p'), c_k.with_momentum('p')}));

    // Setup the Coulomb interaction term: 1/2 sum_{r,p,q} V(q) c^\dagger_{r,sigma} c^\dagger_{p,sigma'}
    // c_{p+q,sigma'} c_{r-q,sigma}. The momentum shifts are encoded directly in the operator constructors.
    const Term H_C(IntFractional(1, 2), Coefficient("V", Momentum('q')),
                   SumContainer{MomentumSum({'r', 'p', 'q'}), IndexSum({Index::Sigma, Index::SigmaPrime})},
                   std::vector<Operator>({
                       Operator('r', 1, false, Index::Sigma, true),
                       Operator('p', 1, false, Index::SigmaPrime, true),
                       Operator(std::vector<MomentumSymbol>({MomentumSymbol(1, 'p'), MomentumSymbol(-1, 'q')}),
                                Index::SigmaPrime, false),
                       Operator(std::vector<MomentumSymbol>({MomentumSymbol(1, 'r'), MomentumSymbol(1, 'q')}),
                                Index::Sigma, false),
                   }));

    // Setup the background-density term: -sum_q rho c^\dagger_{q,sigma} c_{q,sigma}.
    // Here the coefficient is a constant and the operator is created from a momentum expression.
    const Term H_BG(-1, Coefficient("\\rho"), SumContainer{MomentumSum({'q'}), Index::Sigma},
                    std::vector<Operator>(
                        {Operator(Momentum("q"), Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)}));

    // Hamiltonian of the continuum  system
    const std::vector<Term> H({H_Kin, H_Ph, H_C, H_BG});

    // Define the operator strings that will be used as the commutation targets.
    // Each Term is just a product of operators here, without any coefficient or sums.
    const std::vector<Term> base_term({Term(1, std::vector<Operator>({c_minus_k, c_k})),
                                       Term(1, std::vector<Operator>({c_k_dagger, c_minus_k_dagger}))});

    // Wick templates
    const std::vector<WickOperatorTemplate> templates(
        {WickOperatorTemplate{{SC_Comparison}, Momentum(), OperatorType::SC},
         WickOperatorTemplate{{Num_Comparison}, Momentum(), OperatorType::Number}});

    // Applicable symmetries
    std::vector<std::unique_ptr<WickSymmetry>> symmetries;
    symmetries.push_back(std::make_unique<SpinSymmetry>());
    symmetries.push_back(std::make_unique<InversionSymmetry>());
    symmetries.push_back(std::make_unique<PhaseSymmetry<OperatorType::SC>>());

    sym_op_test::SymOpTest tester(COMPARE_DIR);
    if (std::filesystem::exists(COMPARE_DIR)) {
        std::cout << "Found compare dir " << std::filesystem::canonical(COMPARE_DIR) << std::endl;
    } else {
        std::cout << "Did not find compare dir " << COMPARE_DIR << "\nCreating comparison data now..." << std::endl;
    }
    return tester.perform_test(H, base_term, templates, symmetries,
                               (argc > 1) || !(std::filesystem::exists(COMPARE_DIR)));
}