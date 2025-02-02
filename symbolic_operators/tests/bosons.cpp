/**
 * @file bosons.cpp
 * @brief Example code for defining and using bosonic operators.
 *
 * This file demonstrates how to define a Hamiltonian using bosonic operators
 * and how to perform commutation operations on them.
 */

#include "compare_test.hpp"
#include <string>
#include <iostream>
#include <filesystem>

using namespace mrock::symbolic_operators;

const std::string begin_align = "\\begin{align*}\n\t";
const std::string end_align = "\\end{align*}\n";
const std::string file_names[2] = { "first_commutation.bin", "second_commutation.bin" };
const std::string COMPARE_DIR = "../../../symbolic_operators/tests/correct_bosons/";

int main(int argc, char** argv) {
    bool is_baseline = (argc > 1) || !(std::filesystem::exists(COMPARE_DIR));

    // Define the Hamiltonian
    const Term hopping(1, Coefficient::HoneyComb("\\gamma", Momentum('k'), false, false),
        MomentumSum{'k'},
        std::vector<Operator>({
            Operator::Boson(Momentum('k'), Index::TypeA, true),
            Operator::Boson(Momentum('k'), Index::TypeB, false)
        }));

    const Term bogo(IntFractional(1, 2), Coefficient::HoneyComb("\\gamma'", Momentum('k'), false, false),
        MomentumSum{'k'},
        std::vector<Operator>({
            Operator::Boson(Momentum('k'), Index::TypeA, true),
            Operator::Boson(Momentum('k', -1), Index::TypeB, true)
        }));

    const Term chemical_potential(-1, Coefficient::Constant("\\mu", Index::Sigma), 
        SumContainer{ MomentumSum{'k'}, IndexSum{Index::Sigma} },
        std::vector<Operator>({
            Operator::Boson(Momentum('k'), Index::Sigma, true),
            Operator::Boson(Momentum('k'), Index::Sigma, false)
        }));

    const std::vector<Term> hamiltonian { hopping, hopping.hermitian_conjugate(), bogo, bogo.hermitian_conjugate(), chemical_potential };
    
    std::cout << begin_align << "H =" << hamiltonian << end_align << std::endl;

    // Define the commutation targets
    // It is important to choose a different name for the momentum here than in the Hamiltonian
    const Term to_commute_1(1, std::vector<Operator>({
            Operator::Boson(Momentum('l'), Index::TypeA, true),
            Operator::Boson(Momentum('l', -1), Index::TypeB, true)
        }));

    const std::vector<Term> to_commute_2{ 
        Term(1, SumContainer{ MomentumSum{'q'} },
            std::vector<Operator>({
                Operator::Boson(Momentum("l+q"), Index::TypeA, true),
                Operator::Boson(Momentum("l"), Index::TypeA, false)
            })),
        Term(1, SumContainer{ MomentumSum{'q'} },
            std::vector<Operator>({
                Operator::Boson(Momentum("l-q"), Index::TypeB, true),
                Operator::Boson(Momentum("l"), Index::TypeB, false)
            }))
        };

    // Compute the commutators
    std::vector<Term> result_1 = commutator(hamiltonian, to_commute_1);
    clean_up(result_1);
    std::vector<Term> result_2 = commutator(hamiltonian, to_commute_2);
    clean_up(result_2);

    std::cout << begin_align << "[H, " << to_commute_1.to_string_without_prefactor() << "] = " << result_1 << end_align << std::endl;
    std::cout << begin_align << "[H, " << to_string_without_prefactor(to_commute_2) << "]  = " << result_2 << end_align << std::endl;

    // Define tester
    sym_op_test::SymOpTest tester(COMPARE_DIR);

    if (is_baseline) tester.save_as_comparison(file_names[0], result_1);
    if(!tester.load_and_test(file_names[0], result_1)) return 1;

    if (is_baseline) tester.save_as_comparison(file_names[1], result_2);
    if(!tester.load_and_test(file_names[1], result_2)) return 1;

    return 0;
}