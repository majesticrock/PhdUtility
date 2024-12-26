#include "compare_test.hpp"

using namespace mrock::SymbolicOperators;

int main(int argc, char** argv) {
    const std::string COMPARE_DIR = "../../tests/correct_continuum/";
    // Setup
    const Momentum base_k = Momentum('k');

    const Operator c_k = Operator{ base_k, Index::SpinUp, false };
    const Operator c_minus_k = Operator{ -base_k, Index::SpinDown, false };
    const Operator c_k_dagger = Operator{ base_k, Index::SpinUp, true };
    const Operator c_minus_k_dagger = Operator{ -base_k, Index::SpinDown, true };

    const Term H_Kin(1, Coefficient("\\epsilon_0", Momentum('q')), 
        SumContainer{ MomentumSum({ 'q' }), Index::Sigma },
	    std::vector<Operator>({
		    Operator('q', 1, false, Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)
		}));

    const Term H_Ph(-1, Coefficient("g", MomentumList({ 'q', 'p' })), 
        SumContainer{ MomentumSum({'q', 'p'}) }, 
        std::vector<Operator>({
    	    c_k_dagger.with_momentum('q'), c_minus_k_dagger.with_momentum('q'),
    	    c_minus_k.with_momentum('p'), c_k.with_momentum('p')
    	}));

    const Term H_C(IntFractional(1, 2), Coefficient("V", Momentum('q')),
    	SumContainer{ MomentumSum({ 'r', 'p', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
    	std::vector<Operator>({
    		Operator('r', 1, false, Index::Sigma, true),
    		Operator('p', 1, false, Index::SigmaPrime, true),
    		Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
    		Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), Index::Sigma, false),
    	}));

    const Term H_BG(-1, Coefficient("\\rho"), 
        SumContainer{ MomentumSum({ 'q' }), Index::Sigma },
    	std::vector<Operator>({
    		Operator(Momentum("q"), Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)
    	}));

    // Hamiltonian of the continuum  system
    const std::vector<Term> H({ H_Kin, H_Ph, H_C, H_BG });

    // f + f^+
    const std::vector<Term> base_term({
	    Term(1, std::vector<Operator>({ c_minus_k, c_k })),
	    Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
	});

    // Wick templates
    const std::vector<WickOperatorTemplate> templates({
			WickOperatorTemplate{ {IndexComparison{false, Index::SpinDown, Index::SpinUp}}, Momentum(), SC_Type, true },
			WickOperatorTemplate{ {IndexComparison{true}}, Momentum(), Number_Type, false }
		});
    
    // Applicable symmetries
    std::vector<std::unique_ptr<WickSymmetry>> symmetries;
    symmetries.push_back(std::make_unique<SpinSymmetry>());
    symmetries.push_back(std::make_unique<TranslationalSymmetry>());
    symmetries.push_back(std::make_unique<PhaseSymmetry<SC_Type>>());

    sym_op_test::SymOpTest tester(COMPARE_DIR);
    if (std::filesystem::exists(COMPARE_DIR)) {
        std::cout << "Found compare dir " << COMPARE_DIR << std::endl;
    }
    else {
        std::cout << "Did not find compare dir " << COMPARE_DIR << "\nCreating comparison data now..." << std::endl;
    }
    return tester.perform_test(H, base_term, templates, symmetries, (argc > 1) || !(std::filesystem::exists(COMPARE_DIR)));
}