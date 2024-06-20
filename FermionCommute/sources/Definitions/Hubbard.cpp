#include "Hubbard.hpp"

namespace SymbolicOperators {
	const Momentum Hubbard::base_k_Q{ momentum_pairs{ {1, 'k'} }, true };

	const Operator Hubbard::c_k_Q{ base_k_Q, SpinUp, false };
	const Operator Hubbard::c_minus_k_Q{ -base_k_Q, SpinDown, false };
	const Operator Hubbard::c_k_Q_dagger{ base_k_Q, SpinUp, true };
	const Operator Hubbard::c_minus_k_Q_dagger{ -base_k_Q, SpinDown, true };
	const Operator Hubbard::c_k_Q_down_dagger{ base_k_Q, SpinDown, true };
	const Operator Hubbard::c_k_Q_down{ base_k_Q, SpinDown, false };

	std::vector<Term> Hubbard::hamiltonian() const
	{
		const Term H_T(1, Coefficient("\\epsilon_0", Momentum('q')), SumContainer{ MomentumSum({ 'q' }), Sigma },
			std::vector<Operator>({
				Operator('q', 1, false, Sigma, true), Operator('q', 1, false, Sigma, false)
				}));

		const Term H_U(1, Coefficient("\\frac{U}{N}"), MomentumSum({ 'r', 'p', 'q' }), std::vector<Operator>({
			Operator('r', 1, false, SpinUp, true), Operator('p', 1, false, SpinDown, true),
			Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), SpinDown, false),
			Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), SpinUp, false),
			}));

		const Term H_V(1, Coefficient("\\tilde{V}", Momentum('q'), true),
			SumContainer{ MomentumSum({ 'r', 'p', 'q' }), IndexSum({ Sigma, SigmaPrime }) },
			std::vector<Operator>({
				Operator('r', 1, false, Sigma, true),
				Operator('p', 1, false, SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'p'), std::make_pair(-1, 'q') }), SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'r'), std::make_pair(1, 'q') }), Sigma, false),
				}));

		return { H_T, H_U, H_V };
	}
	std::vector<WickOperatorTemplate> Hubbard::templates() const
	{
		return {
			WickOperatorTemplate{ {IndexComparison{false, SpinDown, SpinUp}}, Momentum(), SC_Type, true },
			WickOperatorTemplate{ {IndexComparison{false, SpinDown, SpinUp}}, Momentum(momentum_pairs(), true), Eta_Type, true },
			WickOperatorTemplate{ {IndexComparison{true}}, Momentum(), Number_Type, false },
			WickOperatorTemplate{ {IndexComparison{true}}, Momentum(momentum_pairs(), true), CDW_Type, false }
		};
	}
	std::vector<std::vector<Term>> Hubbard::XP_basis() const
	{
		return {
			// 0: f + f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k, c_k })),
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
				}),
			// 1: eta + eta^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_Q, c_k })),
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
				}),
			// 2/3: g_up/down +
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q })),
				Term(1, std::vector<Operator>({ c_k_Q_dagger, c_k }))
				}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q })),
				Term(1, std::vector<Operator>({ c_minus_k_Q_dagger, c_minus_k }))
				}),
			// 4: transversal magnon, hermitian
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q_down })),
				Term(1, std::vector<Operator>({ c_k_Q_down_dagger, c_k }))
				}),
			// 5/6: n_up/down
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k }))
				}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
				}),
			// 7: f - f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k, c_k })),
				Term(-1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
				}),
			// 8: eta - eta^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_Q, c_k })),
				Term(-1, std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
				}),
			// 9/10: g_up/down -
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q })),
				Term(-1, std::vector<Operator>({ c_k_Q_dagger, c_k }))
				}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q })),
				Term(-1, std::vector<Operator>({ c_minus_k_Q_dagger, c_minus_k }))
				}),
			// 11: transversal magnon, antihermitian
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q_down })),
				Term(-1, std::vector<Operator>({ c_k_Q_down_dagger, c_k }))
				})
		};
	}
	std::vector<std::vector<Term>> Hubbard::STD_basis() const
	{
		return {
			// f, f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k, c_k }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_dagger }))
			}),
			// n_up/down
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k }))
			}),
			// g_up/down
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_dagger, c_minus_k_Q }))
			}),
			// eta, eta^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_minus_k_Q, c_k }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_minus_k_Q_dagger }))
			}),
			// transversal magnon
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_dagger, c_k_Q_down }))
			}),
			std::vector<Term>({
				Term(1, std::vector<Operator>({ c_k_Q_down_dagger, c_k }))
			})
		};
	}
	std::vector<std::unique_ptr<WickSymmetry>> Hubbard::symmetries() const
	{
		std::vector<std::unique_ptr<WickSymmetry>> ret;
		ret.reserve(3);
		ret.push_back(std::make_unique<SpinSymmetry>());
		ret.push_back(std::make_unique<TranslationalSymmetry>());
		ret.push_back(std::make_unique<PhaseSymmetry<SC_Type, CDW_Type>>());
		return ret;
	}
	std::string Hubbard::get_subfolder() const
	{
		return "hubbard/";
	}
}