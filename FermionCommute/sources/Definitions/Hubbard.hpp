#pragma once
#include "StandardOperators.hpp"

namespace SymbolicOperators {
	struct Hubbard : public StandardOperators {
		static const Momentum base_k_Q;

		static const Operator c_k_Q;
		static const Operator c_minus_k_Q;

		static const Operator c_k_Q_dagger;
		static const Operator c_minus_k_Q_dagger;

		// transversal magnon
		static const Operator c_k_Q_down_dagger;
		static const Operator c_k_Q_down;

		virtual std::vector<Term> hamiltonian() const override;
		virtual std::vector<WickOperatorTemplate> templates() const override;
		virtual std::vector<std::vector<Term>> XP_basis() const override;
		virtual std::vector<std::vector<Term>> STD_basis() const override;
		virtual std::vector<std::unique_ptr<WickSymmetry>> symmetries() const override;

		virtual std::string get_subfolder() const override;
	};
}