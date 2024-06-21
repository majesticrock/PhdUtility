#pragma once

#include <SymbolicOperators/Momentum.hpp>
#include <SymbolicOperators/Operator.hpp>
#include <SymbolicOperators/Term.hpp>
#include <SymbolicOperators/WickOperatorTemplate.hpp>
#include <SymbolicOperators/WickSymmetry.hpp>

namespace SymbolicOperators {
	struct StandardOperators {
		static const Momentum base_k;
		static const Momentum base_x;

		static const Operator c_k; // c_{k up}
		static const Operator c_minus_k; // c_{-k down}

		static const Operator c_k_dagger; // c_{k up}^+
		static const Operator c_minus_k_dagger; // c_{-k down}^+

		virtual std::vector<Term> hamiltonian() const = 0;
		virtual std::vector<WickOperatorTemplate> templates() const = 0;
		virtual std::vector<std::vector<Term>> XP_basis() const = 0;
		virtual std::vector<std::vector<Term>> STD_basis() const = 0;
		virtual std::vector<std::unique_ptr<WickSymmetry>> symmetries() const = 0;

		virtual std::string get_subfolder() const = 0;
	};
}