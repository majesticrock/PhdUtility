#pragma once
#include "StandardOperators.hpp"

namespace SymbolicOperators {
	class Continuum : public StandardOperators
	{
		virtual std::vector<Term> hamiltonian() const override;
		virtual std::vector<WickOperatorTemplate> templates() const override;
		virtual std::vector<std::vector<Term>> XP_basis() const override;
		virtual std::vector<std::vector<Term>> STD_basis() const override;
		virtual std::vector<std::unique_ptr<WickSymmetry>> symmetries() const override;

		virtual std::string get_subfolder() const override;
	};
}