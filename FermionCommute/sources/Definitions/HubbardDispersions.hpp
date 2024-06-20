#pragma once
#include "Hubbard.hpp"

namespace SymbolicOperators {
	struct HubbardDispersions : public Hubbard {
		virtual std::vector<std::vector<Term>> XP_basis() const override;
		virtual std::vector<std::vector<Term>> STD_basis() const override;

		virtual std::string get_subfolder() const override;
	};
}