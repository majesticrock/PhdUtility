#pragma once
#include "WickTerm.hpp"
#include <vector>
#include <string>

namespace SymbolicOperators {
	struct TermLoader {
		std::vector<WickTermCollector> M, N;
		void load(std::string const& folder, bool use_XP, int n_terms, int start_at = 0);
	};
}