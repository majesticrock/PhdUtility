#pragma once
#include "WickTerm.hpp"
#include <vector>
#include <string>

namespace mrock::symbolic_operators {
	struct TermLoader {
		std::vector<WickTermCollector> M, N;
		void load(std::string const& folder, bool use_XP, int n_terms, int start_at = 0);
	};
}