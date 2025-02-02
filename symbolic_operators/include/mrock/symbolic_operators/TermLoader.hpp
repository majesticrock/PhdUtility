/**
 * @file TermLoader.hpp
 * @brief Header file for the TermLoader structure in the symbolic_operators namespace.
 */

#pragma once
#include "WickTerm.hpp"
#include <vector>
#include <string>

namespace mrock::symbolic_operators {

	/**
	 * @struct TermLoader
	 * @brief A structure to load and manage Wick terms.
	 */
	struct TermLoader {
		std::vector<WickTermCollector> M; ///< Vector to store Wick terms of the dynamical matrix M
		std::vector<WickTermCollector> N; ///< Vector to store Wick terms of the norm matrix N

		/**
		 * @brief Loads Wick terms from a specified folder.
		 * 
		 * @param folder The path to the folder containing the terms.
		 * @param use_XP A boolean flag to indicate whether to use the XP basis.
		 * @param n_terms The number of terms to load.
		 * @param start_at The starting index for loading terms (default is 0).
		 */
		void load(std::string const& folder, bool use_XP, int n_terms, int start_at = 0);
	};

} // namespace mrock::symbolic_operators