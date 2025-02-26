/**
 * @file SymbolicSum.hpp
 * @brief Defines the SymbolicSum template struct for symbolic summation operations.
 * 
 * This file contains the definition of the SymbolicSum template struct, which is used
 * to represent and manipulate symbolic summation operations. It provides various
 * constructors, serialization support, and comparison operators.
 */

#pragma once
#include <mrock/utility/VectorWrapper.hpp>
#include <ostream>

namespace mrock::symbolic_operators {

	/**
	 * @brief A struct representing a symbolic summation operation.
	 * 
	 * @tparam SumIndex The type of the summation index.
	 */
	template<class SumIndex>
	struct SymbolicSum {
		std::vector<SumIndex> summations; ///< The vector of summation indices.

		/**
		 * @brief Serializes the SymbolicSum object.
		 * 
		 * @tparam Archive The type of the archive.
		 * @param ar The archive to serialize to.
		 * @param version The version of the serialization format.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->summations;
		}

		/**
		 * @brief Default constructor.
		 */
		SymbolicSum() = default;

		/**
		 * @brief Constructs a SymbolicSum with a single summation index.
		 * 
		 * @param sum_index The summation index.
		 */
		SymbolicSum(SumIndex sum_index) 
			: summations(1U, sum_index) {}

		/**
		 * @brief Constructs a SymbolicSum with a vector of summation indices.
		 * 
		 * @param _indizes The vector of summation indices.
		 */
		SymbolicSum(const std::vector<SumIndex>& _indizes)
			: summations(_indizes) {}

		/**
		 * @brief Constructs a SymbolicSum with a moved vector of summation indices.
		 * 
		 * @param _indizes The vector of summation indices to move.
		 */
		SymbolicSum(std::vector<SumIndex>&& _indizes)
			: summations(std::move(_indizes)) {}

		/**
		 * @brief Constructs a SymbolicSum with an initializer list of summation indices.
		 * 
		 * @param init The initializer list of summation indices.
		 */
		SymbolicSum(std::initializer_list<SumIndex> init)
			 : summations(std::move(init)) {}

		/**
		 * @brief Checks if a given index is part of the summation indices.
		 * 
		 * @param what The index to check.
		 * @return True if the index is part of the summation indices, false otherwise.
		 */
		inline bool is_summed_over(SumIndex what) const {
			for (const auto& _s : this->summations) {
				if (_s == what) return true;
			}
			return false;
		}

		MROCK_VECTOR_WRAPPER_FILL_MEMBERS(SumIndex, summations);

		/**
		 * @brief Compares two SymbolicSum objects.
		 * 
		 * @param rhs The other SymbolicSum to compare with.
		 * @return The result of the comparison.
		 */
		inline auto operator<=>(const SymbolicSum<SumIndex>& rhs) const = default;
	};

	/**
	 * @brief Outputs the SymbolicSum object to an output stream.
	 * 
	 * @tparam SumIndex The type of the summation index.
	 * @param os The output stream.
	 * @param sum The SymbolicSum object to output.
	 * @return The output stream.
	 */
	template<class SumIndex>
	std::ostream& operator<<(std::ostream& os, SymbolicSum<SumIndex> const& sum) {
		if (sum.empty()) return os;
		os << "\\sum_{ ";
		for (const auto& index : sum) {
			os << index << " ";
		}
		os << "} ";
		return os;
	}
} // namespace mrock::symbolic_operators