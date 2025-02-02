/**
 * @file WickOperator.hpp
 * @brief Defines the WickOperator structure used in symbolic operators.
 */

#pragma once
#include <iostream>
#include "Momentum.hpp"
#include "IndexWrapper.hpp"
#include "OperatorType.hpp"

namespace mrock::symbolic_operators {

	/**
	 * @class WickOperator
	 * @brief A structure representing a Wick operator.
	 */
	struct WickOperator {
		OperatorType type{ OperatorType::Undefined_Type }; ///< The type of the operator.
		bool is_daggered{}; ///< Indicates if the operator is daggered.
		Momentum momentum; ///< The momentum associated with the operator.
		IndexWrapper indizes; ///< The indices associated with the operator.

		/**
		 * @brief Serializes the WickOperator object.
		 * 
		 * @tparam Archive The type of the archive.
		 * @param ar The archive object.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& type;
			ar& is_daggered;
			ar& momentum;
			ar& indizes;
		}

		/**
		 * @brief Constructs a WickOperator object.
		 * 
		 * @param _type The type of the operator.
		 * @param _is_daggered Whether the operator is daggered.
		 * @param _momentum The momentum of the operator.
		 * @param _indizes The indices of the operator.
		 */
		WickOperator(const OperatorType& _type, const bool _is_daggered, const Momentum& _momentum, const IndexWrapper& _indizes = IndexWrapper());

		/**
		 * @brief Constructs a WickOperator object.
		 * 
		 * @param _type The type of the operator.
		 * @param _is_daggered Whether the operator is daggered.
		 * @param _momentum The momentum of the operator.
		 * @param _index The index of the operator.
		 */
		WickOperator(const OperatorType& _type, const bool _is_daggered, const Momentum& _momentum, const Index _index);

		/**
		 * @brief Default constructor for WickOperator.
		 */
		WickOperator() = default;

		/**
		 * @brief Constructs a WickOperator object from a string expression.
		 * 
		 * @param expression The string expression.
		 */
		WickOperator(const std::string& expression);

		/**
		 * @brief Checks if the operator uses a specific index.
		 * 
		 * @param index The index to check.
		 * @return true if the operator uses the index.
		 * @return false otherwise.
		 */
		inline bool uses_index(const Index index) const noexcept;

		/**
		 * @brief Checks if the operator depends on a specific momentum.
		 * 
		 * @param momentum The momentum to check.
		 * @return true if the operator depends on the momentum.
		 * @return false otherwise.
		 */
		inline bool depends_on(const MomentumSymbol::name_type momentum) const noexcept;

		/**
		 * @brief Removes a momentum contribution from the operator.
		 * 
		 * @param value The momentum value to remove.
		 */
		inline void remove_momentum_contribution(const MomentumSymbol::name_type value);
	};

	/**
	 * @brief Stream insertion operator for WickOperator.
	 * 
	 * @param os The output stream.
	 * @param op The WickOperator object.
	 * @return std::ostream& The updated output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const WickOperator& op);

	/**
	 * @brief Stream insertion operator for a vector of WickOperator objects.
	 * 
	 * @param os The output stream.
	 * @param ops The vector of WickOperator objects.
	 * @return std::ostream& The updated output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops);

	// Inline definitions
	bool WickOperator::uses_index(const Index index) const noexcept {
		for (const auto& idx : this->indizes) {
			if (idx == index) return true;
		}
		return false;
	}
	inline bool WickOperator::depends_on(const MomentumSymbol::name_type momentum) const noexcept {
		return this->momentum.is_used_at(momentum) != -1;
	}
	inline void WickOperator::remove_momentum_contribution(const MomentumSymbol::name_type value) {
		momentum.remove_contribution(value);
	}
} // namespace mrock::symbolic_operators