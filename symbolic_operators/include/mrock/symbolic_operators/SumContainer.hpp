/**
 * @file SumContainer.hpp
 * @brief Defines the SumContainer structure and related operators for symbolic operations.
 */

#pragma once

#include "SymbolicSum.hpp"
#include <mrock/utility/RangeUtility.hpp>
#include "IndexWrapper.hpp"
#include "MomentumSymbol.hpp"

namespace mrock::symbolic_operators {
	/**
	 * @typedef IndexSum
	 * @brief Typedef for SymbolicSum with Index type.
	 */
	typedef SymbolicSum<Index> IndexSum;

	/**
	 * @typedef MomentumSum
	 * @brief Typedef for SymbolicSum with MomentumSymbol::name_type type.
	 */
	typedef SymbolicSum<MomentumSymbol::name_type> MomentumSum;

	/**
	 * @struct SumContainer
	 * @brief A container for holding symbolic sums of momenta and spins.
	 * 
	 * Sums are contained within the \c SumContainer class. 
	 * It hosts both sums of momenta and sums of spins, each one is accessible via the appropriate class member
	 * and its \c operator[], e.g., \c container.momenta[i].
	 * 
	 * @sa Index, MomentumSymbol, MomentumSymbol::name_type
	 */
	struct SumContainer {
		MomentumSum momenta; ///< Container for momentum sums.
		IndexSum spins; ///< Container for spin sums.

		/**
		 * @brief Serializes the SumContainer object.
		 * @tparam Archive The type of the archive.
		 * @param ar The archive to serialize to.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->momenta;
			ar& this->spins;
		}

		/**
		 * @brief Appends another SumContainer to this one.
		 * @param other The other SumContainer to append.
		 * @return Reference to this SumContainer.
		 */
		SumContainer& append(const SumContainer& other);

		/**
		 * @brief Appends a MomentumSum to this SumContainer.
		 * @param other The MomentumSum to append.
		 * @return Reference to this SumContainer.
		 */
		SumContainer& append(const MomentumSum& other);

		/**
		 * @brief Appends an IndexSum to this SumContainer.
		 * @param other The IndexSum to append.
		 * @return Reference to this SumContainer.
		 */
		SumContainer& append(const IndexSum& other);

		/**
		 * @brief Pushes back a momentum into the momenta container.
		 * @param momentum The momentum to push back.
		 */
		inline void push_back(const MomentumSymbol::name_type momentum);

		/**
		 * @brief Pushes back a spin into the spins container.
		 * @param spin The spin to push back.
		 */
		inline void push_back(const Index spin);

		/**
		 * @brief Checks if the container has any momenta.
		 * @return True if the container has momenta, false otherwise.
		 */
		inline bool has_momentum() const noexcept;

		/**
		 * @brief Checks if the container has any spins.
		 * @return True if the container has spins, false otherwise.
		 */
		inline bool has_spins() const noexcept;
	};

	/**
	 * @brief Equality operator for SumContainer.
	 * @param lhs The left-hand side SumContainer.
	 * @param rhs The right-hand side SumContainer.
	 * @return True if both SumContainers are equal, false otherwise.
	 */
	inline bool operator==(const SumContainer& lhs, const SumContainer& rhs) {
		return (lhs.momenta == rhs.momenta && lhs.spins == rhs.spins);
	}

	/**
	 * @brief Inequality operator for SumContainer.
	 * @param lhs The left-hand side SumContainer.
	 * @param rhs The right-hand side SumContainer.
	 * @return True if both SumContainers are not equal, false otherwise.
	 */
	inline bool operator!=(const SumContainer& lhs, const SumContainer& rhs) {
		return !(lhs == rhs);
	}

	/**
	 * @brief Stream insertion operator for SumContainer.
	 * @param os The output stream.
	 * @param sums The SumContainer to insert into the stream.
	 * @return Reference to the output stream.
	 */
	inline std::ostream& operator<<(std::ostream& os, const SumContainer& sums) {
		os << sums.momenta << sums.spins;
		return os;
	}

	// Inline definitions
	void SumContainer::push_back(const MomentumSymbol::name_type momentum) {
		this->momenta.push_back(momentum);
	}
	void SumContainer::push_back(const Index spin) {
		this->spins.push_back(spin);
	}
	bool SumContainer::has_momentum() const noexcept {
		return !momenta.empty();
	}
	bool SumContainer::has_spins() const noexcept {
		return !spins.empty();
	}
} // namespace mrock::symbolic_operators