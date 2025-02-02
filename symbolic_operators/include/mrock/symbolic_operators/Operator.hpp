/**
 * @file Operator.hpp
 * @brief Defines the Operator struct and related functions for symbolic operators.
 */

#pragma once

#include "Momentum.hpp"
#include "IndexWrapper.hpp"

namespace mrock::symbolic_operators {

	/**
	 * @struct Operator
	 * @brief Represents a symbolic operator with momentum, indices, and properties.
	 * 
	 * This class represents the standard fermionic or bosonic creation and annihilation operators. 
	 * You can specify its momentum, its indizes and whether it is supposed to be daggered (a creation operator) 
	 * or not (an annihilation operator).
	 * 
	 * @sa Momentum, IndexWrapper
	 */
	struct Operator {
		Momentum momentum; ///< The momentum associated with the operator.
		IndexWrapper indizes; ///< Contains all indices, standard: first index = spin, all others arbitrary, e.g., orbitals, bands etc.
		bool is_daggered{}; ///< Indicates if the operator is daggered (conjugate transpose).
		bool is_fermion{ true }; ///< Indicates if the operator is a fermion. This of course impacts the commutation relation [O', O^+]_{+/-} = delta_{O,O'}, where the plus applies to fermions and the minus to bosons.

		/**
		 * @brief Serializes the Operator object.
		 * @tparam Archive The type of the archive.
		 * @param ar The archive to serialize to.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& momentum;
			ar& indizes;
			ar& is_daggered;
			ar& is_fermion;
		}

		/**
		 * @brief Default constructor.
		 */
		Operator() = default;

		/**
		 * @brief Constructs an Operator with specified momentum, indices, daggered state, and fermion state.
		 * @param _momentum The momentum of the operator.
		 * @param _indizes The indices of the operator.
		 * @param _is_daggered The daggered state of the operator.
		 * @param _is_fermion The fermion state of the operator (default is true).
		 */
		Operator(const Momentum& _momentum, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);

		/**
		 * @brief Constructs an Operator with specified momentum symbols, indices, daggered state, and fermion state.
		 * @param _momentum The momentum symbols of the operator.
		 * @param _indizes The indices of the operator.
		 * @param _is_daggered The daggered state of the operator.
		 * @param _is_fermion The fermion state of the operator (default is true).
		 */
		Operator(const momentum_symbols& _momentum, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);

		/**
		 * @brief Constructs an Operator with specified momentum symbol name, addition flag, indices, daggered state, and fermion state.
		 * @param _momentum The name of the momentum symbol.
		 * @param add_Q Flag to indicate if Q should be added. Q has the property 2Q = 0, e.g., (pi,pi) on a unit square lattice.
		 * @param _indizes The indices of the operator.
		 * @param _is_daggered The daggered state of the operator.
		 * @param _is_fermion The fermion state of the operator (default is true).
		 */
		Operator(const MomentumSymbol::name_type _momentum, bool add_Q, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);

		/**
		 * @brief Constructs an Operator with specified momentum symbol name, sign, addition flag, indices, daggered state, and fermion state.
		 * @param _momentum The name of the momentum symbol.
		 * @param sign The sign of the momentum.
		 * @param add_Q Flag to indicate if Q should be added. Q has the property 2Q = 0, e.g., (pi,pi) on a unit square lattice.
		 * @param _indizes The indices of the operator.
		 * @param _is_daggered The daggered state of the operator.
		 * @param _is_fermion The fermion state of the operator (default is true).
		 */
		Operator(const MomentumSymbol::name_type _momentum, int sign, bool add_Q, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);

		/**
		 * @brief Creates a Boson operator with specified momentum and indices.
		 * @param _momentum The momentum of the operator.
		 * @param _indizes The indices of the operator.
		 * @param _is_daggered The daggered state of the operator.
		 * @return A Boson operator.
		 */
		inline static Operator Boson(const Momentum& _momentum, const IndexWrapper _indizes, bool _is_daggered) {
			return Operator(_momentum, _indizes, _is_daggered, false);
		}

		/**
		 * @brief Creates a Boson operator with specified momentum.
		 * @param _momentum The momentum of the operator.
		 * @param _is_daggered The daggered state of the operator.
		 * @return A Boson operator.
		 */
		inline static Operator Boson(const Momentum& _momentum, bool _is_daggered) {
			return Operator(_momentum, IndexWrapper{}, _is_daggered, false);
		}

		/**
		 * @brief Toggles the daggered state of the operator.
		 * @return A reference to *this
		 */
		inline Operator& hermitian_conjugate_inplace();

		/**
		 * @brief Creates hermitian conjugate of this as a new object.
		 * @return Returns the new object.
		 */
		inline Operator hermitian_conjugate() const;

		/**
		 * @brief Creates a new operator with updated momentum.
		 * @param new_momentum The new momentum to set.
		 * @return A new operator with the updated momentum.
		 */
		inline Operator with_momentum(Momentum const& new_momentum) const;

		/**
		 * @brief Creates a new operator with updated momentum symbol name.
		 * @param new_momentum The new momentum symbol name to set.
		 * @return A new operator with the updated momentum symbol name.
		 */
		inline Operator with_momentum(const MomentumSymbol::name_type new_momentum) const;

		/**
		 * @brief Creates a new operator by adding momentum.
		 * @param to_add The momentum to add.
		 * @return A new operator with the added momentum.
		 */
		inline Operator add_momentum(Momentum const& to_add) const;

		/**
		 * @brief Creates a new operator by adding momentum symbol name.
		 * @param to_add The momentum symbol name to add.
		 * @return A new operator with the added momentum symbol name.
		 */
		inline Operator add_momentum(const MomentumSymbol::name_type to_add) const;

		/**
		 * @brief Removes a momentum contribution from the operator.
		 * @param value The momentum symbol name to remove.
		 */
		inline void remove_momentum_contribution(const MomentumSymbol::name_type value);

		/**
		 * @brief Returns the first index of the operator.
		 * @return The first index if the operator has indices, otherwise Index::NoIndex.
		 */
		inline Index first_index() const;

		/**
		 * @brief Sets the first index of the operator.
		 * @param index The index to set as the first index.
		 */
		inline void set_first_index(Index index);
	};

	/**
	 * @brief Equality operator for Operator.
	 * @param lhs The left-hand side operator.
	 * @param rhs The right-hand side operator.
	 * @return True if the operators are equal, false otherwise.
	 */
	inline bool operator==(const Operator& lhs, const Operator& rhs) {
		if (lhs.is_fermion != rhs.is_fermion) return false;
		if (lhs.is_daggered != rhs.is_daggered) return false;
		if (lhs.indizes != rhs.indizes) return false;
		return (lhs.momentum == rhs.momentum);
	}

	/**
	 * @brief Inequality operator for Operator.
	 * @param lhs The left-hand side operator.
	 * @param rhs The right-hand side operator.
	 * @return True if the operators are not equal, false otherwise.
	 */
	inline bool operator!=(const Operator& lhs, const Operator& rhs) {
		return !(lhs == rhs);
	}

	/**
	 * @brief Stream insertion operator for Operator.
	 * @param os The output stream.
	 * @param op The operator to insert into the stream.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const Operator& op);

	/**
	 * @brief Stream insertion operator for a vector of Operators.
	 * @param os The output stream.
	 * @param ops The vector of operators to insert into the stream.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const std::vector<Operator>& ops);


	// Inline definitions
	Operator& Operator::hermitian_conjugate_inplace() {
		this->is_daggered = !(this->is_daggered);
		return *this;
	}
	Operator Operator::hermitian_conjugate()  const {
		Operator copy(*this);
		copy.hermitian_conjugate_inplace();
		return copy;
	}
	Operator Operator::with_momentum(Momentum const& new_momentum) const {
		assert(this->momentum.momentum_list.size() == 1U);
		Operator ret{ *this };
		ret.momentum = this->momentum.momentum_list.front().factor * new_momentum;
		return ret;
	}
	Operator Operator::with_momentum(const MomentumSymbol::name_type new_momentum) const {
		assert(this->momentum.momentum_list.size() == 1U);
		Operator ret{ *this };
		ret.momentum.momentum_list.front().name = new_momentum;
		return ret;
	}
	Operator Operator::add_momentum(Momentum const& to_add) const {
		Operator ret{ *this };
		ret.momentum += to_add;
		return ret;
	}
	Operator Operator::add_momentum(const MomentumSymbol::name_type to_add) const {
		Operator ret{ *this };
		ret.momentum += Momentum(to_add);
		return ret;
	}
	void Operator::remove_momentum_contribution(const MomentumSymbol::name_type value) {
		momentum.remove_contribution(value);
	}
	// Returns the first index, if the operator has an index.
	// Return Index::NoIndex otherwise
	Index Operator::first_index() const {
		return (indizes.empty() ? Index::NoIndex : indizes[0]);
	}
	// Sets the first index, if the operator has an index
	// Does nothing otherwise
	void Operator::set_first_index(Index index) {
		if (!indizes.empty()) {
			indizes[0] = index;
		}
	}
} // namespace mrock::symbolic_operators