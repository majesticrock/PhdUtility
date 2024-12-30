#pragma once
#include "Momentum.hpp"
#include "IndexWrapper.hpp"

namespace mrock::symbolic_operators {
	struct Operator {
		Momentum momentum;
		// Contains all indizes, standard: first index = spin, all others arbitrary, e.g., orbitals, bands etc
		IndexWrapper indizes;
		bool is_daggered{};
		bool is_fermion{ true };

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& momentum;
			ar& indizes;
			ar& is_daggered;
			ar& is_fermion;
		}

		Operator() = default;
		Operator(const Momentum& _momentum, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);
		Operator(const momentum_pairs& _momentum, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);
		Operator(char _momentum, bool add_Q, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);
		Operator(char _momentum, int sign, bool add_Q, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);

		inline static Operator Boson(const Momentum& _momentum, const IndexWrapper _indizes, bool _is_daggered) {
			return Operator(_momentum, _indizes, _is_daggered, false);
		}
		inline static Operator Boson(const Momentum& _momentum, bool _is_daggered) {
			return Operator(_momentum, IndexWrapper{}, _is_daggered, false);
		}

		inline void hermitianConjugate() {
			this->is_daggered = !(this->is_daggered);
		};
		inline Operator with_momentum(Momentum const& new_momentum) const {
			assert(this->momentum.momentum_list.size() == 1U);
			Operator ret{ *this };
			ret.momentum = this->momentum.momentum_list.front().first * new_momentum;
			return ret;
		};
		inline Operator with_momentum(char new_momentum) const {
			assert(this->momentum.momentum_list.size() == 1U);
			Operator ret{ *this };
			ret.momentum.momentum_list.front().second = new_momentum;
			return ret;
		};
		inline Operator add_momentum(Momentum const& to_add) const {
			Operator ret{ *this };
			ret.momentum += to_add;
			return ret;
		}
		inline Operator add_momentum(char to_add) const {
			Operator ret{ *this };
			ret.momentum += Momentum(to_add);
			return ret;
		}
		inline void remove_momentum_contribution(char value) {
			momentum.remove_contribution(value);
		}
		// Returns the first index, if the operator has an index.
		// Return Index::NoIndex otherwise
		inline Index first_index() const {
			return (indizes.empty() ? Index::NoIndex : indizes[0]);
		}
		// Sets the first index, if the operator has an index
		// Does nothing otherwise
		inline void set_first_index(Index index) {
			if (!indizes.empty()) {
				indizes[0] = index;
			}
		}
	};

	inline bool operator==(const Operator& lhs, const Operator& rhs) {
		if (lhs.is_fermion != rhs.is_fermion) return false;
		if (lhs.is_daggered != rhs.is_daggered) return false;
		if (lhs.indizes != rhs.indizes) return false;
		return (lhs.momentum == rhs.momentum);
	}
	inline bool operator!=(const Operator& lhs, const Operator& rhs) {
		return !(lhs == rhs);
	}

	std::ostream& operator<<(std::ostream& os, const Operator& op);
	std::ostream& operator<<(std::ostream& os, const std::vector<Operator>& ops);
}