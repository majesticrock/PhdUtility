#pragma once
#include "Momentum.hpp"
#include "IndexWrapper.hpp"

namespace SymbolicOperators {
	struct Operator {
		Momentum momentum;
		// Contains all indizes, standard: first index = spin, all others arbitrary, e.g., orbitals, bands etc
		IndexWrapper indizes;
		bool is_daggered{};
		bool is_fermion{ true };

		Operator(const Momentum& _momentum, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);
		Operator(const momentum_pairs& _momentum, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);
		Operator(char _momentum, bool add_Q, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);
		Operator(char _momentum, int sign, bool add_Q, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion = true);
		//Operator(char _momentum, int sign, bool add_Q, bool _is_daggered, bool _is_fermion = true);

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
		inline Index first_index() const {
			return (indizes.empty() ? Index::NoIndex : indizes[0]);
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