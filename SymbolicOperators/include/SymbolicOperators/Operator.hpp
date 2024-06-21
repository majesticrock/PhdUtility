#pragma once
#include "Momentum.hpp"
#include "IndexWrapper.hpp"

namespace SymbolicOperators {
	struct Operator {
		Momentum momentum;
		// Contains all indizes, standard: first index = spin, all others arbitrary, e.g., orbitals, bands etc
		IndexWrapper indizes;
		bool isDaggered;

		Operator(const Momentum& _momentum, const IndexWrapper _indizes, bool _isDaggered);
		Operator(const momentum_pairs& _momentum, const IndexWrapper _indizes, bool _isDaggered);
		Operator(char _momentum, bool add_Q, const IndexWrapper _indizes, bool _isDaggered);
		Operator(char _momentum, int sign, bool add_Q, const IndexWrapper _indizes, bool _isDaggered);

		inline void hermitianConjugate() {
			this->isDaggered = !(this->isDaggered);
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
	};

	inline bool operator==(const Operator& lhs, const Operator& rhs) {
		if (lhs.isDaggered != rhs.isDaggered) return false;
		if (lhs.indizes != rhs.indizes) return false;
		return (lhs.momentum == rhs.momentum);
	}
	inline bool operator!=(const Operator& lhs, const Operator& rhs) {
		return !(lhs == rhs);
	}

	std::ostream& operator<<(std::ostream& os, const Operator& op);
	std::ostream& operator<<(std::ostream& os, const std::vector<Operator>& ops);
}