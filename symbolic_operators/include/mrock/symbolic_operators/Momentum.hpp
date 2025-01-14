#pragma once
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/string.hpp>
#include <algorithm>
#include <vector>
#include <utility>
#include <mrock/utility/VectorWrapper.hpp>
#include "MomentumSymbol.hpp"

namespace mrock::symbolic_operators {
	typedef std::vector<MomentumSymbol> momentum_symbols;
	struct Momentum {
		// total momentum is then:    sum_i pair_i.factor * pair_i.symbol
		momentum_symbols momentum_list;
		bool add_Q{};

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& momentum_list;
			ar& add_Q;
		}

		Momentum() : momentum_list(), add_Q(false) {};
		explicit Momentum(const char value, int plus_minus = 1, bool Q = false)
			: momentum_list(1, MomentumSymbol(plus_minus, value)), add_Q(Q) {};
		explicit Momentum(const MomentumSymbol::name_type value, int plus_minus = 1, bool Q = false)
			: momentum_list(1, {plus_minus, value}), add_Q(Q) {};
		explicit Momentum(const momentum_symbols& _momenta, bool Q = false)
			: momentum_list(_momenta), add_Q(Q) {};
		explicit Momentum(MomentumSymbol const& momentum_symbol, bool Q = false)
			: momentum_list{momentum_symbol}, add_Q(Q) {};
		Momentum(const std::string& expression, bool Q = false);
		Momentum(char, char) = delete;

		void sort();
		void remove_contribution(const MomentumSymbol::name_type momentum);

		Momentum& operator+=(const Momentum& rhs);
		Momentum& operator-=(const Momentum& rhs);
		inline Momentum& operator*=(const int rhs) {
			if (!(rhs & 1)) {
				this->add_Q = false;
			}
			for (auto& m : momentum_list) {
				m.factor *= rhs;
			}
			return *this;
		};

		inline void multiply_by(int factor) {
			(*this) *= factor;
		};
		inline void flip_momentum() {
			(*this) *= -1;
		};
		// Returns the position in the momentum_list array of the momentum given in <value>
		// Return -1 if the value is not contained in momentum_list
		inline int isUsed(const MomentumSymbol::name_type value) const noexcept {
			for (int i = 0; i < momentum_list.size(); ++i) {
				if (momentum_list[i].name == value) return i;
			}
			return -1;
		};
		inline bool differs_only_in_Q(Momentum rhs) const {
			if (rhs.add_Q == this->add_Q) return false;
			rhs.add_Q = this->add_Q;
			return (*this == rhs);
		};

		void add_in_place(const Momentum& rhs);
		void replace_occurances(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith);

		// removes entries within the momentum_list that have a 0 prefactor
		void remove_zeros();
		// replaces 'momentum' with -'momentum' if it exists within momentum_list
		void flip_single(const MomentumSymbol::name_type momentum);
		inline bool is_zero() const {
			if(add_Q) return false;
			return momentum_list.empty();
		}

		inline bool uses(const MomentumSymbol::name_type what) const {
			for (const auto& momentum_symbol : momentum_list) {
				if(momentum_symbol.name == what) return true;
			}
			return false;
		}

		inline bool first_momentum_is_negative() const {
			if (momentum_list.empty()) return false;
			return momentum_list.front().factor < 0;
		}
		inline bool first_momentum_is(const MomentumSymbol::name_type what) const {
			if (momentum_list.empty()) return false;
			return momentum_list.front().name == what;
		}
		inline bool last_momentum_is_negative() const {
			if (momentum_list.empty()) return false;
			return momentum_list.back().factor < 0;
		}
		inline bool last_momentum_is(const MomentumSymbol::name_type what) const {
			if (momentum_list.empty()) return false;
			return momentum_list.back().name == what;
		}

		bool operator==(const Momentum& rhs) const;
		inline bool operator!=(const Momentum& rhs) const {
			return !(*this == rhs);
		};

		std::string to_string() const;

		VECTOR_WRAPPER_FILL_MEMBERS(MomentumSymbol, momentum_list);
	};

	bool momentum_order(const Momentum& lhs, const Momentum& rhs);

	inline Momentum operator+(Momentum lhs, const Momentum& rhs) {
		lhs += rhs;
		return lhs;
	}
	inline Momentum operator-(Momentum lhs, const Momentum& rhs) {
		lhs -= rhs;
		return lhs;
	}
	inline Momentum operator*(Momentum lhs, const int rhs) {
		lhs *= rhs;
		return lhs;
	}
	inline Momentum operator*(const int lhs, Momentum rhs) {
		rhs *= lhs;
		return rhs;
	}
	inline Momentum operator-(Momentum rhs) {
		rhs.flip_momentum();
		return rhs;
	}

	bool operator>(const Momentum& lhs, const Momentum& rhs);
	bool operator<(const Momentum& lhs, const Momentum& rhs);

	std::ostream& operator<<(std::ostream& os, const Momentum& momentum);
}