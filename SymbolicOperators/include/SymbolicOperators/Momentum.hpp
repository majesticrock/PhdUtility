#pragma once
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/string.hpp>

namespace SymbolicOperators {
	typedef std::vector<std::pair<int, char>> momentum_pairs;
	struct Momentum {
		// total momentum is then:    sum_i pair_i.first * pair_i.second
		momentum_pairs momentum_list;
		bool add_Q{};

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& momentum_list;
			ar& add_Q;
		}

		Momentum() : momentum_list(), add_Q(false) {};
		explicit Momentum(const char value, int plus_minus = 1, bool Q = false)
			: momentum_list(1, std::make_pair(plus_minus, value)), add_Q(Q) {};
		explicit Momentum(const momentum_pairs& _momenta, bool Q = false)
			: momentum_list(_momenta), add_Q(Q) {};
		Momentum(const std::string& expression, bool Q = false);
		Momentum(char, char) = delete;

		void sort();

		Momentum& operator+=(const Momentum& rhs);
		Momentum& operator-=(const Momentum& rhs);
		inline Momentum& operator*=(const int rhs) {
			if (rhs % 2 == 0) {
				this->add_Q = false;
			}
			for (auto& m : momentum_list) {
				m.first *= rhs;
			}
			return *this;
		};

		inline void multiplyMomentum(int factor) {
			(*this) *= factor;
		};
		inline void flipMomentum() {
			(*this) *= -1;
		};
		// Returns the position in the momentum_list array of the momentum given in <value>
		// Return -1 if the value is not contained in momentum_list
		inline int isUsed(const char value) const noexcept {
			for (int i = 0; i < momentum_list.size(); ++i) {
				if (momentum_list[i].second == value) return i;
			}
			return -1;
		};
		inline bool differsOnlyInQ(Momentum rhs) const {
			if (rhs.add_Q == this->add_Q) return false;
			rhs.add_Q = this->add_Q;
			return (*this == rhs);
		};

		void addInPlace(const Momentum& rhs);
		void replaceOccurances(const char replaceWhat, const Momentum& replaceWith);

		// removes entries within the momentum_list that have a 0 prefactor
		void remove_zeros();
		// replaces 'momentum' with -'momentum' if it exists within momentum_list
		void flip_single(char momentum);

		bool operator==(const Momentum& rhs) const;
		inline bool operator!=(const Momentum& rhs) const {
			return !(*this == rhs);
		};
	};

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
		rhs.flipMomentum();
		return rhs;
	}

	bool operator>(const Momentum& lhs, const Momentum& rhs);
	bool operator<(const Momentum& lhs, const Momentum& rhs);

	std::ostream& operator<<(std::ostream& os, const Momentum& momentum);
}