/**
 * @file Momentum.hpp
 * @brief Defines the Momentum structure and related operations for symbolic manipulation of momentum symbols.
 */

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

	/**
	 * @typedef momentum_symbols
	 * @brief Alias for a vector of MomentumSymbol.
	 */
	typedef std::vector<MomentumSymbol> momentum_symbols;

	/**
	 * @struct Momentum
	 * @brief Represents a collection of momentum symbols with associated operations.
	 */
	struct Momentum {
		momentum_symbols momentum_list; ///< List of momentum symbols.
		bool add_Q{}; ///< Flag indicating additional property Q. Q is a special momentum with the property 2Q = 0 (remeber that momenta are only defined in the first Brillouin zone)

		/**
		 * @brief Serialization function for Boost.
		 * @tparam Archive Type of the archive.
		 * @param ar Archive to serialize to/from.
		 * @param version Version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& momentum_list;
			ar& add_Q;
		}

		/**
		 * @brief Default constructor.
		 */
		Momentum() = default;

		/**
		 * @brief Constructs a Momentum with a single symbol.
		 * @param value Character representing the symbol.
		 * @param plus_minus Factor associated with the symbol.
		 * @param Q Additional property Q.
		 */
		explicit Momentum(const char value, int plus_minus = 1, bool Q = false);

		/**
		 * @brief Constructs a Momentum with a single symbol.
		 * @param value Name type of the symbol.
		 * @param plus_minus Factor associated with the symbol.
		 * @param Q Additional property Q.
		 */
		explicit Momentum(const MomentumSymbol::name_type value, int plus_minus = 1, bool Q = false);

		/**
		 * @brief Constructs a Momentum with a list of symbols.
		 * @param _momenta List of momentum symbols.
		 * @param Q Additional property Q.
		 */
		explicit Momentum(const momentum_symbols& _momenta, bool Q = false);

		/**
		 * @brief Constructs a Momentum with a single symbol.
		 * @param momentum_symbol The momentum symbol.
		 * @param Q Additional property Q.
		 */
		explicit Momentum(MomentumSymbol const& momentum_symbol, bool Q = false);

		/**
		 * @brief Constructs a Momentum from a string expression.
		 * @param expression String representing the momentum expression.
		 * @param Q Additional property Q.
		 */
		Momentum(const std::string& expression, bool Q = false);

		/**
		 * @brief Deleted constructor to prevent usage.
		 */
		Momentum(char, char) = delete;

		/**
		 * @brief Sorts the momentum symbols.
		 */
		void sort();

		/**
		 * @brief Removes a specific momentum contribution.
		 * @param momentum Name of the momentum to remove.
		 */
		void remove_contribution(const MomentumSymbol::name_type momentum);

		/**
		 * @brief Adds another Momentum in place.
		 * @param rhs The other Momentum.
		 */
		void add_in_place(const Momentum& rhs);

		/**
		 * @brief Replaces occurrences of a specific momentum with another Momentum.
		 * @param replaceWhat Name of the momentum to replace.
		 * @param replaceWith The Momentum to replace with.
		 */
		void replace_occurances(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith);

		/**
		 * @brief Removes entries with a zero prefactor.
		 */
		void remove_zeros();

		/**
		 * @brief Flips a specific momentum if it exists.
		 * @param momentum Name of the momentum to flip.
		 */
		void flip_single(const MomentumSymbol::name_type momentum);

		/**
		 * @brief Checks if a specific momentum is used.
		 * @param value Name of the momentum.
		 * @return Position of the momentum in the list, or -1 if not found.
		 */
		int is_used_at(const MomentumSymbol::name_type value) const noexcept;

		/**
		 * @brief Multiplies this Momentum by an integer factor.
		 * @param factor The factor.
		 */
		inline void multiply_by(int factor);

		/**
		 * @brief Flips the momentum by multiplying by -1.
		 */
		inline void flip_momentum();

		/**
		 * @brief Checks if this Momentum differs from another only in the Q property.
		 * @param rhs The other Momentum.
		 * @return True if they differ only in Q, false otherwise.
		 */
		inline bool differs_only_in_Q(Momentum rhs) const;

		/**
		 * @brief Checks if this Momentum is zero.
		 * @return True if zero, false otherwise.
		 */
		inline bool is_zero() const;

		/**
		 * @brief Checks if a specific momentum is used.
		 * @param what Name of the momentum.
		 * @return True if used, false otherwise.
		 */
		inline bool uses(const MomentumSymbol::name_type what) const noexcept;

		/**
		 * @brief Checks if the first momentum is negative.
		 * @return True if negative, false otherwise.
		 */
		inline bool first_momentum_is_negative() const;

		/**
		 * @brief Checks if the first momentum is a specific value.
		 * @param what Name of the momentum.
		 * @return True if it matches, false otherwise.
		 */
		inline bool first_momentum_is(const MomentumSymbol::name_type what) const;

		/**
		 * @brief Checks if the last momentum is negative.
		 * @return True if negative, false otherwise.
		 */
		inline bool last_momentum_is_negative() const;

		/**
		 * @brief Checks if the last momentum is a specific value.
		 * @param what Name of the momentum.
		 * @return True if it matches, false otherwise.
		 */
		inline bool last_momentum_is(const MomentumSymbol::name_type what) const;

		/**
		 * @brief Converts this Momentum to a string representation.
		 * @return String representation of this Momentum.
		 */
		std::string to_string() const;

		/**
		 * @brief Equality operator.
		 * @param rhs The other Momentum.
		 * @return True if equal, false otherwise.
		 */
		bool operator==(const Momentum& rhs) const;

		/**
		 * @brief Inequality operator.
		 * @param rhs The other Momentum.
		 * @return True if not equal, false otherwise.
		 */
		inline bool operator!=(const Momentum& rhs) const;

		/**
		 * @brief Adds another Momentum to this one.
		 * @param rhs The other Momentum.
		 * @return Reference to this Momentum.
		 */
		Momentum& operator+=(const Momentum& rhs);

		/**
		 * @brief Subtracts another Momentum from this one.
		 * @param rhs The other Momentum.
		 * @return Reference to this Momentum.
		 */
		Momentum& operator-=(const Momentum& rhs);

		/**
		 * @brief Multiplies this Momentum by an integer factor.
		 * @param rhs The factor.
		 * @return Reference to this Momentum.
		 */
		inline Momentum& operator*=(const int rhs);

		VECTOR_WRAPPER_FILL_MEMBERS(MomentumSymbol, momentum_list);
	};

	/**
	 * @brief Compares two Momentum objects for ordering.
	 * @param lhs The left-hand side Momentum.
	 * @param rhs The right-hand side Momentum.
	 * @return True if lhs is ordered before rhs, false otherwise.
	 */
	bool momentum_order(const Momentum& lhs, const Momentum& rhs);

	/**
	 * @brief Adds two Momentum objects.
	 * @param lhs The left-hand side Momentum.
	 * @param rhs The right-hand side Momentum.
	 * @return The result of the addition.
	 */
	inline Momentum operator+(Momentum lhs, const Momentum& rhs) {
		lhs += rhs;
		return lhs;
	}

	/**
	 * @brief Subtracts one Momentum from another.
	 * @param lhs The left-hand side Momentum.
	 * @param rhs The right-hand side Momentum.
	 * @return The result of the subtraction.
	 */
	inline Momentum operator-(Momentum lhs, const Momentum& rhs) {
		lhs -= rhs;
		return lhs;
	}

	/**
	 * @brief Multiplies a Momentum by an integer factor.
	 * @param lhs The Momentum.
	 * @param rhs The factor.
	 * @return The result of the multiplication.
	 */
	inline Momentum operator*(Momentum lhs, const int rhs) {
		lhs *= rhs;
		return lhs;
	}

	/**
	 * @brief Multiplies an integer factor by a Momentum.
	 * @param lhs The factor.
	 * @param rhs The Momentum.
	 * @return The result of the multiplication.
	 */
	inline Momentum operator*(const int lhs, Momentum rhs) {
		rhs *= lhs;
		return rhs;
	}

	/**
	 * @brief Negates a Momentum.
	 * @param rhs The Momentum.
	 * @return The negated Momentum.
	 */
	inline Momentum operator-(Momentum rhs) {
		rhs.flip_momentum();
		return rhs;
	}

	/**
	 * @brief Compares two Momentum objects for greater-than ordering.
	 * @param lhs The left-hand side Momentum.
	 * @param rhs The right-hand side Momentum.
	 * @return True if lhs is greater than rhs, false otherwise.
	 */
	bool operator>(const Momentum& lhs, const Momentum& rhs);

	/**
	 * @brief Compares two Momentum objects for less-than ordering.
	 * @param lhs The left-hand side Momentum.
	 * @param rhs The right-hand side Momentum.
	 * @return True if lhs is less than rhs, false otherwise.
	 */
	bool operator<(const Momentum& lhs, const Momentum& rhs);

	/**
	 * @brief Outputs a Momentum to an output stream.
	 * @param os The output stream.
	 * @param momentum The Momentum.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const Momentum& momentum);


	// Inline definitions
	Momentum& Momentum::operator*=(const int rhs) {
		if (!(rhs & 1)) {
			this->add_Q = false;
		}
		for (auto& m : momentum_list) {
			m.factor *= rhs;
		}
		return *this;
	}
	void Momentum::multiply_by(int factor) {
		(*this) *= factor;
	}
	void Momentum::flip_momentum() {
		(*this) *= -1;
	}
	bool Momentum::differs_only_in_Q(Momentum rhs) const {
		if (rhs.add_Q == this->add_Q) return false;
		rhs.add_Q = this->add_Q;
		return (*this == rhs);
	}
	bool Momentum::is_zero() const {
		if(add_Q) return false;
		return momentum_list.empty();
	}
	bool Momentum::uses(const MomentumSymbol::name_type what) const noexcept {
		return is_used_at(what) != -1;
	}
	bool Momentum::first_momentum_is_negative() const {
		if (momentum_list.empty()) return false;
		return momentum_list.front().factor < 0;
	}
	bool Momentum::first_momentum_is(const MomentumSymbol::name_type what) const {
		if (momentum_list.empty()) return false;
		return momentum_list.front().name == what;
	}
	bool Momentum::last_momentum_is_negative() const {
		if (momentum_list.empty()) return false;
		return momentum_list.back().factor < 0;
	}
	bool Momentum::last_momentum_is(const MomentumSymbol::name_type what) const {
		if (momentum_list.empty()) return false;
		return momentum_list.back().name == what;
	}
	bool Momentum::operator!=(const Momentum& rhs) const {
		return !(*this == rhs);
	}
} // namespace mrock::symbolic_operators