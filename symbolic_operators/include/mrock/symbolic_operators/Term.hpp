/**
 * @file Term.hpp
 * @brief Defines the Term class and related functions for symbolic operators.
 */

#pragma once

#include "KroneckerDelta.hpp"
#include "Coefficient.hpp"
#include "SumContainer.hpp"
#include <mrock/utility/Fractional.hpp>
#include <algorithm>
#include <vector>

// Debug option, allows to track an individual term across the entire algorithm
//#define _TRACK_TERM

#ifdef _TRACK_TERM
#define _TERM_TRACKER_PARAMETER bool is_tracked = false
#define _TERM_TRACKER_ATTRIBUTE bool is_tracked{};
#define IF_IS_TERM_TRACKED(statement) if ( is_tracked ) { statement ; }
#define CLEAR_TRACKED(terms) for (auto& __term__ : terms) { __term__.is_tracked = false; }
#else
#define _TERM_TRACKER_PARAMETER 
#define _TERM_TRACKER_ATTRIBUTE 
#define IF_IS_TERM_TRACKED(statement)
#define CLEAR_TRACKED(terms)
#endif

namespace mrock::symbolic_operators {
	using IntFractional = mrock::utility::Fractional<int>;

	/**
	 * @class Term
	 * @brief Represents a term in symbolic operator expressions.
	 */
	class Term {
	public:
		std::vector<Coefficient> coefficients; ///< Coefficients of the term.
		SumContainer sums; ///< Sum container for the term. Contains e.g. \sum_{k,l} \sum_{sigma}
		std::vector<Operator> operators; ///< Operators in the term, if empty the term is considered to contain the identiy operator
		std::vector<KroneckerDelta<Momentum>> delta_momenta; ///< Kronecker delta for momenta.
		std::vector<KroneckerDelta<Index>> delta_indizes; ///< Kronecker delta for indices.
		IntFractional multiplicity; ///< Multiplicity of the term.
		_TERM_TRACKER_ATTRIBUTE; ///< Attribute for tracking terms (if enabled).

		/**
		 * @brief Serializes the term.
		 * @tparam Archive The archive type.
		 * @param ar The archive.
		 * @param version The version.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& coefficients;
			ar& sums;
			ar& operators;
			ar& delta_momenta;
			ar& delta_indizes;
			ar& multiplicity;
		}

		friend struct WickTerm;

		/**
		 * @brief Constructs a Term with a summation over momenta and spins and multiple coefficients
		 * 
		 * @param _multiplicity The _multiplicity of the term
		 * @param _coefficients The coefficients
		 * @param _sums The sums
		 * @param _operators The operators of the term
		 */
		Term(IntFractional _multiplicity, std::vector<Coefficient> _coefficients, const SumContainer& _sums,
			const std::vector<Operator>& _operators = std::vector<Operator>());

		/**
		 * @brief Constructs a Term with a summation over momenta and spins and a coefficient
		 * 
		 * @param _multiplicity The _multiplicity of the term
		 * @param _coefficient The coefficient
		 * @param _sums The sums
		 * @param _operators The operators of the term
		 */
		Term(IntFractional _multiplicity, Coefficient _coefficient, const SumContainer& _sums,
			const std::vector<Operator>& _operators = std::vector<Operator>());

		/**
		 * @brief Constructs a Term with a summation over momenta and a coefficient
		 * 
		 * @param _multiplicity The _multiplicity of the term
		 * @param _coefficient The coefficient
		 * @param _sum_momenta Sum over momenta
		 * @param _operators The operators of the term
		 */
		Term(IntFractional _multiplicity, Coefficient _coefficient, const MomentumSum& _sum_momenta,
			const std::vector<Operator>& _operators = std::vector<Operator>());

		/**
		 * @brief Constructs a Term with a summation over spins (or other indizes) and a coefficient
		 * 
		 * @param _multiplicity The _multiplicity of the term
		 * @param _coefficient The coefficient
		 * @param _sum_spins Sum over spins (or other indizes)
		 * @param _operators The operators of the term
		 */
		Term(IntFractional _multiplicity, Coefficient _coefficient, const IndexSum& _sum_spins,
			const std::vector<Operator>& _operators = std::vector<Operator>());

		/**
		 * @brief Constructs a Term with a coefficient
		 * 
		 * @param _multiplicity The _multiplicity of the term
		 * @param _coefficient The coefficient
		 * @param _operators The operators of the term
		 */
		Term(IntFractional _multiplicity, Coefficient _coefficient, const std::vector<Operator>& _operators = std::vector<Operator>());

		/**
		 * @brief Constructs a Term with a summation over momenta and indizes
		 * 
		 * @param _multiplicity The _multiplicity of the term
		 * @param _sums Sums
		 * @param _operators The operators of the term
		 */
		Term(IntFractional _multiplicity, const SumContainer& _sums, const std::vector<Operator>& _operators = std::vector<Operator>());

		/**
		 * @brief Constructs a Term with a summation over momenta
		 * 
		 * @param _multiplicity The _multiplicity of the term
		 * @param _sum_momenta Sum over momenta
		 * @param _operators The operators of the term
		 */
		Term(IntFractional _multiplicity, const MomentumSum& _sum_momenta, const std::vector<Operator>& _operators = std::vector<Operator>());

		/**
		 * @brief Constructs a Term with a summation over spins (or other indizes)
		 * 
		 * @param _multiplicity The _multiplicity of the term
		 * @param _sum_spins Sum over spins (or other indizes)
		 * @param _operators The operators of the term
		 */
		Term(IntFractional _multiplicity, const IndexSum& _sum_spins, const std::vector<Operator>& _operators = std::vector<Operator>());

		/**
		 * @brief Constructs a Term with only a multiplicity
		 * 
		 * @param _multiplicity The _multiplicity of the term
		 * @param _operators The operators of the term
		 */
		explicit Term(IntFractional _multiplicity, const std::vector<Operator>& _operators = std::vector<Operator>());
		/**
		 * @brief Default constructor.
		 */
		Term() = default;

		/**
		 * @brief Checks if the term is an identity.
		 * @return True if the term is an identity, false otherwise.
		 */
		inline bool is_identity() const;

		/**
		 * @brief Checks if the term contains a boson.
		 * @return True if the term contains a boson, false otherwise.
		 */
		inline bool contains_boson() const;

		/**
		 * @brief Checks if the term contains a fermion.
		 * @return True if the term contains a fermion, false otherwise.
		 */
		inline bool contains_fermion() const;

		/**
		 * @brief Counts the number of bosons in the term.
		 * @return The number of bosons.
		 */
		inline int count_bosons() const;

		/**
		 * @brief Counts the number of fermions in the term.
		 * @return The number of fermions.
		 */
		inline int count_fermions() const;

		/**
		 * @brief Prints the term.
		 */
		void print() const;

		/**
		 * @brief Flips the sign of the term.
		 */
		inline void flip_sign();

		/**
		 * @brief Swaps two operators in the term. Does NOT consider possible additional terms spawned by this operation due to non-commutivity!
		 * @param lhs The first operator.
		 * @param rhs The second operator.
		 */
		inline void perform_operator_swap(Operator& lhs, Operator& rhs);

		/**
		 * @brief Gets the operators in the term.
		 * @return The operators.
		 */
		inline const std::vector<Operator>& get_operators() const;

		/**
		 * @brief Sets the Kronecker deltas in the term.
		 * @return True if successful, false otherwise.
		 */
		bool set_deltas();

		/**
		 * @brief Computes the sums in the term.
		 * @return True if successful, false otherwise.
		 */
		bool compute_sums();

		/**
		 * @brief Discards zero momenta in the term.
		 */
		void discard_zero_momenta();

		/**
		 * @brief Sorts the term.
		 */
		void sort();

		/**
		 * @brief Renames the sum indices in the term.
		 */
		void rename_sums();

		/**
		 * @brief Checks if the term is equal to another term (excluding multiplicity).
		 * @param other The other term.
		 * @return True if equal, false otherwise.
		 */
		bool is_equal(const Term& other) const;

		/**
		 * @brief Checks if the term is in normal order.
		 * @return True if in normal order, false otherwise.
		 */
		bool is_normal_ordered() const;

		/**
		 * @brief Converts the term to a string without the prefactor.
		 * @return The string representation.
		 */
		std::string to_string_without_prefactor() const;

		/**
		 * @brief Applies the Hermitian conjugate to the term.
		 * @return A reference to *this
		 */
		Term& hermitian_conjugate_inplace();

		/**
		 * @brief Creates hermitian conjugate of this as a new object.
		 * @return Returns the new object.
		 */
		Term hermitian_conjugate() const;

		/**
		 * @brief Renames indices in the term.
		 * @param what The index to rename.
		 * @param to The new index.
		 */
		void rename_indizes(const Index what, const Index to);

		/**
		 * @brief Renames momenta in the term.
		 * @param what The momentum to rename.
		 * @param to The new momentum.
		 */
		void rename_momenta(const MomentumSymbol::name_type what, const MomentumSymbol::name_type to);

		/**
		 * @brief Swaps two momenta in the term.
		 * @param a The first momentum.
		 * @param b The second momentum.
		 */
		void swap_momenta(const MomentumSymbol::name_type a, const MomentumSymbol::name_type b);

		/**
		 * @brief Transforms a momentum sum in the term.
		 * @param what The momentum to transform.
		 * @param to The new momentum.
		 * @param new_sum_index The new sum index.
		 */
		void transform_momentum_sum(const MomentumSymbol::name_type what, const Momentum to, const MomentumSymbol::name_type new_sum_index);

		/**
		 * @brief Inverts a momentum in the term.
		 * @param what The momentum to invert.
		 */
		void invert_momentum(const MomentumSymbol::name_type what);

		/**
		 * @brief Inverts a momentum sum in the term.
		 * @param what The momentum to invert.
		 */
		void invert_momentum_sum(const MomentumSymbol::name_type what);

		/**
		 * @brief Removes a momentum contribution from the term.
		 * @param value The momentum to remove.
		 */
		void remove_momentum_contribution(const MomentumSymbol::name_type value);

		friend void normal_order(std::vector<Term>& terms);
		friend std::vector<Term> commutator(const Term& left, const Term& right);
		friend std::ostream& operator<<(std::ostream& os, const Term& term);
	};

	/**
	 * @brief Computes the commutator of two sets of terms.
	 * @param left The left-hand side terms.
	 * @param right The right-hand side terms.
	 * @return The result of [left, right]
	 */
	std::vector<Term> commutator(const std::vector<Term>& left, const std::vector<Term>& right);

	/**
	 * @brief Computes the commutator of a term and a set of terms.
	 * @param left The left-hand side term.
	 * @param right The right-hand side terms.
	 * @return The result of [left, right]
	 */
	inline std::vector<Term> commutator(const Term& left, const std::vector<Term>& right);

	/**
	 * @brief Computes the commutator of a set of terms and a term.
	 * @param left The left-hand side terms.
	 * @param right The right-hand side term.
	 * @return The result of [left, right]
	 */
	inline std::vector<Term> commutator(const std::vector<Term>& left, const Term& right);

	/**
	 * @brief Checks if two terms are equal.
	 * @param lhs The left-hand side term.
	 * @param rhs The right-hand side term.
	 * @return True if equal, false otherwise.
	 */
	inline bool operator==(const Term& lhs, const Term& rhs) { return lhs.is_equal(rhs); }

	/**
	 * @brief Checks if two terms are not equal.
	 * @param lhs The left-hand side term.
	 * @param rhs The right-hand side term.
	 * @return True if not equal, false otherwise.
	 */
	inline bool operator!=(const Term& lhs, const Term& rhs) { return !(lhs == rhs); }

	/**
	 * @brief Outputs a coefficient to a stream.
	 * @param os The output stream.
	 * @param coeff The coefficient.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const Coefficient& coeff);

	/**
	 * @brief Outputs a vector of coefficients to a stream.
	 * @param os The output stream.
	 * @param coeffs The coefficients.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const std::vector<Coefficient>& coeffs);

	/**
	 * @brief Outputs a vector of terms to a stream.
	 * @param os The output stream.
	 * @param terms The terms.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const std::vector<Term>& terms);

	/**
	 * @brief Clears duplicate terms from a vector.
	 * @param terms The vector of terms.
	 */
	void clear_duplicates(std::vector<Term>& terms);

	/**
	 * @brief Cleans up a vector of terms.
	 * @param terms The vector of terms.
	 */
	void clean_up(std::vector<Term>& terms);

	/**
	 * @brief Applies the Hermitian conjugate to a vector of terms.
	 * @param terms The vector of terms.
	 */
	inline void hermitian_conjugate(std::vector<Term>& terms);

	/**
	 * @brief Renames momenta in a vector of terms.
	 * @param terms The vector of terms.
	 * @param what The momentum to rename.
	 * @param to The new momentum.
	 */
	inline void rename_momenta(std::vector<Term>& terms, const MomentumSymbol::name_type what, const MomentumSymbol::name_type to);

	/**
	 * @brief Converts a vector of terms to a string without the prefactor.
	 * @param terms The vector of terms.
	 * @return The string representation.
	 */
	std::string to_string_without_prefactor(const std::vector<Term>& terms);

	// Inline definitions
	bool Term::is_identity() const {
		return this->operators.empty();
	}
	bool Term::contains_boson() const {
		return std::any_of(operators.begin(), operators.end(), [](Operator const& op) { return (!op.is_fermion); });
	}
	bool Term::contains_fermion() const {
		return std::any_of(operators.begin(), operators.end(), [](Operator const& op) { return op.is_fermion; });
	}
	int Term::count_bosons() const {
		return std::count_if(operators.begin(), operators.end(), [](Operator const& op) { return (!op.is_fermion); });
	}
	int Term::count_fermions() const {
		return std::count_if(operators.begin(), operators.end(), [](Operator const& op) { return op.is_fermion; });
	}
	void Term::flip_sign() {
		this->multiplicity *= -1;
	}
	void Term::perform_operator_swap(Operator& lhs, Operator& rhs) {
		if (lhs.is_fermion && rhs.is_fermion) {
			this->multiplicity *= -1;
		}
		std::swap(lhs, rhs);
	}
	const std::vector<Operator>& Term::get_operators() const {
		return this->operators;
	}

	// Non-member inlines
	std::vector<Term> commutator(const Term& left, const std::vector<Term>& right) {
		const std::vector<Term> buffer = { left };
		return commutator(buffer, right);
	}
	std::vector<Term> commutator(const std::vector<Term>& left, const Term& right) {
		const std::vector<Term> buffer = { right };
		return commutator(left, buffer);
	}
	void hermitian_conjugate(std::vector<Term>& terms) {
		for (auto& t : terms) {
			t.hermitian_conjugate_inplace();
		}
	}
	void rename_momenta(std::vector<Term>& terms, const MomentumSymbol::name_type what, const MomentumSymbol::name_type to) {
		for (auto& t : terms) {
			t.rename_momenta(what, to);
		}
	}
} // namespace mrock::symbolic_operators