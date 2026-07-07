/**
 * @file Term.hpp
 * @brief Defines the Term class and related functions for symbolic operators.
 */

#pragma once

#include "AbstractTerm.hpp"
#include <algorithm>

namespace mrock::symbolic_operators {
	/**
	 * @class Term
	 * @brief Represents a term in symbolic operator expressions.
	 * 
	 * This class represents a Term. It has various kind of constructors that allow setting coefficient(s), sums, operators and deltas.
	 * Using IntFractional, the term can have rational prefactors, e.g., 1/2.
	 * 
	 * A Hamiltonian (or any other summation of operators) is characterized as \c std::vector<Term>.
	 * It can consist of any number of individual terms. 
	 * For a few practical examples, see the files in the tests folder.
	 * See bosons.cpp, continuum.cpp, and compare_test.hpp.
	 * My own projects using this library are, e.g., 
	 * https://github.com/majesticrock/FermionCommute and https://github.com/majesticrock/FlowCommutators.
	 * 
	 * After creating atleast two \c Terms (or \c std::vector<Term>), you may commute them by calling
	 * \code
	 * std::vector<Term> result = commutator(A, B);
	 * clean_up(result);
	 * \endcode
	 * After calling the commutator, you should pretty much always call mrock::symbolic_operators::clean_up(std::vector<Term>)
	 * because commutator performs the normal ordering procedure, however, does not attempt to beautify the result.
	 * clean_up then sorts the terms, adds identical ones togeFther and removes those that are equal to 0.
	 * 
	 * Similarly, a double commutator \f$ [C, [A, B]] \f$ can be evaluated by
	 * \code
	 * std::vector<Term> inner_result = commutator(A, B);
	 * clean_up(inner_result);
	 * std::vector<Term> result = commutator(C, inner_result);
	 * clean_up(result);
	 * \endcode
	 * 
	 * To output the results, an overload of \c operator<< is provided for both \c Term and \c std::vector<Term>.
	 * The out put is formatted so that it can be used within an align-environment within LaTeX.
	 * 
	 * @sa Coefficient, SumContainer, Operator, KroneckerDelta
	 */
	class Term : public AbstractTerm<Operator> {
	public:
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
		 * @brief Swaps two operators in the term. Does NOT consider possible additional terms spawned by this operation due to non-commutivity!
		 * @param lhs The first operator.
		 * @param rhs The second operator.
		 */
		inline void perform_operator_swap(Operator& lhs, Operator& rhs);

		/**
		 * @brief Resolves the Kronecker deltas in the term ( calls \c resolve_momentum_deltas() and \c resolve_index_deltas() )
		 * @return True if successful, false otherwise.
		 */
		bool resolve_deltas();

		/**
		 * @brief Discards zero momenta in the term.
		 */
		void discard_zero_momenta();

		/**
		 * @brief Sorts the term.
		 */
		void sort();

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
		 * @brief Normal orders the terms by using the canoncical (anti-)commutation relations
		 * The result is stored in the input vector.
		 * A simple example is
		 * \f$ b b^\dagger = 1 \pm b^\dagger b \f$, where the + applies to bosons and the minus to fermions.
		 * @param terms The terms to normal order.
		 */
		friend void normal_order(std::vector<Term>& terms);

		/**
		 * @brief Computes the commutator of two terms: \f$ [A, B] = AB - BA \f$.
		 * @param left The left term.
		 * @param right The right term.
		 * @return The commutation result.
		 */
		friend std::vector<Term> commutator(const Term& left, const Term& right);

		/**
	 	 * @brief Overloads the stream insertion operator for the Term class.
	 	 * 
	 	 * @param os The output stream.
	 	 * @param term The Term object to insert into the stream.
	 	 * @return The output stream.
	 	 */
		friend std::ostream& operator<<(std::ostream& os, const Term& term);

		/**
		 * @brief Multiplies this by rhs
		 * IMPOARTANT: The result will not be normal ordered! If you require
		 * a normal ordered expression, please call normal_order!
		 * Note, that doing so may create additional terms, and the result must therefore
		 * be std::vector<Term> and cannot be implemented as Term::operator*=
		 * 
		 * @param rhs the right-hand side Term object.
		 * @return Reference to this containing the result.
		 */
		Term& operator*=(const Term& rhs);
	};

	/**
	 * @brief Computes the commutator of two sets of terms: \f$ [A, B] = AB - BA \f$.
	 * @param left The left-hand side terms.
	 * @param right The right-hand side terms.
	 * @return The result of [left, right]
	 */
	std::vector<Term> commutator(const std::vector<Term>& left, const std::vector<Term>& right);

	/**
	 * @brief Computes the commutator of a term and a set of terms: \f$ [A, B] = AB - BA \f$.
	 * @param left The left-hand side term.
	 * @param right The right-hand side terms.
	 * @return The result of [left, right]
	 */
	inline std::vector<Term> commutator(const Term& left, const std::vector<Term>& right);

	/**
	 * @brief Computes the commutator of a set of terms and a term: \f$ [A, B] = AB - BA \f$.
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
	 * @brief Multiplies a Term by another Term.
	 * IMPOARTANT: The result will not be normal ordered! If you require
	 * a normal ordered expression, please call normal_order!
	 * See Term::operator*=(const Term& rhs)
	 * 
	 * @param lhs The left-hand side Term.
	 * @param rhs The right-hand side Term.
	 * @return The result of the multiplication.
	 */
	inline Term operator*(Term lhs, const Term& rhs) {
		lhs *= rhs;
		return lhs;
	}

	/**
	 * @brief Negates each term in the vector.
	 *
	 * @param terms The terms to negate.
	 * @return A vector containing the negated terms.
	 */
	std::vector<Term> operator-(std::vector<Term> terms);

	/**
	 * @brief Adds the right-hand side terms to the left-hand side terms.
	 *
	 * @param lhs The left-hand side terms that are modified in place.
	 * @param rhs The right-hand side terms to add.
	 * @return A reference to the modified left-hand side terms.
	 */
	std::vector<Term>& operator+=(std::vector<Term>& lhs, const std::vector<Term>& rhs);
	
	/**
	 * @brief Subtracts the right-hand side terms from the left-hand side terms.
	 *
	 * @param lhs The left-hand side terms that are modified in place.
	 * @param rhs The right-hand side terms to subtract.
	 * @return A reference to the modified left-hand side terms.
	 */
	std::vector<Term>& operator-=(std::vector<Term>& lhs, const std::vector<Term>& rhs);

	/**
	 * @brief Multiplies the left-hand side term vector by the right-hand side term vector.
	 *
	 * @param lhs The left-hand side terms that are modified in place.
	 * @param rhs The right-hand side terms to multiply by.
	 * @return A reference to the modified left-hand side terms.
	 */
	std::vector<Term>& operator*=(std::vector<Term>& lhs, const std::vector<Term>& rhs);

	/**
	 * @brief Adds two vectors of terms and returns the result.
	 *
	 * @param lhs The left-hand side terms.
	 * @param rhs The right-hand side terms.
	 * @return A vector containing the sum of the operands.
	 */
	inline std::vector<Term> operator+(std::vector<Term> lhs, const std::vector<Term>& rhs) {
		lhs += rhs;
		return lhs;
	}

	/**
	 * @brief Subtracts the right-hand side term vector from the left-hand side term vector.
	 *
	 * @param lhs The left-hand side terms.
	 * @param rhs The right-hand side terms.
	 * @return A vector containing the difference of the operands.
	 */
	inline std::vector<Term> operator-(std::vector<Term> lhs, const std::vector<Term>& rhs) {
		lhs -= rhs;
		return lhs;
	}

	/**
	 * @brief Multiplies two vectors of terms and returns the result.
	 *
	 * @param lhs The left-hand side terms.
	 * @param rhs The right-hand side terms.
	 * @return A vector containing the product of the operands.
	 */
	inline std::vector<Term> operator*(std::vector<Term> lhs, const std::vector<Term>& rhs) {
		lhs *= rhs;
		return lhs;
	}

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
	 * @brief Sorts the terms, adds identical ones together and removes those that are equal to 0.
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
	void Term::perform_operator_swap(Operator& lhs, Operator& rhs) {
		if (lhs.is_fermion && rhs.is_fermion) {
			this->multiplicity *= -1;
		}
		std::swap(lhs, rhs);
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