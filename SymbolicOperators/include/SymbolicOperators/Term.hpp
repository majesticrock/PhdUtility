#pragma once
#include "KroneckerDelta.hpp"
#include "Coefficient.hpp"
#include "SymbolicSum.hpp"
#include <Utility/Fractional.hpp>
#include <algorithm>
#include <vector>

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

namespace SymbolicOperators {
	using IntFractional = Utility::Fractional<int>;

	class Term {
	public:
		std::vector<Coefficient> coefficients;
		SumContainer sums;
		std::vector<Operator> operators;
		std::vector<KroneckerDelta<Momentum>> delta_momenta;
		std::vector<KroneckerDelta<Index>> delta_indizes;
		IntFractional multiplicity;
		_TERM_TRACKER_ATTRIBUTE;

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
		Term(IntFractional _multiplicity, std::vector<Coefficient> _coefficients, const SumContainer& _sums,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(IntFractional _multiplicity, Coefficient _coefficient, const SumContainer& _sums,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(IntFractional _multiplicity, Coefficient _coefficient, const MomentumSum& _sum_momenta,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(IntFractional _multiplicity, Coefficient _coefficient, const IndexSum& _sum_spins,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(IntFractional _multiplicity, Coefficient _coefficient,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(IntFractional _multiplicity, const SumContainer& _sums,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(IntFractional _multiplicity, const MomentumSum& _sum_momenta,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term(IntFractional _multiplicity, const IndexSum& _sum_spins,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		explicit Term(IntFractional _multiplicity,
			const std::vector<Operator>& _operators = std::vector<Operator>());
		Term() = default;

		inline bool isIdentity() const {
			return this->operators.empty();
		}
		inline bool contains_boson() const {
			return std::any_of(operators.begin(), operators.end(), [](Operator const& op) { return (!op.is_fermion); });
		}
		inline bool contains_fermion() const {
			return std::any_of(operators.begin(), operators.end(), [](Operator const& op) { return op.is_fermion; });
		}
		inline int count_bosons() const {
			return std::count_if(operators.begin(), operators.end(), [](Operator const& op) { return (!op.is_fermion); });
		}
		inline int count_fermions() const {
			return std::count_if(operators.begin(), operators.end(), [](Operator const& op) { return op.is_fermion; });
		}
		
		void print() const;
		inline void flipSign() {
			this->multiplicity *= -1;
		}
		// Swaps the two operators lhs and rhs. 
		// Does NOT consider possible additional terms spawned by this operation due to non-commutivity!
		inline void perform_operator_swap(Operator& lhs, Operator& rhs) {
			if (lhs.is_fermion && rhs.is_fermion) {
				this->multiplicity *= -1;
			}
			std::swap(lhs, rhs);
		}
		inline const std::vector<Operator>& getOperators() const {
			return this->operators;
		}

		bool setDeltas();
		bool computeSums();
		void discardZeroMomenta();
		void sort();
		// Unifies the sum indizes
		void renameSums();
		//void wick(WickTermCollector& reciever) const;

		// Checks for equality of everything except of multiplicity
		inline bool is_equal(const Term& other) const {
			if (this->coefficients != other.coefficients) return false;
			if (this->sums != other.sums) return false;
			if (this->delta_indizes != other.delta_indizes) return false;
			if (this->delta_momenta != other.delta_momenta) return false;
			if (this->operators != other.operators) return false;
			return true;
		};

		bool is_normal_ordered() const;

		std::string toStringWithoutPrefactor() const;

		inline void hermitianConjugate() {
			std::reverse(this->operators.begin(), this->operators.end());
			for (auto& op : this->operators) {
				op.hermitianConjugate();
			}
		}

		void rename_indizes(Index what, Index to);
		void rename_momenta(char what, char to);

		inline void swap_momenta(char a, char b) {
			this->rename_momenta(a, '_');
			this->rename_momenta(b, a);
			this->rename_momenta('_', b);
		}

		void transform_momentum_sum(char what, Momentum to, char new_sum_index);
		// Inverts a momenta, e.g., q -> -q
		void invert_momentum(char what);
		// Same as invert_momentum, but performs a check, whether 'what' is actually being summed over
		void invert_momentum_sum(char what);
		void remove_momentum_contribution(char value);

		friend void normalOrder(std::vector<Term>& terms);
		friend void commutator(std::vector<Term>& reciever, const Term& left, const Term& right);
		friend std::ostream& operator<<(std::ostream& os, const Term& term);
	};

	void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const std::vector<Term>& right);
	inline void commutator(std::vector<Term>& reciever, const Term& left, const std::vector<Term>& right) {
		const std::vector<Term> buffer = { left };
		commutator(reciever, buffer, right);
	};
	inline void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const Term& right) {
		const std::vector<Term> buffer = { right };
		commutator(reciever, left, buffer);
	};

	inline bool operator==(const Term& lhs, const Term& rhs) {
		return lhs.is_equal(rhs);
	};
	inline bool operator!=(const Term& lhs, const Term& rhs) {
		return !(lhs.is_equal(rhs));
	};

	std::ostream& operator<<(std::ostream& os, const Coefficient& coeff);
	std::ostream& operator<<(std::ostream& os, const std::vector<Coefficient>& coeffs);
	std::ostream& operator<<(std::ostream& os, const std::vector<Term>& terms);

	void clear_duplicates(std::vector<Term>& terms);
	void cleanUp(std::vector<Term>& terms);
	inline void hermitianConjugate(std::vector<Term>& terms) {
		for (auto& t : terms) {
			t.hermitianConjugate();
		}
	};
	inline void rename_momenta(std::vector<Term>& terms, char what, char to) {
		for (auto& t : terms) {
			t.rename_momenta(what, to);
		}
	};
	inline std::string toStringWithoutPrefactor(const std::vector<Term>& terms) {
		std::string ret = "";
		for (size_t i = 0; i < terms.size(); i++)
		{
			if (terms[i].multiplicity < 0) {
				ret += "-";
			}
			else if (i > 0) {
				ret += "+";
			}
			ret += terms[i].toStringWithoutPrefactor();
		}
		return ret;
	}
}