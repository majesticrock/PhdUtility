#pragma once
#include "KroneckerDelta.hpp"
#include "Coefficient.hpp"
#include "SymbolicSum.hpp"
#include <Utility/Fractional.hpp>

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

		friend struct WickTerm;
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
		void print() const;
		inline void flipSign() {
			this->multiplicity *= -1;
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
		inline bool isEqual(const Term& other) const {
			if (this->coefficients != other.coefficients) return false;
			if (this->sums != other.sums) return false;
			if (this->delta_indizes != other.delta_indizes) return false;
			if (this->delta_momenta != other.delta_momenta) return false;
			if (this->operators != other.operators) return false;

			return true;
		};

		std::string toStringWithoutPrefactor() const;

		friend void normalOrder(std::vector<Term>& terms);
		friend void commutator(std::vector<Term>& reciever, const Term& left, const Term& right);
		friend std::ostream& operator<<(std::ostream& os, const Term& term);

		inline void hermitianConjugate() {
			std::reverse(this->operators.begin(), this->operators.end());
			for (auto& op : this->operators) {
				op.hermitianConjugate();
			}
		};
		inline void renameMomenta(char what, char to) {
			for (auto& op : operators)
			{
				op.momentum.replaceOccurances(what, Momentum(to));
			}
		};
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
		return lhs.isEqual(rhs);
	};
	inline bool operator!=(const Term& lhs, const Term& rhs) {
		return !(lhs.isEqual(rhs));
	};

	std::ostream& operator<<(std::ostream& os, const Coefficient& coeff);
	std::ostream& operator<<(std::ostream& os, const std::vector<Coefficient>& coeffs);
	std::ostream& operator<<(std::ostream& os, const std::vector<Term>& terms);

	void cleanUp(std::vector<Term>& terms);
	inline void hermitianConjugate(std::vector<Term>& terms) {
		for (auto& t : terms) {
			t.hermitianConjugate();
		}
	};
	inline void renameMomenta(std::vector<Term>& terms, char what, char to) {
		for (auto& t : terms) {
			t.renameMomenta(what, to);
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