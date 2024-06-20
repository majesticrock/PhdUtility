#pragma once
#include "Term.hpp"
#include "WickOperator.hpp"
#include "WickOperatorTemplate.hpp"
#include <algorithm>
#include <Utility/Fractional.hpp>

namespace SymbolicOperators {
	class Term;
	struct WickTerm
	{
	private:
		void string_parser(std::string&& expression);
	public:
		IntFractional multiplicity{};
		std::vector<Coefficient> coefficients;
		SumContainer sums;
		std::vector<WickOperator> operators;

		// symbolises the Kronecker delta
		std::vector<KroneckerDelta<Momentum>> delta_momenta;
		std::vector<KroneckerDelta<Index>> delta_indizes;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& multiplicity;
			ar& coefficients;
			ar& sums;
			ar& operators;
			ar& delta_momenta;
			ar& delta_indizes;
		}

		std::vector<Operator> temporary_operators;

		explicit WickTerm(const Term* base);
		explicit WickTerm(const Term& base);
		WickTerm() = default;
		WickTerm(const WickTerm& base, const TemplateResult::SingleResult& result);

		explicit WickTerm(const std::string& expression);

		inline bool includesType(const OperatorType operator_type) const {
			return std::any_of(this->operators.begin(), this->operators.end(),
				[operator_type](const WickOperator& op) { return op.type == operator_type; });
		};
		inline bool hasSingleCoefficient() const noexcept {
			return this->coefficients.size() == 1U;
		};
		inline bool usesIndex(const Index index) const noexcept {
			for (const auto& op : operators) {
				if (op.usesIndex(index)) return true;
			}
			for (const auto& coeff : coefficients) {
				if (coeff.usesIndex(index)) return true;
			}
			return false;
		}
		inline bool isIdentity() const noexcept {
			return this->operators.empty();
		}
		inline bool isBilinear() const noexcept {
			return this->operators.size() == 1U;
		};
		inline bool isQuartic() const noexcept {
			return this->operators.size() == 2U;
		}
		// Returns this->multiplicity, but always as a double
		inline double getFactor() const noexcept {
			return static_cast<double>(this->multiplicity);
		}
		// Returns the position of the first operator that depends on 'momentum'
		// Returns -1 if no operators depend on 'momentum'.
		inline int whichOperatorDependsOn(char momentum) const noexcept {
			for (int i = 0U; i < operators.size(); ++i)
			{
				if (operators[i].dependsOn(momentum)) return i;
			}
			return -1;
		}
		inline const Coefficient& getFirstCoefficient() const {
			assert(!(this->coefficients.empty()));
			return this->coefficients.front();
		};

		inline bool handled() const noexcept {
			if (this->temporary_operators.empty()) return true;
			return !(this->operators.empty());
		}

		// returns false if there is atleast one delta
		// or a combination of deltas, that can never be achieved
		// for example delta_k,k+Q, as k can never be equal to k+Q
		bool setDeltas();
		// May call setDeltas. If setDeltas returns false this functions also returns false
		// In all other cases it returns true
		bool computeSums();
		void discardZeroMomenta();
		void renameSums();
		void sort();
		void includeTemplateResult(const TemplateResult::SingleResult& result);
	};

	inline bool operator==(const WickOperator& lhs, const WickOperator& rhs) {
		if (lhs.type != rhs.type) return false;
		if (lhs.isDaggered != rhs.isDaggered) return false;
		if (lhs.momentum != rhs.momentum) return false;
		return (lhs.indizes == rhs.indizes);
	};
	inline bool operator!=(const WickOperator& lhs, const WickOperator& rhs) {
		return !(lhs == rhs);
	};
	inline bool operator==(const WickTerm& lhs, const WickTerm& rhs) {
		if (lhs.coefficients != rhs.coefficients) return false;
		if (lhs.sums != rhs.sums) return false;
		if (lhs.delta_indizes != rhs.delta_indizes) return false;
		if (lhs.delta_momenta != rhs.delta_momenta) return false;
		if (lhs.operators != rhs.operators) return false;
		return true;
	};
	inline bool operator!=(const WickTerm& lhs, const WickTerm& rhs) {
		return !(lhs == rhs);
	};

	struct WickTermCollector : public Utility::VectorWrapper<WickTerm> {
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& _vector;
		};
	};

	WickTermCollector& operator+=(WickTermCollector& lhs, const WickTerm& rhs);
	WickTermCollector& operator-=(WickTermCollector& lhs, const WickTerm& rhs);
	WickTermCollector& operator+=(WickTermCollector& lhs, const WickTermCollector& rhs);
	WickTermCollector& operator-=(WickTermCollector& lhs, const WickTermCollector& rhs);
	inline WickTermCollector operator+(WickTermCollector lhs, const WickTerm& rhs) { lhs += rhs; return lhs; };
	inline WickTermCollector operator-(WickTermCollector lhs, const WickTerm& rhs) { lhs -= rhs; return lhs; };
	inline WickTermCollector operator+(const WickTerm& lhs, WickTermCollector rhs) { rhs += lhs; return rhs; };
	inline WickTermCollector operator-(const WickTerm& lhs, WickTermCollector rhs) { rhs -= lhs; return rhs; };
	inline WickTermCollector operator+(WickTermCollector lhs, const WickTermCollector& rhs) { lhs += rhs; return lhs; };
	inline WickTermCollector operator-(WickTermCollector lhs, const WickTermCollector& rhs) { lhs -= rhs; return lhs; };

	std::ostream& operator<<(std::ostream& os, const WickTerm& term);
	std::ostream& operator<<(std::ostream& os, const WickTermCollector& terms);
}