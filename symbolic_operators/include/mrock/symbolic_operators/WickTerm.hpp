#pragma once
#include "Term.hpp"
#include "WickOperator.hpp"
#include "WickOperatorTemplate.hpp"
#include <algorithm>
#include <mrock/utility/Fractional.hpp>

namespace mrock::symbolic_operators {
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

		inline bool includes_type(const OperatorType operator_type) const {
			return std::any_of(this->operators.begin(), this->operators.end(),
				[operator_type](const WickOperator& op) { return op.type == operator_type; });
		};
		inline bool has_single_coefficient() const noexcept {
			return this->coefficients.size() == 1U;
		};
		inline bool uses_index(const Index index) const noexcept {
			for (const auto& op : operators) {
				if (op.uses_index(index)) return true;
			}
			for (const auto& coeff : coefficients) {
				if (coeff.uses_index(index)) return true;
			}
			return false;
		}
		inline bool is_identity() const noexcept {
			return this->operators.empty();
		}
		inline bool is_bilinear() const noexcept {
			return this->operators.size() == 1U;
		};
		inline bool is_quartic() const noexcept {
			return this->operators.size() == 2U;
		}
		// Returns this->multiplicity, but always as a double
		inline double get_factor() const noexcept {
			return static_cast<double>(this->multiplicity);
		}
		// Returns the position of the first operator that depends on 'momentum'
		// Returns -1 if no operators depend on 'momentum'.
		inline int which_operator_depends_on(const MomentumSymbol::name_type momentum) const noexcept {
			for (int i = 0U; i < operators.size(); ++i)
			{
				if (operators[i].depends_on(momentum)) return i;
			}
			return -1;
		}
		inline const Coefficient& get_first_coefficient() const {
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
		void discard_zero_momenta();
		void renameSums();
		void sort();
		void include_template_result(const TemplateResult::SingleResult& result);

		// Inverts a momenta, e.g., q -> -q
		void invert_momentum(const MomentumSymbol::name_type what);
		// Same as invert_momentum, but performs a check, whether 'what' is actually being summed over
		void invert_momentum_sum(const MomentumSymbol::name_type what);

		inline void remove_momentum_contribution(const MomentumSymbol::name_type value) {
			for (auto& coeff : coefficients) {
				coeff.remove_momentum_contribution(value);
			}
			for (auto& op : operators) {
				op.remove_momentum_contribution(value);
			}
			for (auto& delta : delta_momenta) {
				delta.first.remove_contribution(value);
				delta.second.remove_contribution(value);
			}
			std::erase_if(sums.momenta._vector, [&](const MomentumSymbol::name_type sum_idx) { return sum_idx == value; });
		}
	};

	inline bool operator==(const WickOperator& lhs, const WickOperator& rhs) {
		if (lhs.type != rhs.type) return false;
		if (lhs.is_daggered != rhs.is_daggered) return false;
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

	struct WickTermCollector : public mrock::utility::VectorWrapper<WickTerm> {
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

	class bad_term_exception : public std::runtime_error {
	protected:
		const WickTerm _term;
	public:
		bad_term_exception(const std::string& what_arg, const WickTerm& term) : std::runtime_error(what_arg), _term(term) {};
		bad_term_exception(const char* what_arg, const WickTerm& term) : std::runtime_error(what_arg), _term(term) {};

		const WickTerm& which_term() const noexcept { return this->_term; };
	};
}