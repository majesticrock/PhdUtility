/**
 * @file WickTerm.hpp
 * @brief Defines the WickTerm structure and related functions.
 */

#pragma once
#include "Term.hpp"
#include "WickOperator.hpp"
#include "WickOperatorTemplate.hpp"
#include <algorithm>
#include <mrock/utility/Fractional.hpp>

namespace mrock::symbolic_operators {

	/**
	 * @class WickTerm
	 * @brief A structure representing a Wick term.
	 */
	struct WickTerm {
	private:
		/**
		 * @brief Parses a string expression to initialize the WickTerm.
		 * 
		 * @param expression The string expression.
		 */
		void string_parser(std::string&& expression);

	public:
		IntFractional multiplicity{}; ///< The multiplicity of the term.
		std::vector<Coefficient> coefficients; ///< The coefficients of the term.
		SumContainer sums; ///< The sums in the term.
		std::vector<WickOperator> operators; ///< The operators in the term.

		// symbolises the Kronecker delta
		std::vector<KroneckerDelta<Momentum>> delta_momenta; ///< The momentum deltas.
		std::vector<KroneckerDelta<Index>> delta_indizes; ///< The index deltas.

		/**
		 * @brief Serializes the WickTerm object.
		 * 
		 * @tparam Archive The type of the archive.
		 * @param ar The archive object.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& multiplicity;
			ar& coefficients;
			ar& sums;
			ar& operators;
			ar& delta_momenta;
			ar& delta_indizes;
		}

		std::vector<Operator> temporary_operators; ///< Temporary operators used in the term.

		/**
		 * @brief Constructs a WickTerm object from a base Term pointer.
		 * 
		 * @param base The base Term pointer.
		 */
		explicit WickTerm(const Term* base);

		/**
		 * @brief Constructs a WickTerm object from a base Term reference.
		 * 
		 * @param base The base Term reference.
		 */
		explicit WickTerm(const Term& base);

		/**
		 * @brief Default constructor for WickTerm.
		 */
		WickTerm() = default;

		/**
		 * @brief Constructs a WickTerm object from a base WickTerm and a TemplateResult::SingleResult.
		 * 
		 * @param base The base WickTerm.
		 * @param result The TemplateResult::SingleResult.
		 */
		WickTerm(const WickTerm& base, const TemplateResult::SingleResult& result);

		/**
		 * @brief Constructs a WickTerm object from a string expression.
		 * 
		 * @param expression The string expression.
		 */
		explicit WickTerm(const std::string& expression);

		/**
		 * @brief Checks if the term includes a specific operator type.
		 * 
		 * @param operator_type The operator type to check.
		 * @return true if the term includes the operator type.
		 * @return false otherwise.
		 */
		inline bool includes_type(const OperatorType operator_type) const;

		/**
		 * @brief Checks if the term has a single coefficient.
		 * 
		 * @return true if the term has a single coefficient.
		 * @return false otherwise.
		 */
		inline bool has_single_coefficient() const noexcept;

		/**
		 * @brief Checks if the term uses a specific index.
		 * 
		 * @param index The index to check.
		 * @return true if the term uses the index.
		 * @return false otherwise.
		 */
		inline bool uses_index(const Index index) const noexcept;

		/**
		 * @brief Checks if the term is an identity term.
		 * 
		 * @return true if the term is an identity term.
		 * @return false otherwise.
		 */
		inline bool is_identity() const noexcept;

		/**
		 * @brief Checks if the term is bilinear.
		 * 
		 * @return true if the term is bilinear.
		 * @return false otherwise.
		 */
		inline bool is_bilinear() const noexcept;

		/**
		 * @brief Checks if the term is quartic.
		 * 
		 * @return true if the term is quartic.
		 * @return false otherwise.
		 */
		inline bool is_quartic() const noexcept;

		/**
		 * @brief Returns the multiplicity as a double.
		 * 
		 * @return double The multiplicity as a double.
		 */
		inline double get_factor() const noexcept;

		/**
		 * @brief Returns the position of the first operator that depends on a specific momentum.
		 * 
		 * @param momentum The momentum to check.
		 * @return int The position of the first operator that depends on the momentum, or -1 if none.
		 */
		inline int which_operator_depends_on(const MomentumSymbol::name_type momentum) const noexcept;

		/**
		 * @brief Returns the first coefficient in the term.
		 * 
		 * @return const Coefficient& The first coefficient.
		 */
		inline const Coefficient& get_first_coefficient() const;

		/**
		 * @brief Checks if the term has been handled.
		 * 
		 * @return true if the term has been handled.
		 * @return false otherwise.
		 */
		inline bool handled() const noexcept;

		/**
		 * @brief Sets the deltas in the term.
		 * 
		 * @return true if the deltas were set successfully.
		 * @return false otherwise.
		 */
		bool set_deltas();

		/**
		 * @brief Computes the sums in the term.
		 * 
		 * @return true if the sums were computed successfully.
		 * @return false otherwise.
		 */
		bool compute_sums();

		/**
		 * @brief Discards zero momenta in the term.
		 */
		void discard_zero_momenta();

		/**
		 * @brief Renames the sums in the term.
		 */
		void rename_sums();

		/**
		 * @brief Sorts the elements in the term.
		 */
		void sort();

		/**
		 * @brief Includes a template result in the term.
		 * 
		 * @param result The TemplateResult::SingleResult to include.
		 */
		void include_template_result(const TemplateResult::SingleResult& result);

		/**
		 * @brief Inverts a momentum in the term.
		 * 
		 * @param what The momentum to invert.
		 */
		void invert_momentum(const MomentumSymbol::name_type what);

		/**
		 * @brief Inverts a momentum sum in the term.
		 * 
		 * @param what The momentum sum to invert.
		 */
		void invert_momentum_sum(const MomentumSymbol::name_type what);

		/**
		 * @brief Removes a momentum contribution from the term.
		 * 
		 * @param value The momentum value to remove.
		 */
		inline void remove_momentum_contribution(const MomentumSymbol::name_type value);
	};

	/**
	 * @brief Equality operator for WickOperator.
	 * 
	 * @param lhs The left-hand side WickOperator.
	 * @param rhs The right-hand side WickOperator.
	 * @return true if the two WickOperator objects are equal.
	 * @return false otherwise.
	 */
	inline bool operator==(const WickOperator& lhs, const WickOperator& rhs);

	/**
	 * @brief Inequality operator for WickOperator.
	 * 
	 * @param lhs The left-hand side WickOperator.
	 * @param rhs The right-hand side WickOperator.
	 * @return true if the two WickOperator objects are not equal.
	 * @return false otherwise.
	 */
	inline bool operator!=(const WickOperator& lhs, const WickOperator& rhs);

	/**
	 * @brief Equality operator for WickTerm.
	 * 
	 * @param lhs The left-hand side WickTerm.
	 * @param rhs The right-hand side WickTerm.
	 * @return true if the two WickTerm objects are equal.
	 * @return false otherwise.
	 */
	inline bool operator==(const WickTerm& lhs, const WickTerm& rhs);

	/**
	 * @brief Inequality operator for WickTerm.
	 * 
	 * @param lhs The left-hand side WickTerm.
	 * @param rhs The right-hand side WickTerm.
	 * @return true if the two WickTerm objects are not equal.
	 * @return false otherwise.
	 */
	inline bool operator!=(const WickTerm& lhs, const WickTerm& rhs);

	/**
	 * @class WickTermCollector
	 * @brief A wrapper for a vector of WickTerm objects.
	 */
	struct WickTermCollector : public mrock::utility::VectorWrapper<WickTerm> {
		/**
		 * @brief Serializes the WickTermCollector object.
		 * 
		 * @tparam Archive The type of the archive.
		 * @param ar The archive object.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& _vector;
		};
	};

	/**
	 * @brief Addition assignment operator for WickTermCollector and WickTerm.
	 * 
	 * @param lhs The left-hand side WickTermCollector.
	 * @param rhs The right-hand side WickTerm.
	 * @return WickTermCollector& The updated WickTermCollector.
	 */
	WickTermCollector& operator+=(WickTermCollector& lhs, const WickTerm& rhs);

	/**
	 * @brief Subtraction assignment operator for WickTermCollector and WickTerm.
	 * 
	 * @param lhs The left-hand side WickTermCollector.
	 * @param rhs The right-hand side WickTerm.
	 * @return WickTermCollector& The updated WickTermCollector.
	 */
	WickTermCollector& operator-=(WickTermCollector& lhs, const WickTerm& rhs);

	/**
	 * @brief Addition assignment operator for two WickTermCollector objects.
	 * 
	 * @param lhs The left-hand side WickTermCollector.
	 * @param rhs The right-hand side WickTermCollector.
	 * @return WickTermCollector& The updated WickTermCollector.
	 */
	WickTermCollector& operator+=(WickTermCollector& lhs, const WickTermCollector& rhs);

	/**
	 * @brief Subtraction assignment operator for two WickTermCollector objects.
	 * 
	 * @param lhs The left-hand side WickTermCollector.
	 * @param rhs The right-hand side WickTermCollector.
	 * @return WickTermCollector& The updated WickTermCollector.
	 */
	WickTermCollector& operator-=(WickTermCollector& lhs, const WickTermCollector& rhs);

	/**
	 * @brief Addition operator for WickTermCollector and WickTerm.
	 * 
	 * @param lhs The left-hand side WickTermCollector.
	 * @param rhs The right-hand side WickTerm.
	 * @return WickTermCollector The resulting WickTermCollector.
	 */
	inline WickTermCollector operator+(WickTermCollector lhs, const WickTerm& rhs) { lhs += rhs; return lhs; };

	/**
	 * @brief Subtraction operator for WickTermCollector and WickTerm.
	 * 
	 * @param lhs The left-hand side WickTermCollector.
	 * @param rhs The right-hand side WickTerm.
	 * @return WickTermCollector The resulting WickTermCollector.
	 */
	inline WickTermCollector operator-(WickTermCollector lhs, const WickTerm& rhs) { lhs -= rhs; return lhs; };

	/**
	 * @brief Addition operator for WickTerm and WickTermCollector.
	 * 
	 * @param lhs The left-hand side WickTerm.
	 * @param rhs The right-hand side WickTermCollector.
	 * @return WickTermCollector The resulting WickTermCollector.
	 */
	inline WickTermCollector operator+(const WickTerm& lhs, WickTermCollector rhs) { rhs += lhs; return rhs; };

	/**
	 * @brief Subtraction operator for WickTerm and WickTermCollector.
	 * 
	 * @param lhs The left-hand side WickTerm.
	 * @param rhs The right-hand side WickTermCollector.
	 * @return WickTermCollector The resulting WickTermCollector.
	 */
	inline WickTermCollector operator-(const WickTerm& lhs, WickTermCollector rhs) { rhs -= lhs; return rhs; };

	/**
	 * @brief Addition operator for two WickTermCollector objects.
	 * 
	 * @param lhs The left-hand side WickTermCollector.
	 * @param rhs The right-hand side WickTermCollector.
	 * @return WickTermCollector The resulting WickTermCollector.
	 */
	inline WickTermCollector operator+(WickTermCollector lhs, const WickTermCollector& rhs) { lhs += rhs; return lhs; };

	/**
	 * @brief Subtraction operator for two WickTermCollector objects.
	 * 
	 * @param lhs The left-hand side WickTermCollector.
	 * @param rhs The right-hand side WickTermCollector.
	 * @return WickTermCollector The resulting WickTermCollector.
	 */
	inline WickTermCollector operator-(WickTermCollector lhs, const WickTermCollector& rhs) { lhs -= rhs; return lhs; };

	/**
	 * @brief Stream insertion operator for WickTerm.
	 * 
	 * @param os The output stream.
	 * @param term The WickTerm object.
	 * @return std::ostream& The updated output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const WickTerm& term);

	/**
	 * @brief Stream insertion operator for WickTermCollector.
	 * 
	 * @param os The output stream.
	 * @param terms The WickTermCollector object.
	 * @return std::ostream& The updated output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const WickTermCollector& terms);

	/**
	 * @class bad_term_exception
	 * @brief An exception class for bad terms.
	 */
	class bad_term_exception : public std::runtime_error {
	protected:
		const WickTerm _term; ///< The bad term.

	public:
		/**
		 * @brief Constructs a bad_term_exception object.
		 * 
		 * @param what_arg The error message.
		 * @param term The bad term.
		 */
		bad_term_exception(const std::string& what_arg, const WickTerm& term) : std::runtime_error(what_arg), _term(term) {};

		/**
		 * @brief Constructs a bad_term_exception object.
		 * 
		 * @param what_arg The error message.
		 * @param term The bad term.
		 */
		bad_term_exception(const char* what_arg, const WickTerm& term) : std::runtime_error(what_arg), _term(term) {};

		/**
		 * @brief Returns the bad term.
		 * 
		 * @return const WickTerm& The bad term.
		 */
		const WickTerm& which_term() const noexcept { return this->_term; };
	};

	// Inline definitions
	inline bool WickTerm::includes_type(const OperatorType operator_type) const {
		return std::any_of(this->operators.begin(), this->operators.end(),
			[operator_type](const WickOperator& op) { return op.type == operator_type; });
	}
	inline bool WickTerm::has_single_coefficient() const noexcept {
		return this->coefficients.size() == 1U;
	}
	inline bool WickTerm::uses_index(const Index index) const noexcept {
		for (const auto& op : operators) {
			if (op.uses_index(index)) return true;
		}
		for (const auto& coeff : coefficients) {
			if (coeff.uses_index(index)) return true;
		}
		return false;
	}
	inline bool WickTerm::is_identity() const noexcept {
		return this->operators.empty();
	}
	inline bool WickTerm::is_bilinear() const noexcept {
		return this->operators.size() == 1U;
	}
	inline bool WickTerm::is_quartic() const noexcept {
		return this->operators.size() == 2U;
	}
	inline double WickTerm::get_factor() const noexcept {
		return static_cast<double>(this->multiplicity);
	}
	inline int WickTerm::which_operator_depends_on(const MomentumSymbol::name_type momentum) const noexcept {
		for (int i = 0U; i < operators.size(); ++i)
		{
			if (operators[i].depends_on(momentum)) return i;
		}
		return -1;
	}
	inline const Coefficient& WickTerm::get_first_coefficient() const {
		assert(!(this->coefficients.empty()));
		return this->coefficients.front();
	}
	inline bool WickTerm::handled() const noexcept {
		if (this->temporary_operators.empty()) return true;
		return !(this->operators.empty());
	}
	inline void WickTerm::remove_momentum_contribution(const MomentumSymbol::name_type value) {
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
		std::erase_if(sums.momenta.summations, [&](const MomentumSymbol::name_type sum_idx) { return sum_idx == value; });
	}
	inline bool operator==(const WickOperator& lhs, const WickOperator& rhs) {
		if (lhs.type != rhs.type) return false;
		if (lhs.is_daggered != rhs.is_daggered) return false;
		if (lhs.momentum != rhs.momentum) return false;
		return (lhs.indizes == rhs.indizes);
	}
	inline bool operator!=(const WickOperator& lhs, const WickOperator& rhs) {
		return !(lhs == rhs);
	}
	inline bool operator==(const WickTerm& lhs, const WickTerm& rhs) {
		if (lhs.coefficients != rhs.coefficients) return false;
		if (lhs.sums != rhs.sums) return false;
		if (lhs.delta_indizes != rhs.delta_indizes) return false;
		if (lhs.delta_momenta != rhs.delta_momenta) return false;
		if (lhs.operators != rhs.operators) return false;
		return true;
	}
	inline bool operator!=(const WickTerm& lhs, const WickTerm& rhs) {
		return !(lhs == rhs);
	}
} // namespace mrock::symbolic_operators