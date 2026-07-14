#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICKTERM_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICKTERM_HPP
/**
 * @file WickTerm.hpp
 * @brief Defines the WickTerm class and related functions.
 */

#include "Term.hpp"
#include "WickOperator.hpp"
#include "WickOperatorTemplate.hpp"
#include "detail/vector_macro.hpp"
#include "OperatorType.hpp"
#include "AbstractTerm.hpp"
#include "Coefficient.hpp"
#include "IndexWrapper.hpp"
#include "KroneckerDelta.hpp"
#include "Momentum.hpp"
#include "MomentumSymbol.hpp"
#include "Operator.hpp"
#include "SumContainer.hpp"

#include <algorithm>
#include <cassert>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mrock::symbolic_operators {
	/**
	 * @class WickTerm
	 * @brief A class representing a term consisting of expectation values, represented via WickOperator objects.
	 * 
	 * Prerequisite: The terms you want to apply Wick's theorem on are saved in an \c std::vector<Term>.
	 * 
	 * Applying Wick's theorem often involves omiting certain expecation values because you know them to be 0 for symmtry reasons. 
	 * Therefore the class \c WickOperatorTemplate exists. Here, you specify, which kind of expectation values will be finite.
	 * In the following, the meaning of the different attributes is listed:
	 * 
	 * \b std::vector< \b IndexComparison \b > \b indexComparison} \n
	 * If \c any_identical is \c true, any two identical indizes are considered valid. 
	 * An example would be in the number operator \f$c_{k, \sigma}^\dagger c_{k, \sigma'}*\f$: 
	 * No matter what \f$ \sigma\f$ is, as long as \f$sigma = sigma'\f$ the expectation value will be finite.
	 * If \c any_identical is \c false, the members \c base and \c other become relevant: 
	 * They define what the indizes need to be, e.g., for a pair annihlation operator \f$c_{-k down} c_{k up}\f$
	 * one would set \c base to \f$ \downarrow \f$ and \c other to \f$ \uparrow\f$.
	 * 
	 * Note, once one operator is set as a template, it is not necessary to set its Hermitian conjugate.
	 * 
	 * \b Momentum \b momentum_difference \n
	 * Defines the allowed difference in momentum, e.g., for a number operator, this would be 0.
	 * Note, this also applies to a standard pair creation/annihilation operator, because in total, 
	 * these operators create/annihilate a particle with \f$ -k \f$ and one with \f$ k \f$, resulting in 0 net momentum.
	 * 
	 * \b OperatorType \b type \n
	 * Specifies what kind of WickOperator will be the result, see \c enum \c OperatorType in WickOperator.hpp.
	 * 
	 * \b bool \b is_sc_type \n
	 * Specifies whether the operator is a pair creation/annihilation operator or a standard \f$c^\dagger c\f$ type term.
	 * 
	 * \b Apply \b Wick's \b theorem
	 * Create an instance of \c WickTermCollector. Then simply call 
	 * \code
	 * WickTermCollector wicks;
	 * wicks_theorem(terms, templates, wicks);
	 * clean_wicks(wicks);
	 * \endcode
	 * Similar to how we worked with the Term class and commutators, it is strongly recommended to call clean_wicks() after applying Wick's theorem.
	 * 
	 * clean_wicks() will also make use of polymorphism to apply symmetries to the term, e.g., inversion symmetry.
	 * For details, see \c WickSymmetry.
	 * 
	 * You can print the the result to the console or utilize boost's serialization to load it later (or within another program).
	 * Serialization can be achieved via this code
	 * \code
	 * std::ofstream ofs("path/to/file.txt");
	 * boost::archive::text_oarchive oa(ofs);
	 * oa << wicks;
	 * ofs.close();
	 * \endcode
	 * To later on load the output use
	 * \code
	 * std::ifstream ifs("path/to/file.txt");
	 * boost::archive::text_iarchive ia(ifs);
	 * target.clear();
	 * ia >> target;
	 * ifs.close();
	 * \endcode
	 * 
	 * or if you want to use this code of the iEoM, there is the class \c TermLoader for easy use.
	 * It loads the terms for the matrices M and N and saves same as class members.
	 * 
	 * @sa WickTermCollector, Coefficient, SumContainer, WickOperator, KroneckerDelta, Momentum, Index, clean_wicks(), wicks_theorem(), TermLoader
	 */
	class WickTerm : public AbstractTerm<WickOperator> {
	private:
		/**
		 * @brief Parses a string expression to initialize the WickTerm.
		 * 
		 * @param expression The string expression.
		 */
		void string_parser(std::string&& expression);

	public:
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
		 * @brief Resolves the Kronecker deltas in the term ( calls \c resolve_momentum_deltas() and \c resolve_index_deltas() )
		 * @return True if successful, false otherwise.
		 */
		bool resolve_deltas();

		/**
		 * @brief Renames the sums in the term.
		 */
		void rename_sums();

		/**
		 * @brief Discards zero momenta in the term.
		 */
		void discard_zero_momenta();


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
		 * @brief Checks if the term can be finite.
		 * Works under the assumption that the initial term was normal ordered with respect to the vacuum, e.g., c^+ c^+ c c.
		 * Then, for instance, <n_k> <n_k> cannot be finite because it must have originated from
		 * c_k^+ c_k^+ c_k c_k, which must be 0 due to the Pauli principle.
		 * IMPORTANT: This only works if the original term structure was c^+ c^+ c c, since c_k^+ c_k c_k^+ c_k can be finite!
		 * 
		 * @return true if the term is forbidden and false otherwise.
		 */
		bool is_pauli_forbidden() const;
	}; // WickTerm

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
	struct WickTermCollector {
		/**
		 * @brief Serializes the WickTermCollector object.
		 * 
		 * @tparam Archive The type of the archive.
		 * @param ar The archive object.
		 * @param version The version of the serialization.
		 */
		std::vector<WickTerm> terms; ///< The collected \c WickTerm objects

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& terms;
		};

		MROCK_VECTOR_WRAPPER_FILL_MEMBERS(WickTerm, terms);
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
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_WICKTERM_HPP
