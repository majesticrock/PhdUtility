/**
 * @file WickSymmetry.hpp
 * @brief Defines symmetries for Wick terms.
 */

#pragma once
#include "OperatorType.hpp"
#include "WickTerm.hpp"
#include <type_traits>

namespace mrock::symbolic_operators {

	/**
	 * @class WickSymmetry
	 * @brief An abstract base class for Wick symmetries.
	 * 
	 * There may be some symmetries that simplify your results, e.g., \f$ \langle O^\dagger \rangle = \langle O \rangle \f$.
	 * These symmetries can be implemented by inheriting from the \c WickSymmetry class 
	 * and defining the member function \c virtual \c void \c apply_to(WickTerm& \c term) \c const.
	 * Then create a \c std::vector<std::unique_ptr<WickSymmetry>> \c symmetries
	 * and make use of polymorphism by calling \c clean_wicks(wicks,symmetries).
	 * There are the following predefined symmetry operations:
	 * 
	 * \b SpinSymmetry \n
	 * Changes all spins of the operators in \c term to \f$ \uparrow \f$.
	 * 
	 * \b InversionSymmetry \n
	 * Flips the momenta in such a way, that the first momentum in a term is always positive, i.e., 
	 * \f$ -k+l \f$ is changed to \f$ k-l\f$ while \f$ k-l\f$ would stay unmodified.
	 * 
	 * \b PhaseSymmetry \n
	 * Takes a list of \c OperatorType as template arguments.
	 * Removes any dagger from all operators with a type from the list.
	 * Example:
	 * \c PhaseSymmetry<SC_Type,CDW_Type removes the dagger from \c SC_Type and \c CDW_Type operators.
	 * 
	 * @sa SpinSymmetry, InversionSymmetry, PhaseSymmetry
	 */
	struct WickSymmetry {
		/**
		 * @brief Applies the symmetry to a Wick term.
		 * 
		 * @param term The Wick term to apply the symmetry to.
		 */
		virtual void apply_to(WickTerm& term) const = 0;

		/**
		 * @brief Virtual destructor for WickSymmetry.
		 */
		virtual ~WickSymmetry() = default;
	};

	/**
	 * @class SpinSymmetry
	 * @brief A symmetry where expectation values for spin up and down are the same.
	 */
	struct SpinSymmetry : public WickSymmetry {
		/**
		 * @brief Applies the spin symmetry to a Wick term.
		 * 
		 * @param term The Wick term to apply the symmetry to.
		 */
		void apply_to(WickTerm& term) const override;
	};

	/**
	 * @class InversionSymmetry
	 * @brief A symmetry where expectation values for k and -k are the same.
	 */
	struct InversionSymmetry : public WickSymmetry {
		/**
		 * @brief Applies the inversion symmetry to a Wick term.
		 * 
		 * @param term The Wick term to apply the symmetry to.
		 */
		void apply_to(WickTerm& term) const override;
	};

	/**
	 * @class PhaseSymmetry
	 * @brief A symmetry where <operator^+> = <operator>.
	 * 
	 * @tparam operators The operator types to which the symmetry applies.
	 */
	template<OperatorType... operators>
	struct PhaseSymmetry : public WickSymmetry {
		/**
		 * @brief Applies the phase symmetry to a Wick term.
		 * 
		 * @param term The Wick term to apply the symmetry to.
		 */
		void apply_to(WickTerm& term) const override {
			for (auto& op : term.operators) {
				auto equals_op = [&op](OperatorType comp) {
					return op.type == comp;
				};
				if ((equals_op(operators) || ...)) {
					op.is_daggered = false;
				}
			}
		};
	};

} // namespace mrock::symbolic_operators