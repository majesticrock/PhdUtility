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