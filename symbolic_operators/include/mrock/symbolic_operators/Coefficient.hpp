/**
 * @file Coefficient.hpp
 * @brief Defines the Coefficient structure used in symbolic operators.
 */

#pragma once

#include "Operator.hpp"
#include "IndexWrapper.hpp"
#include "MomentumList.hpp"
#include <optional>
#include <functional>

namespace mrock::symbolic_operators {

	/**
	 * @struct Coefficient
	 * @brief Represents a coefficient. Various symmetries are pre defined (e.g. inversion symmetry) and can be toggled on or off
	 * 		A custom symmetry can also be provided
	 */
	struct Coefficient {
		std::string name; ///< Name of the coefficient.
		MomentumList momenta; ///< List of momenta associated with the coefficient.
		IndexWrapper indizes; ///< Contains all indices, standard: first index = spin, all others arbitrary, e.g., orbitals, bands, etc.
		std::optional<std::function<void(Coefficient&)>> custom_symmetry = std::nullopt; ///< Optional custom symmetry function.
		bool inversion_symmetry{ true }; ///< Indicates if V(k) = V(-k).
		bool is_symmetrized_interaction{ }; ///< Indicates if the interaction is symmetrized, i.e., V(k, k', q) = V(k', k, -q)
		bool Q_changes_sign{}; ///< Indicates if V(k+Q) = -V(k).
		bool is_real{ true }; ///< Indicates if V^* = V. Default is true.
		bool is_daggered{}; ///< Indicates if the coefficient is daggered.

		/**
		 * @brief Serializes the Coefficient object.
		 * @tparam Archive The archive type.
		 * @param ar The archive object.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& name;
			ar& momenta;
			ar& indizes;
			ar& is_daggered;
		}

		/**
		 * @brief Default constructor.
		 */
		Coefficient() = default;

		/**
		 * @brief Constructs a Coefficient with a given name
		 * 
		 * @param _name The name of the coefficient
		 */
		explicit Coefficient(const std::string& _name);

		/**
		 * @brief Constructs a Coefficient with a given name, a single momentum, and a set of indizes (can be of size 1)
		 * 
		 * @param _name The name of the coefficient
		 * @param _momentum The momentum of the coefficient
		 * @param _indizes The indizes of the coefficient
		 * @param _Q_changes_sign Toggles the V(k+Q) = -V(k) symmetry. Default is false.
		 * @param _inversion_symmetry Toggles inversion symmetry, V(k) = V(-k). Default is true.
		 * @param _is_daggered Toggles whether the coefficient is a complex conjugate or not. Default is false.
		 */
		Coefficient(const std::string& _name, const Momentum& _momentum, const IndexWrapper& _indizes, bool _Q_changes_sign = false, bool _inversion_symmetry = true, bool _is_daggered = false);

		/**
		 * @brief Constructs a Coefficient with a given name, and a single momentum
		 * 
		 * @param _name The name of the coefficient
		 * @param _momentum The momentum of the coefficient
		 * @param _Q_changes_sign Toggles the V(k+Q) = -V(k) symmetry. Default is false.
		 * @param _inversion_symmetry Toggles inversion symmetry, V(k) = V(-k). Default is true.
		 * @param _is_daggered Toggles whether the coefficient is a complex conjugate or not. Default is false.
		 */
		Coefficient(const std::string& _name, const Momentum& _momentum, bool _Q_changes_sign = false, bool _inversion_symmetry = true, bool _is_daggered = false);

		/**
		 * @brief Constructs a Coefficient with a given name, multiple momenta, and a set of indizes (can be of size 1)
		 * 
		 * @param _name The name of the coefficient
		 * @param _momenta The momenta of the coefficient, in order. Usually occurs for interactions, i.e., V(k, k', q)
		 * @param _indizes The indizes of the coefficient
		 * @param _Q_changes_sign Toggles the V(k+Q) = -V(k) symmetry. Default is false.
		 * @param _inversion_symmetry Toggles inversion symmetry, V(k) = V(-k). Default is true.
		 * @param _is_daggered Toggles whether the coefficient is a complex conjugate or not. Default is false.
		 */
		Coefficient(const std::string& _name, const MomentumList& _momenta, const IndexWrapper& _indizes = IndexWrapper{}, bool _Q_changes_sign = false, bool _inversion_symmetry = true, bool _is_daggered = false);

		/**
		 * @brief Generates a real and inversion symmetric coefficient.
		 * @param _name The name of the coefficient.
		 * @param _momenta The list of momenta.
		 * @param _custom_symmetry Optional custom symmetry function.
		 * @return A real and inversion symmetric Coefficient.
		 */
		static Coefficient RealInversionSymmetric(const std::string& name, const MomentumList& momenta, const std::optional<std::function<void(Coefficient&)>>& custom_symmetry = std::nullopt);

		/**
		 * @brief Generates a real Coefficient with V(k, k', q) = V(k', k, -q).
		 * @param name The name of the coefficient.
		 * @param momenta The list of momenta.
		 * @param custom_symmetry Optional custom symmetry function.
		 * @return A real interaction Coefficient.
		 */
		static Coefficient RealInteraction(const std::string& name, const MomentumList& momenta, const std::optional<std::function<void(Coefficient&)>>& custom_symmetry = std::nullopt);

		/**
		 * @brief Generates a Coefficient as they occur on a honeycomb lattice.
		 * @param name The name of the coefficient.
		 * @param momentum The momentum.
		 * @param daggered Indicates if the coefficient is daggered.
		 * @param is_real Indicates if the coefficient is real. Default is true.
		 * @param custom_symmetry Optional custom symmetry function.
		 * @return A Coefficient with the symmetries of a honeycomb lattice.
		 */
		static Coefficient HoneyComb(const std::string& name, const Momentum& momentum, bool daggered, bool is_real = true, const std::optional<std::function<void(Coefficient&)>>& custom_symmetry = std::nullopt);

		/**
		 * @brief Generates a Coefficient that does not depend on any momentum
		 * @param name The name of the coefficient.
		 * @param indizes The indizes of the coefficient. Default is no index.
		 * @param daggered Indicates if the coefficient is daggered. Default is false.
		 * @param is_real Indicates if the coefficient is real. Default is true.
		 * @return A Coefficient with the symmetries of a honeycomb lattice.
		 */
		static Coefficient Constant(const std::string& name, const IndexWrapper& indizes = IndexWrapper{}, bool is_real = true, bool daggered = false);

		/**
		 * @brief Parses a string to create a Coefficient.
		 * @param expression The string expression.
		 * @param _Q_changes_sign Indicates if Q changes sign.
		 * @param _inversion_symmetry Indicates if inversion symmetry is present.
		 * @return A parsed Coefficient.
		 */
		static Coefficient parse_string(const std::string& expression, bool _Q_changes_sign = false, bool _inversion_symmetry = true);

		/**
		 * @brief Parses a string to create a standard interaction Coefficient.
		 * @param expression The string expression.
		 * @return A parsed interaction Coefficient.
		 */
		static Coefficient parse_interaction_string(const std::string& expression);

		/**
		 * @brief Checks if the coefficient uses a specific index.
		 * @param index The index to check.
		 * @return True if the index is used, false otherwise.
		 */
		inline bool uses_index(const Index index) const noexcept;

		/**
		 * @brief Checks if the coefficient depends on momentum.
		 * @return True if it depends on momentum, false otherwise.
		 */
		inline bool depends_on_momentum() const noexcept;

		/**
		 * @brief Checks if the coefficient depends on a specific momentum.
		 * @param momentum The momentum to check.
		 * @return True if it depends on the momentum, false otherwise.
		 */
		inline bool depends_on(const MomentumSymbol::name_type momentum) const noexcept;

		/**
		 * @brief Checks if the coefficient depends on two momenta, e.g, k-l.
		 * 		The Current implementation is restricted to a MomentumList of size 1, i.e., V(k) or V(k-l) but not V(k, l)
		 * @return True if it depends on two momenta, false otherwise.
		 */
		inline bool depends_on_two_momenta() const noexcept;

		/**
		 * @brief Toggles the daggered state of the operator.
		 * @return A reference to *this
		 */
		inline Coefficient& hermitian_conjugate_inplace();

		/**
		 * @brief Creates hermitian conjugate of this as a new object.
		 * @return Returns the new object.
		 */
		inline Coefficient hermitian_conjugate() const;

		/**
		 * @brief Inverts the momentum of the coefficient.
		 * @param what The momentum to invert.
		 */
		void invert_momentum(const MomentumSymbol::name_type what);

		/**
		 * @brief Utilizes V(k, k', q) = V(k', k, -q) symmetry.
		 */
		void use_symmetric_interaction_exchange();

		/**
		 * @brief Utilizes V(k, k', q) = V(-k, -k', -q) symmetry.
		 */
		void use_symmetric_interaction_inversion();

		/**
		 * @brief Removes a momentum contribution from the coefficient.
		 * @param value The momentum value to remove.
		 */
		void remove_momentum_contribution(const MomentumSymbol::name_type value);

		/**
		 * @brief Applies the custom symmetry function if it exists.
		 */
		void apply_custom_symmetry();
	};

	/**
	 * @brief Equality operator for Coefficient.
	 * @param lhs The left-hand side Coefficient.
	 * @param rhs The right-hand side Coefficient.
	 * @return True if the coefficients are equal, false otherwise.
	 */
	inline bool operator==(const Coefficient& lhs, const Coefficient& rhs) {
		if (lhs.name != rhs.name) return false;
		if (lhs.momenta != rhs.momenta) return false;
		if (lhs.is_daggered != rhs.is_daggered) return false;
		return (lhs.indizes == rhs.indizes);
	}

	/**
	 * @brief Inequality operator for Coefficient.
	 * @param lhs The left-hand side Coefficient.
	 * @param rhs The right-hand side Coefficient.
	 * @return True if the coefficients are not equal, false otherwise.
	 */
	inline bool operator!=(const Coefficient& lhs, const Coefficient& rhs) { return !(lhs == rhs); }

	// Inline definitions
	bool Coefficient::uses_index(const Index index) const noexcept {
		for (const auto& idx : indizes) {
			if (idx == index) return true;
		}
		return false;
	}
	bool Coefficient::depends_on_momentum() const noexcept {
		if (this->momenta.empty()) return false;
		return std::any_of(this->momenta.begin(), this->momenta.end(), [](const Momentum& momentum) {
			return !momentum.momentum_list.empty();
			});
	}
	bool Coefficient::depends_on(const MomentumSymbol::name_type momentum) const noexcept {
		if (this->momenta.empty()) return false;
		return std::any_of(this->momenta.begin(), this->momenta.end(), [momentum](const Momentum& mom) {
			return mom.uses(momentum);
			});
	}
	// This function determines whether the coefficient depends on something like k-l
	// Currently, this only makes sense if the coefficient does not depend on
	bool Coefficient::depends_on_two_momenta() const noexcept {
		assert(momenta.size() == 1U);
		return this->momenta.front().momentum_list.size() == 2U;
	}
	Coefficient& Coefficient::hermitian_conjugate_inplace() {
		if (is_real) return *this;
		is_daggered = !is_daggered;
		return *this;
	}
	Coefficient Coefficient::hermitian_conjugate() const {
		Coefficient copy(*this);
		copy.hermitian_conjugate_inplace();
		return copy;
	}
} // namespace mrock::symbolic_operators