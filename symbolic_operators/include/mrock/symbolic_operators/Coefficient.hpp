#pragma once
#include "Operator.hpp"
#include "IndexWrapper.hpp"
#include "MomentumList.hpp"
#include <optional>
#include <functional>

namespace mrock::symbolic_operators {
	struct Coefficient {
		std::string name;
		MomentumList momenta;
		// Contains all indizes, standard: first index = spin, all others arbitrary, e.g. orbitals, bands etc
		IndexWrapper indizes;
		std::optional<std::function<void(Coefficient&)>> custom_symmetry = std::nullopt;
		// if Coeff(k) = Coeff(-k)
		bool inversion_symmetry{ true };
		// valid for most interactions that depend on 3 momenta
		// allows V(k, k', q) = V(k', k, -q) and V(k, k', q) = V(-k, -k', -q)
		bool is_symmetrized_interaction{ };
		// if Coeff(k+Q) = -Coeff(k)
		bool Q_changes_sign{};
		bool is_daggered{};

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& name;
			ar& momenta;
			ar& indizes;
			ar& is_daggered;
		}

		Coefficient() = default;
		explicit Coefficient(const std::string& _name);
		Coefficient(const std::string& _name, const Momentum& _momentum, const IndexWrapper& _indizes, bool _Q_changes_sign = false, bool _inversion_symmetry = true, bool _is_daggered = false);
		Coefficient(const std::string& _name, const Momentum& _momentum, bool _Q_changes_sign = false, bool _inversion_symmetry = true, bool _is_daggered = false);
		Coefficient(const std::string& _name, const MomentumList& _momenta, const IndexWrapper& _indizes = IndexWrapper{}, bool _Q_changes_sign = false, bool _inversion_symmetry = true, bool _is_daggered = false);

		// Generates a real and inversion symmetric coefficient
		static Coefficient RealInversionSymmetric(const std::string& _name, const MomentumList& _momenta, const std::optional<std::function<void(Coefficient&)>>& _custom_symmetry = std::nullopt);
		// Generates a real Coefficient with V(k, k', q) = V(k', k, -q)
		static Coefficient RealInteraction(const std::string& _name, const MomentumList& _momenta, const std::optional<std::function<void(Coefficient&)>>& _custom_symmetry = std::nullopt);
		// Generates a Coefficient as they occur on a honeycomb lattice, that means, it is not inversion symmetric and it can be complex
		static Coefficient HoneyComb(const std::string& _name, const Momentum& _momentum, bool daggered, const std::optional<std::function<void(Coefficient&)>>& _custom_symmetry = std::nullopt);

		static Coefficient parse_string(const std::string& expression, bool _Q_changes_sign = false, bool _inversion_symmetry = true);
		// Calls parse_string and sets the symmetries to conform with a standard interaction
		static Coefficient parse_interaction_string(const std::string& expression);

		inline bool uses_index(const Index index) const noexcept {
			for (const auto& idx : indizes) {
				if (idx == index) return true;
			}
			return false;
		}
		inline bool depends_on_momentum() const noexcept {
			if (this->momenta.empty()) return false;
			return std::any_of(this->momenta.begin(), this->momenta.end(), [](const Momentum& momentum) {
				return !momentum.momentum_list.empty();
				});
		};
		inline bool depends_on(const MomentumSymbol::name_type momentum) const noexcept {
			if (this->momenta.empty()) return false;
			return std::any_of(this->momenta.begin(), this->momenta.end(), [momentum](const Momentum& mom) {
				return mom.isUsed(momentum) != -1;
				});
		}
		// This function determines whether the coefficient depends on something like k-l
		// Currently, this only makes sense if the coefficient does not depend on
		inline bool depends_on_two_momenta() const noexcept {
			assert(momenta.size() == 1U);
			return this->momenta.front().momentum_list.size() == 2U;
		};

		void invert_momentum(const MomentumSymbol::name_type what);

		// Utilizes V(k, k', q) = V(k', k, -q)
		void use_symmetric_interaction_exchange();
		// Utilizes V(k, k', q) = V(-k, -k', -q)
		void use_symmetric_interaction_inversion();

		void remove_momentum_contribution(const MomentumSymbol::name_type value);

		// If a function is stored in custom_symmetry, it is applied to this
		void apply_custom_symmetry();
	};

	inline bool operator==(const Coefficient& lhs, const Coefficient& rhs) {
		if (lhs.name != rhs.name) return false;
		if (lhs.momenta != rhs.momenta) return false;
		if (lhs.is_daggered != rhs.is_daggered) return false;
		return (lhs.indizes == rhs.indizes);
	}
	inline bool operator!=(const Coefficient& lhs, const Coefficient& rhs) {
		return !(lhs == rhs);
	}
}