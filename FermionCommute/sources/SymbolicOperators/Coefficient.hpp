#pragma once
#include "Operator.hpp"
#include "IndexWrapper.hpp"
#include "MomentumList.hpp"

namespace SymbolicOperators {
	struct Coefficient {
		std::string name;
		MomentumList momenta;
		// Contains all indizes, standard: first index = spin, all others arbitrary, e.g. orbitals, bands etc
		IndexWrapper indizes;
		// if Coeff(k) = Coeff(-k)
		bool translationalInvariance = true;
		// if Coeff(k+Q) = -Coeff(k)
		bool Q_changes_sign{};
		bool isDaggered{};

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& name;
			ar& momenta;
			ar& indizes;
			ar& isDaggered;
			ar& translationalInvariance;
			ar& Q_changes_sign;
		}

		Coefficient();
		explicit Coefficient(std::string _name);
		Coefficient(std::string _name, const Momentum& _momentum, const IndexWrapper& _indizes, bool _Q_changes_sign = false, bool _isDaggered = false);
		Coefficient(std::string _name, const Momentum& _momentum, bool _Q_changes_sign = false, bool _isDaggered = false);
		Coefficient(std::string _name, const MomentumList& _momenta, const IndexWrapper& _indizes = IndexWrapper{}, bool _Q_changes_sign = false, bool _isDaggered = false);

		static Coefficient parse_string(const std::string& expression);

		inline bool usesIndex(const Index index) const noexcept {
			for (const auto& idx : indizes) {
				if (idx == index) return true;
			}
			return false;
		}
		inline bool dependsOnMomentum() const noexcept {
			if (this->momenta.empty()) return false;
			return std::any_of(this->momenta.begin(), this->momenta.end(), [](const Momentum& momentum) {
				return !momentum.momentum_list.empty();
				});
		};
		inline bool dependsOn(char momentum) const noexcept {
			if (this->momenta.empty()) return false;
			return std::any_of(this->momenta.begin(), this->momenta.end(), [momentum](const Momentum& mom) {
				return mom.isUsed(momentum) != -1;
				});
		}
		// This function determines whether the coefficient depends on something like k-l
		// Currently, this only makes sense if the coefficient does not depend on
		inline bool dependsOnTwoMomenta() const noexcept {
			assert(momenta.size() == 1U);
			return this->momenta.front().momentum_list.size() == 2U;
		};
	};

	inline bool operator==(const Coefficient& lhs, const Coefficient& rhs) {
		if (lhs.name != rhs.name) return false;
		if (lhs.momenta != rhs.momenta) return false;
		if (lhs.isDaggered != rhs.isDaggered) return false;
		return (lhs.indizes == rhs.indizes);
	}
	inline bool operator!=(const Coefficient& lhs, const Coefficient& rhs) {
		return !(lhs == rhs);
	}
}