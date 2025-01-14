#pragma once
#include <iostream>
#include "Momentum.hpp"
#include "IndexWrapper.hpp"
#include "OperatorType.hpp"

namespace mrock::symbolic_operators {
	struct WickOperator {
		OperatorType type;
		bool is_daggered{};
		Momentum momentum;
		IndexWrapper indizes;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& type;
			ar& is_daggered;
			ar& momentum;
			ar& indizes;
		}

		WickOperator(const OperatorType& _type, const bool _is_daggered, const Momentum& _momentum, const IndexWrapper& _indizes = IndexWrapper());
		WickOperator(const OperatorType& _type, const bool _is_daggered, const Momentum& _momentum, const Index _index);
		WickOperator();
		WickOperator(const std::string& expression);

		inline bool uses_index(const Index index) const noexcept {
			for (const auto& idx : this->indizes) {
				if (idx == index) return true;
			}
			return false;
		};
		inline bool depends_on(const MomentumSymbol::name_type momentum) const noexcept {
			return this->momentum.isUsed(momentum) != -1;
		}
		inline void remove_momentum_contribution(const MomentumSymbol::name_type value) {
			momentum.remove_contribution(value);
		}
	};

	std::ostream& operator<<(std::ostream& os, const WickOperator& op);
	std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops);
}