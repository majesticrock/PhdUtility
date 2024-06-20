#pragma once
#include "KroneckerDelta.hpp"
#include "Momentum.hpp"
#include "IndexWrapper.hpp"
#include <Utility/defines_arithmetic_operators.hpp>

namespace SymbolicOperators {
	template <class T>
	void remove_delta_squared(std::vector<KroneckerDelta<T>>& deltas) {
		using predicate_type = std::conditional_t<Utility::is_linearly_combinable_v<T>(), KroneckerDelta<T>, const KroneckerDelta<T>&>;
		auto new_end = std::remove_if(deltas.begin(), deltas.end(), [](predicate_type delta) {
			if constexpr (Utility::is_linearly_combinable_v<T>()) {
				delta.first -= delta.second;
				delta.second = T{};
			}
			return delta.first == delta.second;
			});
		deltas.erase(new_end, deltas.end());
	}

	template <class T>
	void remove_delta_is_one(std::vector<KroneckerDelta<T>>& deltas) {
		auto new_end = std::remove_if(deltas.begin(), deltas.end(), [](const KroneckerDelta<T>& delta) {
			return delta.isOne();
			});
		deltas.erase(new_end, deltas.end());
	}

	inline bool is_always_zero(const std::vector<KroneckerDelta<Index>>& deltas) {
		return std::any_of(deltas.begin(), deltas.end(), [](const KroneckerDelta<Index>& delta) {
			return (delta.first != delta.second && (!is_mutable(delta.first) && !is_mutable(delta.second)));
			});
	}
	inline bool is_always_zero(const std::vector<KroneckerDelta<Momentum>>& deltas) {
		return std::any_of(deltas.begin(), deltas.end(), [](const KroneckerDelta<Momentum>& delta) {
			return delta.first.differsOnlyInQ(delta.second);
			});
	}
}