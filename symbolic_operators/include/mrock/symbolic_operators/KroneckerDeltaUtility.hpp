#pragma once
#include "KroneckerDelta.hpp"
#include "Momentum.hpp"
#include "IndexWrapper.hpp"
#include <mrock/utility/defines_arithmetic_operators.hpp>

namespace mrock::symbolic_operators {
	template <class T>
	void remove_delta_squared(std::vector<KroneckerDelta<T>>& deltas) {
		using predicate_type = std::conditional_t<mrock::utility::is_linearly_combinable_v<T>(), KroneckerDelta<T>, const KroneckerDelta<T>&>;
		auto new_end = std::remove_if(deltas.begin(), deltas.end(), [](predicate_type delta) {
			if constexpr (mrock::utility::is_linearly_combinable_v<T>()) {
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
			return delta.first.differs_only_in_Q(delta.second);
			});
	}
	inline void remove_double_occurances(KroneckerDelta<Momentum>& delta) {
		if(delta.first.add_Q) {
			delta.first.add_Q = false;
			delta.second.add_Q = !(delta.second.add_Q);
		}
		for (auto it = delta.first.momentum_list.begin(); it != delta.first.momentum_list.end();) {
			const int idx = delta.second.is_used_at(it->name);
			if(idx < 0) {
				++it;
				continue;
			}

			delta.second.momentum_list[idx].factor -= it->factor;
			it = delta.first.momentum_list.erase(it);
			if(delta.second.momentum_list[idx].factor == 0) {
				delta.second.momentum_list.erase(delta.second.momentum_list.begin() + idx);
			}
		}
	}
} // namespace mrock::symbolic_operators