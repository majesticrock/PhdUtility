#pragma once
#include <utility>
#include <iostream>

namespace SymbolicOperators {
	template<typename T>
	struct KroneckerDelta {
		T first{};
		T second{};

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& first;
			ar& second;
		}

		constexpr bool isOne() const {
			return first == second;
		};
	};

	template<typename T>
	constexpr auto make_delta(const T& first, const T& second) {
		return KroneckerDelta<T>{first, second};
	}

	template<typename T>
	constexpr auto make_delta(std::decay_t<T>&& first, std::decay_t<T>&& second) {
		return KroneckerDelta<std::decay_t<T>>{std::move(first), std::move(second)};
	}

	template<typename T>
	bool operator==(const KroneckerDelta<T>& lhs, const KroneckerDelta<T>& rhs) {
		if (lhs.first == rhs.first && lhs.second == rhs.second) return true;
		if (lhs.first == rhs.second && lhs.second == rhs.first) return true;
		return false;
	};
	template<typename T>
	bool operator!=(const KroneckerDelta<T>& lhs, const KroneckerDelta<T>& rhs) {
		return !(lhs == rhs);
	};

	template<typename T>
	inline std::ostream& operator<<(std::ostream& os, const KroneckerDelta<T>& delta) {
		os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		return os;
	};
}