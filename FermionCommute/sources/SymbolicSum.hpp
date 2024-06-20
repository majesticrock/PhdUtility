#pragma once

#include <Utility/VectorWrapper.hpp>
#include <Utility/RangeUtility.hpp>
#include "IndexWrapper.hpp"

namespace SymbolicOperators {
	template<class SumIndex>
	struct SymbolicSum : public Utility::VectorWrapper<SumIndex>
	{
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->_vector;
		}

		SymbolicSum() = default;
		SymbolicSum(SumIndex sum_index) : Utility::VectorWrapper<SumIndex>(1U, sum_index) {};
		SymbolicSum(const std::vector<SumIndex>& _indizes)
			: Utility::VectorWrapper<SumIndex>(_indizes) {};
		SymbolicSum(std::vector<SumIndex>&& _indizes)
			: Utility::VectorWrapper<SumIndex>(std::move(_indizes)) {};
	};

	template<class SumIndex>
	std::ostream& operator<<(std::ostream& os, SymbolicSum<SumIndex> const& sum) {
		if (sum.empty()) return os;
		os << "\\sum_{ ";
		for (const auto& index : sum) {
			os << index << " ";
		}
		os << "} ";
		return os;
	}

	typedef SymbolicSum<Index> IndexSum;
	typedef SymbolicSum<char> MomentumSum;

	struct SumContainer {
		MomentumSum momenta;
		IndexSum spins;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->momenta;
			ar& this->spins;
		};

		inline SumContainer& append(const SumContainer& other) {
			Utility::append_vector(this->momenta, other.momenta);
			Utility::append_vector(this->spins, other.spins);
			return *this;
		}
		inline SumContainer& append(const MomentumSum& other) {
			Utility::append_vector(this->momenta, other);
			return *this;
		}
		inline SumContainer& append(const IndexSum& other) {
			Utility::append_vector(this->spins, other);
			return *this;
		}
		inline void push_back(const char momentum) {
			this->momenta.push_back(momentum);
		}
		inline void push_back(const Index spin) {
			this->spins.push_back(spin);
		}
	};

	inline bool operator==(const SumContainer& lhs, const SumContainer& rhs) {
		return (lhs.momenta == rhs.momenta && lhs.spins == rhs.spins);
	}
	inline bool operator!=(const SumContainer& lhs, const SumContainer& rhs) {
		return !(lhs == rhs);
	}

	inline std::ostream& operator<<(std::ostream& os, const SumContainer& sums) {
		os << sums.momenta << sums.spins;
		return os;
	}
}