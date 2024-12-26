#pragma once

#include <mrock/Utility/VectorWrapper.hpp>
#include <mrock/Utility/RangeUtility.hpp>
#include "IndexWrapper.hpp"

namespace mrock::SymbolicOperators {
	template<class SumIndex>
	struct SymbolicSum : public mrock::Utility::VectorWrapper<SumIndex>
	{
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->_vector;
		}

		SymbolicSum() = default;
		SymbolicSum(SumIndex sum_index) : mrock::Utility::VectorWrapper<SumIndex>(1U, sum_index) {};
		SymbolicSum(const std::vector<SumIndex>& _indizes)
			: mrock::Utility::VectorWrapper<SumIndex>(_indizes) {};
		SymbolicSum(std::vector<SumIndex>&& _indizes)
			: mrock::Utility::VectorWrapper<SumIndex>(std::move(_indizes)) {};
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
			mrock::Utility::append_vector(this->momenta, other.momenta);
			mrock::Utility::append_vector(this->spins, other.spins);
			return *this;
		}
		inline SumContainer& append(const MomentumSum& other) {
			mrock::Utility::append_vector(this->momenta, other);
			return *this;
		}
		inline SumContainer& append(const IndexSum& other) {
			mrock::Utility::append_vector(this->spins, other);
			return *this;
		}
		inline void push_back(const char momentum) {
			this->momenta.push_back(momentum);
		}
		inline void push_back(const Index spin) {
			this->spins.push_back(spin);
		}

		inline bool has_momentum() const noexcept {
			return !momenta.empty();
		}
		inline bool has_spins() const noexcept {
			return !spins.empty();
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