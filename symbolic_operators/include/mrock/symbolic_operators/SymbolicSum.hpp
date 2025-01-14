#pragma once

#include <concepts>
#include <mrock/utility/VectorWrapper.hpp>
#include <mrock/utility/RangeUtility.hpp>
#include "IndexWrapper.hpp"
#include "MomentumSymbol.hpp"

namespace mrock::symbolic_operators {
	template<class SumIndex>
	struct SymbolicSum : public mrock::utility::VectorWrapper<SumIndex>
	{
	public:
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->_vector;
		}

		SymbolicSum() = default;
		SymbolicSum(SumIndex sum_index) : mrock::utility::VectorWrapper<SumIndex>(1U, sum_index) {}
		SymbolicSum(const std::vector<SumIndex>& _indizes)
			: mrock::utility::VectorWrapper<SumIndex>(_indizes) {}
		SymbolicSum(std::vector<SumIndex>&& _indizes)
			: mrock::utility::VectorWrapper<SumIndex>(std::move(_indizes)) {}
		
		template<std::convertible_to<SumIndex> input_type> requires (!std::same_as<SumIndex, input_type>)
		SymbolicSum(std::initializer_list<input_type> init) { 
			this->_vector.reserve(init.size());
			for (const auto& x : init) {
				this->_vector.emplace_back(x);
			}
		}
		SymbolicSum(std::initializer_list<SumIndex> init) : mrock::utility::VectorWrapper<SumIndex>(std::move(init)) {}

		inline bool is_summed_over(SumIndex what) const {
			for (const auto& _s : this->_vector) {
				if (_s == what) return true;
			}
			return false;
		}
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
	typedef SymbolicSum<MomentumSymbol::name_type> MomentumSum;

	struct SumContainer {
		MomentumSum momenta;
		IndexSum spins;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->momenta;
			ar& this->spins;
		};

		inline SumContainer& append(const SumContainer& other) {
			mrock::utility::append_vector(this->momenta, other.momenta);
			mrock::utility::append_vector(this->spins, other.spins);
			return *this;
		}
		inline SumContainer& append(const MomentumSum& other) {
			mrock::utility::append_vector(this->momenta, other);
			return *this;
		}
		inline SumContainer& append(const IndexSum& other) {
			mrock::utility::append_vector(this->spins, other);
			return *this;
		}
		inline void push_back(const MomentumSymbol::name_type momentum) {
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