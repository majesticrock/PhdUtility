#pragma once
#include "Momentum.hpp"
#include <Utility/VectorWrapper.hpp>
#include <algorithm>

namespace SymbolicOperators {
	struct MomentumList : public Utility::VectorWrapper<Momentum>
	{
	private:
		using _parent = Utility::VectorWrapper<Momentum>;

	public:
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->_vector;
		};

		MomentumList() : _parent() {};
		explicit MomentumList(const Momentum& momentum)
			: _parent{ momentum } {};
		MomentumList(const Momentum& first, const Momentum& second)
			: _parent{ first, second } {};
		MomentumList(std::initializer_list<char> const& init)
			: _parent{ std::vector<Momentum>(init.size()) }
		{
			for (size_t i = 0U; i < init.size(); i++)
			{
				this->_vector[i] = Momentum(init.begin()[i]);
			}
		};

		inline MomentumList& operator*=(const int rhs) {
			for (auto& mom : _vector) {
				mom *= rhs;
			}
			return *this;
		};
		inline void multiplyMomentum(int factor) {
			(*this) *= factor;
		};
		inline void flipMomentum() {
			(*this) *= -1;
		};

		void replaceOccurances(const char replaceWhat, const Momentum& replaceWith);
		void remove_zeros();
		void flip_single(char momentum);
		inline void sort() {
			std::sort(this->begin(), this->end());
		};
	};

	std::ostream& operator<<(std::ostream& os, const MomentumList& momenta);
}