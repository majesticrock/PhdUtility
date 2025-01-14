#pragma once
#include "Momentum.hpp"
#include <mrock/utility/VectorWrapper.hpp>
#include <algorithm>

namespace mrock::symbolic_operators {
	struct MomentumList : public mrock::utility::VectorWrapper<Momentum>
	{
	private:
		using _parent = mrock::utility::VectorWrapper<Momentum>;

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
		MomentumList(std::initializer_list<Momentum> const& init)
			: _parent(init) {};
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
		inline void multiply_by(int factor) {
			(*this) *= factor;
		};
		inline void flip_momentum() {
			(*this) *= -1;
		};

		void replace_occurances(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith);
		void remove_zeros();
		void flip_single(const MomentumSymbol::name_type momentum);
		inline void sort() {
			std::sort(this->begin(), this->end());
		};
	};

	std::ostream& operator<<(std::ostream& os, const MomentumList& momenta);
}