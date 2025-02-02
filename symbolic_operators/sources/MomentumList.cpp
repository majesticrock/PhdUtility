#include <mrock/symbolic_operators/MomentumList.hpp>

namespace mrock::symbolic_operators {
	void MomentumList::replace_occurances(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith) {
		for (auto& mom : _vector) {
			mom.replace_occurances(replaceWhat, replaceWith);
		}
	}

	void MomentumList::remove_zeros() {
		for (auto& mom : _vector) {
			mom.remove_zeros();
		}
	}

	void MomentumList::flip_single(const MomentumSymbol::name_type momentum) {
		for (auto& mom : _vector) {
			mom.flip_single(momentum);
		}
	}

	std::ostream& operator<<(std::ostream& os, const MomentumList& momenta) {
		if (momenta.empty()) return os;
		os << "( ";
		for (size_t i = 0U; i < momenta.size(); ++i)
		{
			os << momenta[i];
			if (i + 1U < momenta.size()) {
				os << ", ";
			}
		}
		os << " )";
		return os;
	}

	MomentumList::MomentumList() : _parent() {}

	MomentumList::MomentumList(const Momentum& momentum)
		: _parent{ momentum } {}

	MomentumList::MomentumList(const Momentum& first, const Momentum& second)
		: _parent{ first, second } {}

	MomentumList::MomentumList(std::initializer_list<Momentum> init)
		: _parent(init) {}

	MomentumList::MomentumList(std::initializer_list<char> init)
		: _parent{ std::vector<Momentum>(init.size()) }
	{
		for (size_t i = 0U; i < init.size(); i++)
		{
			this->_vector[i] = Momentum(init.begin()[i]);
		}
	}
}