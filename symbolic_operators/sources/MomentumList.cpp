#include <mrock/symbolic_operators/MomentumList.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/MomentumSymbol.hpp>
#include <cstddef>
#include <ostream>

namespace mrock::symbolic_operators {
	void MomentumList::replace_occurances(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith) {
		for (auto& mom : momenta) {
			mom.replace_occurances(replaceWhat, replaceWith);
		}
	}

	void MomentumList::remove_zeros() {
		for (auto& mom : momenta) {
			mom.remove_zeros();
		}
	}

	void MomentumList::flip_single(const MomentumSymbol::name_type momentum) {
		for (auto& mom : momenta) {
			mom.flip_single(momentum);
		}
	}

	std::ostream& operator<<(std::ostream& os, const MomentumList& momenta) {
		if (momenta.empty()) return os;
		os << "( ";
		for (std::size_t i = 0U; i < momenta.size(); ++i)
		{
			os << momenta[i];
			if (i + 1U < momenta.size()) {
				os << ", ";
			}
		}
		os << " )";
		return os;
	}

	MomentumList::MomentumList() {}

	MomentumList::MomentumList(const Momentum& momentum)
		: momenta{ momentum } {}

	MomentumList::MomentumList(const Momentum& first, const Momentum& second)
		: momenta{ first, second } {}

	MomentumList::MomentumList(std::initializer_list<Momentum> init)
		: momenta(init) {}

	MomentumList::MomentumList(std::initializer_list<char> init)
		: momenta{ std::vector<Momentum>(init.size()) }
	{
		for (std::size_t i = 0U; i < init.size(); i++)
		{
			this->momenta[i] = Momentum(init.begin()[i]);
		}
	}
}