#include <mrock/symbolic_operators/MomentumList.hpp>

namespace mrock::symbolic_operators {
	void MomentumList::replace_occurances(const char replaceWhat, const Momentum& replaceWith) {
		for (auto& mom : _vector) {
			mom.replace_occurances(replaceWhat, replaceWith);
		}
	}

	void MomentumList::remove_zeros() {
		for (auto& mom : _vector) {
			mom.remove_zeros();
		}
	}

	void MomentumList::flip_single(char momentum) {
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
}