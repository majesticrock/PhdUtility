#include <SymbolicOperators/Operator.hpp>

namespace SymbolicOperators {
	std::ostream& operator<<(std::ostream& os, const Operator& op)
	{
		os << "c_{ " << op.momentum << ", " << op.indizes << "}";
		if (op.is_daggered) {
			os << "^\\dagger ";
		}
		return os;
	}

	std::ostream& operator<<(std::ostream& os, const std::vector<Operator>& ops)
	{
		for (const auto& op : ops) {
			os << op;
		}
		return os;
	}

	Operator::Operator(const Momentum& _momentum, const IndexWrapper _indizes, bool _is_daggered)
		: momentum(_momentum), indizes(_indizes), is_daggered(_is_daggered) {}

	Operator::Operator(const momentum_pairs& _momentum, const IndexWrapper _indizes, bool _is_daggered)
		: momentum(_momentum), indizes(_indizes), is_daggered(_is_daggered) {}

	Operator::Operator(char _momentum, bool add_Q, const IndexWrapper _indizes, bool _is_daggered)
		: momentum(_momentum, add_Q), indizes(_indizes), is_daggered(_is_daggered) {}

	Operator::Operator(char _momentum, int sign, bool add_Q, const IndexWrapper _indizes, bool _is_daggered)
		: momentum(_momentum, sign, add_Q), indizes(_indizes), is_daggered(_is_daggered) {}
}