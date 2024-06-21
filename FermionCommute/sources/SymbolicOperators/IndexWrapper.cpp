#include "IndexWrapper.hpp"

namespace SymbolicOperators {
	std::ostream& operator<<(std::ostream& os, const Index index)
	{
		switch (index)
		{
		case SpinUp:
			os << "\\uparrow";
			break;
		case SpinDown:
			os << "\\downarrow";
			break;
		case Sigma:
			os << "\\sigma";
			break;
		case SigmaPrime:
			os << "\\sigma'";
			break;
		default:
			os << "ERROR_INDEX";
			break;
		}
		return os;
	}

	std::ostream& operator<<(std::ostream& os, const IndexWrapper& indizes) {
		for (const auto& idx : indizes) {
			os << idx << " ";
		}
		return os;
	}
}