#include <SymbolicOperators/IndexWrapper.hpp>

namespace SymbolicOperators {
	std::ostream& operator<<(std::ostream& os, const Index index)
	{
		switch (index)
		{
		case Index::SpinUp:
			os << "\\uparrow";
			break;
		case Index::SpinDown:
			os << "\\downarrow";
			break;
		case Index::Sigma:
			os << "\\sigma";
			break;
		case Index::SigmaPrime:
			os << "\\sigma'";
			break;
		case Index::GeneralSpin_S:
			os << "S";
			break;
		case Index::GeneralSpin_SPrime:
			os << "S'";
			break;
		case Index::NoIndex:
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