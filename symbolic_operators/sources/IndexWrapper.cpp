#include <mrock/symbolic_operators/IndexWrapper.hpp>

namespace mrock::symbolic_operators {
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
		case Index::BosonA:
			os << "A";
			break;
		case Index::BosonB:
			os << "B";
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