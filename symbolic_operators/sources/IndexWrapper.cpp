#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <cassert>

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
		case Index::TypeA:
			os << "A";
			break;
		case Index::TypeB:
			os << "B";
			break;
		case Index::TypeC:
			os << "C";
			break;
		case Index::UndefinedIndex:
			os << "UNDEFINED INDEX";
			break;
		case Index::NoIndex:
			break;
		default:
			os << static_cast<char>(index);
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