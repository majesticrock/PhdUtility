#include <mrock/symbolic_operators/OperatorType.hpp>

namespace mrock::symbolic_operators {
    std::ostream& operator<<(std::ostream& os, const OperatorType op)
	{
		switch (op) {
		case SC_Type:
			os << "\\hat{f}";
			break;
		case Eta_Type:
			os << "\\hat{\\eta}";
			break;
		case CDW_Type:
			os << "\\hat{g}";
			break;
		case Number_Type:
			os << "\\hat{n}";
			break;
		default:
			os << "ERROR_OPERATOR";
		}
		return os;
	}
}