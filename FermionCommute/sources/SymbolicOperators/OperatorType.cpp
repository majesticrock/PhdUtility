#include "OperatorType.hpp"

namespace SymbolicOperators {
    std::ostream& operator<<(std::ostream& os, const OperatorType op)
	{
		switch (op) {
		case SC_Type:
			os << "f";
			break;
		case Eta_Type:
			os << "\\eta";
			break;
		case CDW_Type:
			os << "g";
			break;
		case Number_Type:
			os << "n";
			break;
		default:
			os << "ERROR_OPERATOR";
		}
		return os;
	}
}