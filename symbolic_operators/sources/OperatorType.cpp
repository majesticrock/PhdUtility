#include <mrock/symbolic_operators/OperatorType.hpp>

namespace mrock::symbolic_operators {
std::ostream& operator<<(std::ostream& os, const OperatorType op) {
    switch (op) {
        case OperatorType::SC:
            os << "\\hat{f}";
            break;
        case OperatorType::Eta:
            os << "\\hat{\\eta}";
            break;
        case OperatorType::CDW:
            os << "\\hat{g}";
            break;
        case OperatorType::Number:
            os << "\\hat{n}";
            break;
        default:
            os << "ERROR_OPERATOR";
    }
    return os;
}
}  // namespace mrock::symbolic_operators