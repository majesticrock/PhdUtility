#include <mrock/symbolic_operators/IndexWrapper.hpp>

namespace mrock::symbolic_operators {
void IndexWrapper::replace_index(Index target, Index replace_with) {
    for (auto& index : indizes) {
        if (index == target) {
            index = replace_with;
        }
    }
}

std::ostream& operator<<(std::ostream& os, const Index index) {
    switch (index) {
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

        case Index::UndefinedIndex:
            os << "UNDEFINED INDEX";
            break;
        case Index::NoIndex:
            break;
        default:
            index_base c = static_cast<index_base>(index);
            if (is_mutable(index)) {
                os << "S";
                c -= static_cast<index_base>(Index::GeneralSpin_S);
                while (c-- > 0) {
                    os << "'";
                }
            } else {
                os << static_cast<index_base>(index);
            }
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
}  // namespace mrock::symbolic_operators