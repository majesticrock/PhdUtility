#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/OperatorType.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/detail/string_helper.hpp>

#include <algorithm>
#include <cassert>
#include <map>
#include <stdexcept>
#include <utility>

namespace mrock::symbolic_operators {
WickOperator::WickOperator(const OperatorType& _type,
                           const bool _is_daggered,
                           const Momentum& _momentum,
                           const IndexWrapper& _indices)
    : type(_type), is_daggered(_is_daggered), momentum(_momentum), indices(_indices) {}
WickOperator::WickOperator(const OperatorType& _type,
                           const bool _is_daggered,
                           const Momentum& _momentum,
                           const Index _index)
    : type(_type), is_daggered(_is_daggered), momentum(_momentum), indices(_index) {}

WickOperator::WickOperator(const std::string& expression) {
    // Syntax    type{Momentum_expression;index1,index2,...}(^+)

    this->type = string_to_wick.at(expression.substr(0U, expression.find('{')));
    std::vector<std::string> momentum_strings = extract_elements(expression, '{', ';');
    std::vector<std::string> index_strings = extract_elements(expression, ';', '}');

    assert(momentum_strings.size() == 1U);
    this->momentum = Momentum(momentum_strings.front());

    this->indices.reserve(index_strings.size());
    for (const auto& arg : index_strings) {
        this->indices.push_back(string_to_index.at(arg));
    }

    this->is_daggered = expression.find("^+") != std::string::npos;
}

std::vector<Operator> WickOperator::to_operator_expression() const {
    std::vector<Operator> result{Operator(this->momentum, this->indices, false),
                                 Operator(this->momentum, this->indices, false)};

    switch (this->type) {
        case OperatorType::Eta:
            result[0].momentum.add_PI = true;
            [[fallthrough]];
        case OperatorType::SC:
            result[0].momentum.flip_momentum();
            result[0].indices.insert(result[0].indices.begin(), Index::SpinDown);
            result[1].indices.insert(result[1].indices.begin(), Index::SpinUp);

            result[0].is_daggered = this->is_daggered;
            result[1].is_daggered = this->is_daggered;
            if (this->is_daggered) {
                std::swap(result[0], result[1]);
            }
            break;

        case OperatorType::CDW:
            result[0].momentum.add_PI = true;
            [[fallthrough]];
        case OperatorType::Number:
            if (this->is_daggered) {
                std::swap(result[0], result[1]);
            }
            result[0].is_daggered = true;
            result[1].is_daggered = false;
            break;

        default:
            throw std::runtime_error("Operator type not handled!");
    }

    return result;
}

bool operator==(const WickOperator& lhs, const WickOperator& rhs) {
    if (lhs.type != rhs.type)
        return false;
    if (lhs.is_daggered != rhs.is_daggered)
        return false;
    if (lhs.momentum != rhs.momentum)
        return false;
    return (lhs.indices == rhs.indices);
}

bool operator!=(const WickOperator& lhs, const WickOperator& rhs) {
    return !(lhs == rhs);
}

bool operator>(const WickOperator& lhs, const WickOperator& rhs) {
    return !(lhs <= rhs);
}

bool operator<(const WickOperator& lhs, const WickOperator& rhs) {
    if (lhs.type < rhs.type)
        return true;
    if (lhs.type > rhs.type)
        return false;

    if (lhs.indices.empty() || rhs.indices.empty())
        return false;
    if (lhs.indices[0] < rhs.indices[0])
        return true;
    if (lhs.indices[0] > rhs.indices[0])
        return false;

    return lhs.momentum < rhs.momentum;
}

bool operator>=(const WickOperator& lhs, const WickOperator& rhs) {
    return (lhs > rhs || lhs == rhs);
}

bool operator<=(const WickOperator& lhs, const WickOperator& rhs) {
    return (lhs < rhs || lhs == rhs);
}

std::ostream& operator<<(std::ostream& os, const WickOperator& op) {
    os << "\\langle " << op.type << "_{ " << op.momentum << ", ";
    for (const auto& index : op.indices) {
        os << index << " ";
    }
    os << "}";
    if (op.is_daggered) {
        os << "^\\dagger";
    }
    os << " \\rangle";
    return os;
}
std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops) {
    for (const auto& op : ops) {
        os << op << " ";
    }
    return os;
}
}  // namespace mrock::symbolic_operators