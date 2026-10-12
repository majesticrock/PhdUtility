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

Coefficient WickOperator::to_coefficient(const std::string& base) const {
    const auto type_name = wick_to_string.at(this->type);
    Coefficient result(base + "_{(" + type_name + ")}", this->momentum, this->indices, false, false,
                      this->is_daggered);
    result.custom_symmetry = std::nullopt;
    result.is_real = false;
    result.is_symmetrized_interaction = false;
    return result;
}

WickOperator WickOperator::from_coefficient(const Coefficient& coeff, const std::string& base) {
    if (coeff.momenta.size() != 1U) {
        throw std::runtime_error("WickOperator conversion requires exactly one momentum.");
    }

    const auto match_suffix = [&](const std::string& suffix) {
        if (coeff.name.size() < suffix.size() || coeff.name.compare(coeff.name.size() - suffix.size(), suffix.size(), suffix) != 0) {
            return false;
        }
        return coeff.name.substr(0U, coeff.name.size() - suffix.size()) == base;
    };

    for (const auto& [operator_type, type_string] : wick_to_string) {
        const std::string with_underscore = "_{(" + type_string + ")}";
        if (match_suffix(with_underscore)) {
            return WickOperator(operator_type, coeff.is_daggered,
                                coeff.momenta.empty() ? Momentum{} : coeff.momenta.front(), coeff.indices);
        }

        const std::string without_underscore = "{(" + type_string + ")}";
        if (match_suffix(without_underscore)) {
            return WickOperator(operator_type, coeff.is_daggered,
                                coeff.momenta.empty() ? Momentum{} : coeff.momenta.front(), coeff.indices);
        }
    }

    return WickOperator(OperatorType::Undefined, coeff.is_daggered,
                        coeff.momenta.empty() ? Momentum{} : coeff.momenta.front(), coeff.indices);
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