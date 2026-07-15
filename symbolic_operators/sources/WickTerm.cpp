#include <mrock/symbolic_operators/AbstractTerm.hpp>
#include <mrock/symbolic_operators/Coefficient.hpp>
#include <mrock/symbolic_operators/Fractional.hpp>
#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <mrock/symbolic_operators/KroneckerDelta.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/MomentumList.hpp>
#include <mrock/symbolic_operators/MomentumSymbol.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/OperatorType.hpp>
#include <mrock/symbolic_operators/SumContainer.hpp>
#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/symbolic_operators/detail/container_helper.hpp>
#include <mrock/symbolic_operators/detail/string_helper.hpp>

#include <cassert>
#include <cctype>
#include <compare>
#include <cstddef>
#include <map>
#include <ostream>
#include <utility>

#define LEFT temporary_operators[i]
#define RIGHT temporary_operators[i + 1]
#define L_SPIN temporary_operators[i].first_index()
#define R_SPIN temporary_operators[i + 1].first_index()

namespace mrock::symbolic_operators {
// Constructors
WickTerm::WickTerm(const Term* base)
    : AbstractTerm<WickOperator>(base->multiplicity,
                                 base->coefficients,
                                 base->sums,
                                 base->delta_momenta,
                                 base->delta_indizes,
                                 std::vector<WickOperator>()),
      temporary_operators() {}

WickTerm::WickTerm(const Term& base)
    : AbstractTerm<WickOperator>(base.multiplicity,
                                 base.coefficients,
                                 base.sums,
                                 base.delta_momenta,
                                 base.delta_indizes,
                                 std::vector<WickOperator>()),
      temporary_operators() {}

WickTerm::WickTerm(const WickTerm& base, const TemplateResult::SingleResult& result)
    : AbstractTerm<WickOperator>(base.multiplicity,
                                 base.coefficients,
                                 base.sums,
                                 base.delta_momenta,
                                 base.delta_indizes,
                                 base.operators),
      temporary_operators() {
    this->operators.push_back(result.op);
    this->delta_indizes.insert(this->delta_indizes.end(), result.index_deltas.begin(), result.index_deltas.end());
}

WickTerm::WickTerm(const std::string& expression) : AbstractTerm<WickOperator>(1) {
    // Syntax
    // [factor] [index_sum] [momentum_sum] [coefficients...] [momentum_deltas...] [index_deltas...] [operators...]
    /* factor needs to be an integer
     *  Sums must be "sum:index{index1,index2,...}" or "sum:momentum{momentum_name1,momentum_name2,...}"
     *  coefficient must be "c:name{Momentum_expression1,...;index1,index2,...}"
     *  deltas must be "delta:momentum{Momentum_expression,Momentum_expression}" or "delta:index{Index,Index}"
     *  operators must be "o:type{Momentum_expression;index1,index2,...}(^+)"
     */

    std::size_t pos{};
    if (std::isdigit(expression.front()) || expression.front() == '-' || expression.front() == '+') {
        pos = expression.find(' ');
        if (pos != std::string::npos) {
            this->multiplicity = std::stoi(expression.substr(0U, pos));
        }
    }
    std::size_t new_pos{};
    ++pos;
    while (new_pos != std::string::npos) {
        new_pos = expression.find(' ', pos);
        string_parser(expression.substr(pos, new_pos - pos));
        pos = new_pos + 1;
    }
}

// Member functions
void WickTerm::string_parser(std::string&& expression) {
    std::size_t forward_pos{expression.find(' ')};
    if (forward_pos == std::string::npos)
        forward_pos = expression.size();

    const std::string sub = expression.substr(0U, forward_pos);
    const std::size_t sub_delimiter = sub.find(':');

    if (sub_delimiter == std::string::npos)
        throw std::invalid_argument("Did not find ':' in " + expression);

    if (sub.substr(0U, sub_delimiter) == "sum") {
        const std::string type =
            expression.substr(sub_delimiter + 1, expression.find('{', sub_delimiter) - sub_delimiter - 1);
        const std::vector<std::string> argument_list = extract_elements(expression);

        if (type == "index") {
            this->sums.spins.reserve(argument_list.size());
            for (const auto& arg : argument_list) {
                this->sums.spins.push_back(string_to_index.at(arg));
            }
        } else if (type == "momentum") {
            this->sums.momenta.reserve(argument_list.size());
            for (const auto& arg : argument_list) {
                assert(arg.size() == 1U);
                this->sums.momenta.push_back(arg.front());
            }
        } else {
            throw std::invalid_argument("Sum type not recognized " + type + " in expression " + expression);
        }
    } else if (sub.substr(0U, sub_delimiter) == "delta") {
        const std::string type =
            expression.substr(sub_delimiter + 1, expression.find('{', sub_delimiter) - sub_delimiter - 1);
        const std::vector<std::string> argument_list = extract_elements(expression);
        assert(argument_list.size() == 2U);

        if (type == "index") {
            this->delta_indizes.push_back(
                make_delta(string_to_index.at(argument_list[0]), string_to_index.at(argument_list[1])));
        } else if (type == "momentum") {
            this->delta_momenta.push_back(make_delta(Momentum(argument_list[0]), Momentum(argument_list[1])));
        } else {
            throw std::invalid_argument("Delta type not recognized " + type + " in expression " + expression);
        }
    } else if (sub.substr(0U, sub_delimiter) == "c") {
        this->coefficients.push_back(Coefficient::parse_string(sub.substr(sub_delimiter + 1)));
    } else if (sub.substr(0U, sub_delimiter) == "o") {
        this->operators.push_back(WickOperator(sub.substr(sub_delimiter + 1)));
    } else {
        throw std::invalid_argument("Did not parse expression <" + expression + "> at <" + sub + "> with delimiter " +
                                    std::to_string(sub_delimiter));
    }
}

bool WickTerm::resolve_deltas() {
    if (!resolve_momentum_deltas())
        return false;
    if (!resolve_index_deltas())
        return false;

    return true;
}

void WickTerm::rename_sums() {
    AbstractTerm<WickOperator>::rename_sums();

    for (const auto& sum : sums.momenta) {
        for (auto& op : operators) {
            int index = op.momentum.is_used_at(sum);
            if (index < 0)
                continue;
            if (op.momentum.momentum_list.size() == 1)
                break;

            Momentum buffer = op.momentum;
            if (buffer.momentum_list[index].factor > 0)
                buffer.flip_momentum();
            buffer.momentum_list[index].factor *= -1;
            buffer.momentum_list[index].name = buffer_list[0];

            for (auto& op2 : operators) {
                op2.momentum.replace_occurances(sum, buffer);
                op2.momentum.replace_occurances(buffer_list[0], Momentum(sum));
            }
            for (auto& coeff : coefficients) {
                coeff.momenta.replace_occurances(sum, buffer);
                coeff.momenta.replace_occurances(buffer_list[0], Momentum(sum));
            }
        }
    }
    discard_zero_momenta();
}

void WickTerm::discard_zero_momenta() {
    for (auto& op : operators) {
        op.momentum.remove_zeros();
    }
    for (auto& coeff : coefficients) {
        coeff.momenta.remove_zeros();
    }
}

void WickTerm::sort() {
    for (auto& delta : delta_momenta) {
        if (delta.first.momentum_list.size() == 1 && delta.second.momentum_list.size() == 1) {
            // This comparison is well defined because we save the momentum as char i.e. byte
            // which is easily comparable
            if (delta.first.momentum_list[0].name < delta.second.momentum_list[0].name) {
                std::swap(delta.first, delta.second);
                if (delta.first.momentum_list[0].factor < 0) {
                    delta.first.flip_momentum();
                    delta.second.flip_momentum();
                }
                if (delta.first.add_Q) {
                    delta.first.add_Q = false;
                    delta.second.add_Q = !(delta.second.add_Q);
                }
            }
            for (auto& op : operators) {
                op.momentum.replace_occurances(delta.first.momentum_list[0].name, delta.second);
            }
            for (auto& coeff : coefficients) {
                coeff.momenta.replace_occurances(delta.first.momentum_list[0].name, delta.second);
            }
        }
    }

    for (auto& op : operators) {
        if (op.type == OperatorType::CDW && op.momentum.add_Q) {
            op.momentum.add_Q = false;
            op.is_daggered = !(op.is_daggered);
        }
    }

    for (std::size_t i = 0U; i < operators.size(); ++i) {
        for (std::size_t j = i + 1U; j < operators.size(); ++j) {
            if (operators[i].type > operators[j].type) {
                std::swap(operators[i], operators[j]);
            } else if (operators[i].type == operators[j].type) {
                if (momentum_order(operators[i].momentum, operators[j].momentum)) {
                    std::swap(operators[i], operators[j]);
                }
            }
        }
    }

    for (auto& coeff : coefficients) {
        for (auto& momentum : coeff.momenta) {
            momentum.sort();

            if (coeff.inversion_symmetry && !momentum.momentum_list.empty()) {
                if (momentum.momentum_list[0].factor < 0) {
                    momentum.flip_momentum();
                }
            }
            if (coeff.Q_changes_sign && momentum.add_Q) {
                momentum.add_Q = false;
                this->multiplicity *= -1;
            }
        }
    }

    for (auto& coeff : coefficients) {
        if (coeff.momenta.size() == 3U) {
            Momentum* first_momentum = &coeff.momenta.front();
            if (first_momentum != nullptr && first_momentum->empty()) {
                if (coeff.momenta.size() > 1U)
                    first_momentum = &coeff.momenta[1];
            }
            if ((first_momentum != nullptr) && (!first_momentum->empty())) {
                if (!first_momentum->first_momentum_is('k')) {
                    coeff.use_symmetric_interaction_exchange();
                }
                if (sums.momenta.empty()) {
                    if (coeff.momenta.back().first_momentum_is_negative()) {
                        coeff.use_symmetric_interaction_inversion();
                    }
                } else if ((!sums.momenta.is_summed_over(first_momentum->front().name)) &&
                           first_momentum->first_momentum_is_negative()) {
                    coeff.use_symmetric_interaction_inversion();
                }
            }
        }

        for (auto& momentum : coeff.momenta) {
            if (momentum.empty())
                continue;
            if ((!sums.momenta.is_summed_over(momentum.front().name)) && momentum.front().factor < 0) {
                if (coeff.inversion_symmetry)
                    momentum.flip_momentum();
            }

            for (const auto& sum : sums.momenta) {
                int idx = momentum.is_used_at(sum);
                if (idx < 0)
                    continue;

                if (momentum.momentum_list[idx].factor < 0) {
                    invert_momentum_sum(sum);
                }
            }
        }
    }
}

void WickTerm::include_template_result(const TemplateResult::SingleResult& result) {
    this->delta_indizes.insert(this->delta_indizes.begin(), result.index_deltas.begin(), result.index_deltas.end());
    this->operators.push_back(result.op);
    this->multiplicity *= result.factor;
}

bool WickTerm::is_pauli_forbidden() const {
    // Cannot be forbidden if there is only one WickOperator
    if (operators.size() < 2U)
        return false;
    // IMPORTANT: This function assumes that the term prior to the application of
    // Wick's theorem was normal ordered!

    std::vector<Operator> transformed_operators;
    transformed_operators.reserve(2 * operators.size());

    for (const auto& wick_op : operators) {
        append_vector(transformed_operators, wick_op.to_operator_expression());
    }

    for (std::size_t i = 0U; i < transformed_operators.size(); ++i) {
        for (std::size_t j = i + 1; j < transformed_operators.size(); ++j) {
            if (transformed_operators[i] == transformed_operators[j]) {
                return true;
            }
        }
    }

    return false;
}

int WickTerm::which_operator_depends_on(const MomentumSymbol::name_type momentum) const noexcept {
    for (std::size_t i = 0U; i < operators.size(); ++i) {
        if (operators[i].depends_on(momentum))
            return i;
    }
    return -1;
}

bool WickTerm::uses_index(const Index index) const noexcept {
    for (const auto& op : operators) {
        if (op.uses_index(index))
            return true;
    }
    for (const auto& coeff : coefficients) {
        if (coeff.uses_index(index))
            return true;
    }
    return false;
}

bool WickTerm::includes_type(const OperatorType operator_type) const {
    return std::any_of(this->operators.begin(), this->operators.end(),
                       [operator_type](const WickOperator& op) { return op.type == operator_type; });
}
bool WickTerm::has_single_coefficient() const noexcept {
    return this->coefficients.size() == 1U;
}
bool WickTerm::is_bilinear() const noexcept {
    return this->operators.size() == 1U;
}
bool WickTerm::is_quartic() const noexcept {
    return this->operators.size() == 2U;
}
double WickTerm::get_factor() const noexcept {
    return static_cast<double>(this->multiplicity);
}
const Coefficient& WickTerm::get_first_coefficient() const {
    assert(!(this->coefficients.empty()));
    return this->coefficients.front();
}
bool WickTerm::handled() const noexcept {
    if (this->temporary_operators.empty())
        return true;
    return !(this->operators.empty());
}

bool operator==(const WickOperator& lhs, const WickOperator& rhs) {
    if (lhs.type != rhs.type)
        return false;
    if (lhs.is_daggered != rhs.is_daggered)
        return false;
    if (lhs.momentum != rhs.momentum)
        return false;
    return (lhs.indizes == rhs.indizes);
}
bool operator!=(const WickOperator& lhs, const WickOperator& rhs) {
    return !(lhs == rhs);
}
bool operator==(const WickTerm& lhs, const WickTerm& rhs) {
    if (lhs.coefficients != rhs.coefficients)
        return false;
    if (lhs.sums != rhs.sums)
        return false;
    if (lhs.delta_indizes != rhs.delta_indizes)
        return false;
    if (lhs.delta_momenta != rhs.delta_momenta)
        return false;
    if (lhs.operators != rhs.operators)
        return false;
    return true;
}
bool operator!=(const WickTerm& lhs, const WickTerm& rhs) {
    return !(lhs == rhs);
}

// Operator overloads
std::ostream& operator<<(std::ostream& os, const WickTerm& term) {
    if (term.multiplicity > 0) {
        os << "+";
    }
    os << term.multiplicity << " ";
    os << term.sums;
    os << term.coefficients << " ";
    for (const auto& delta : term.delta_momenta) {
        os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
    }
    for (const auto& delta : term.delta_indizes) {
        os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
    }
    if (term.is_identity()) {
        os << " \\hat{1} ";
        return os;
    }
    for (const auto& op : term.operators) {
        os << op << " ";
    }
    return os;
}
std::ostream& operator<<(std::ostream& os, const WickTermCollector& terms) {
    for (WickTermCollector::const_iterator it = terms.begin(); it != terms.end(); ++it) {
        os << "\t&" << *it;
        if (it != terms.end() - 1) {
            os << " \\\\";
        }
        os << "\n";
    }
    return os;
}
WickTermCollector& operator+=(WickTermCollector& lhs, const WickTerm& rhs) {
    for (auto it = lhs.begin(); it != lhs.end(); ++it) {
        if (*it == rhs) {
            it->multiplicity += rhs.multiplicity;
            if (it->multiplicity == 0)
                lhs.erase(it);
            return lhs;
        }
    }
    lhs.push_back(rhs);
    return lhs;
}
WickTermCollector& operator-=(WickTermCollector& lhs, const WickTerm& rhs) {
    for (auto it = lhs.begin(); it != lhs.end(); ++it) {
        if (*it == rhs) {
            it->multiplicity -= rhs.multiplicity;
            if (it->multiplicity == 0)
                lhs.erase(it);
            return lhs;
        }
    }
    lhs.push_back(rhs);
    return lhs;
}
WickTermCollector& operator+=(WickTermCollector& lhs, const WickTermCollector& rhs) {
    for (const auto& term : rhs) {
        lhs += term;
    }
    return lhs;
}
WickTermCollector& operator-=(WickTermCollector& lhs, const WickTermCollector& rhs) {
    for (const auto& term : rhs) {
        lhs -= term;
    }
    return lhs;
}
}  // namespace mrock::symbolic_operators