#include <mrock/symbolic_operators/ExpectationValues>

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace mrock::symbolic_operators;

namespace {
struct NumericState {
    std::map<unsigned char, double> momentum_values;
    std::map<std::pair<double, Index>, double> number_values;
    std::map<double, double> pair_values;
};

double momentum_value(const Momentum& momentum, const NumericState& state) {
    if (momentum.add_PI)
        throw std::invalid_argument("This test evaluator does not assign a value to Pi momentum.");

    double value{};
    for (const auto& symbol : momentum) {
        value += symbol.factor * state.momentum_values.at(static_cast<unsigned char>(symbol.name));
    }
    return value;
}

double number_value(const WickOperator& op, const NumericState& state) {
    return state.number_values.at({momentum_value(op.momentum, state), op.indices[0]});
}

double expectation_value(const WickOperator& op, const NumericState& state) {
    switch (op.type) {
        case OperatorType::Number:
            return number_value(op, state);
        case OperatorType::HoleNumber:
            return 1.0 - number_value(op, state);
        case OperatorType::SC:
            return state.pair_values.at(momentum_value(op.momentum, state));
        default:
            throw std::invalid_argument("No numerical value assigned to this Wick operator type.");
    }
}

double evaluate(const WickTermCollector& terms, const NumericState& state) {
    double result{};
    for (const auto& term : terms) {
        bool satisfies_deltas = true;
        for (const auto& delta : term.delta_momenta) {
            if (momentum_value(delta.first, state) != momentum_value(delta.second, state)) {
                satisfies_deltas = false;
                break;
            }
        }
        if (!satisfies_deltas)
            continue;

        double contribution = term.multiplicity.value<double>();
        for (const auto& op : term.operators) {
            contribution *= expectation_value(op, state);
        }
        result += contribution;
    }
    return result;
}

WickTermCollector normal_order_then_wick(const TermCollector& terms,
                                         const std::vector<WickOperatorTemplate>& templates) {
    TermCollector normal_ordered = terms;
    normal_ordered.normal_order();
    normal_ordered.clean_up();
    return wicks_theorem(normal_ordered, templates, true);
}
}  // namespace

int main() {
    const std::vector<WickOperatorTemplate> templates(
        {WickOperatorTemplate{{Num_Comparison}, Momentum(), OperatorType::Number},
         WickOperatorTemplate{{SC_Comparison}, Momentum(), OperatorType::SC}});

    const Momentum k('k');
    const Momentum p('p');
    const Operator c_k{k, Index::SpinUp, false};
    const Operator c_minus_k{-k, Index::SpinDown, false};
    const Operator c_p{p, Index::SpinUp, false};
    const Operator c_minus_p{-p, Index::SpinDown, false};
    const Operator c_k_dagger{k, Index::SpinUp, true};
    const Operator c_p_dagger{p, Index::SpinUp, true};
    const Operator c_minus_p_dagger{-p, Index::SpinDown, true};

    const TermCollector annihilation_then_creation(
        {Term(1, std::vector<Operator>{c_minus_k, c_k, c_minus_p_dagger, c_p_dagger})});
    const TermCollector interleaved({Term(1, std::vector<Operator>{c_k_dagger, c_p, c_p_dagger, c_k})});

    const auto annihilation_direct = wicks_theorem(annihilation_then_creation, templates, true);
    const auto annihilation_normal_ordered = normal_order_then_wick(annihilation_then_creation, templates);
    const auto interleaved_direct = wicks_theorem(interleaved, templates, true);
    const auto interleaved_normal_ordered = normal_order_then_wick(interleaved, templates);

    NumericState equal_momenta{
        {{'k', 1.0}, {'p', 1.0}}, {{{1.0, Index::SpinUp}, 0.2}, {{-1.0, Index::SpinDown}, 0.4}}, {{1.0, 0.3}}};
    NumericState distinct_momenta{{{'k', 1.0}, {'p', 2.0}},
                                  {{{1.0, Index::SpinUp}, 0.2}, {{2.0, Index::SpinUp}, 0.35}},
                                  {{1.0, 0.3}, {2.0, 0.25}}};

    const double annihilation_direct_value = evaluate(annihilation_direct, equal_momenta);
    const double annihilation_normal_ordered_value = evaluate(annihilation_normal_ordered, equal_momenta);
    const double interleaved_direct_value = evaluate(interleaved_direct, distinct_momenta);
    const double interleaved_normal_ordered_value = evaluate(interleaved_normal_ordered, distinct_momenta);

    std::cout << "\\begin{align*}\n"
              << "\\langle cc c^{\\dagger} c^{\\dagger} \\rangle_{\\mathrm{direct}} &= " << annihilation_direct_value
              << " \\\\\n"
              << "\\langle cc c^{\\dagger} c^{\\dagger} \\rangle_{\\mathrm{normal\\ order\\ first}} &= "
              << annihilation_normal_ordered_value << " \\\\\n"
              << "\\langle c^{\\dagger} c c^{\\dagger} c \\rangle_{\\mathrm{direct}} &= " << interleaved_direct_value
              << " \\\\\n"
              << "\\langle c^{\\dagger} c c^{\\dagger} c \\rangle_{\\mathrm{normal\\ order\\ first}} &= "
              << interleaved_normal_ordered_value << "\n\\end{align*}" << std::endl;

    constexpr double tolerance = 1e-12;
    if (std::abs(annihilation_direct_value - annihilation_normal_ordered_value) > tolerance ||
        std::abs(annihilation_direct_value - (-0.57)) > tolerance ||
        std::abs(interleaved_direct_value - interleaved_normal_ordered_value) > tolerance ||
        std::abs(interleaved_direct_value - 0.13) > tolerance) {
        return 1;
    }
    return 0;
}