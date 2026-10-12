#include <mrock/symbolic_operators/experimental/WickOrderedTerm.hpp>

#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/symbolic_operators/detail/container_helper.hpp>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <variant>

namespace mrock::symbolic_operators::experimental {
void WickOrderedTerm::for_each_momentum_except_deltas(const std::function<void(Momentum&)>& f) 
{
    AbstractTerm<WickOperator>::for_each_momentum_except_deltas(f);
    for (auto& op : wick_expression) {
        f(op.momentum);
    }
}

void WickOrderedTerm::for_each_index_except_deltas(const std::function<void(IndexWrapper&)>& f) 
{
    AbstractTerm<WickOperator>::for_each_index_except_deltas(f);
    for (auto& op : wick_expression) {
        f(op.indices);
    }
}

WickOrderedTerm::WickOrderedTerm(const Term& base) 
    : AbstractTerm<WickOperator>(base.multiplicity,
                                 base.coefficients,
                                 base.sums,
                                 base.delta_momenta,
                                 base.delta_indices,
                                 std::vector<WickOperator>()),
      wick_expression()
{ }

WickOrderedTerm::WickOrderedTerm(const WickTerm& base) 
    : AbstractTerm<WickOperator>(base.multiplicity,
                                 base.coefficients,
                                 base.sums,
                                 base.delta_momenta,
                                 base.delta_indices,
                                 base.operators),
      wick_expression()
{ }

WickOrderedTerm WickOrderedTerm::from_wick_term_skip_wick_expression(const WickOrderedTerm& base)
{
    WickOrderedTerm result;
    result.multiplicity = base.multiplicity;
    result.coefficients = base.coefficients;
    result.sums = base.sums;
    result.delta_momenta = base.delta_momenta;
    result.delta_indices = base.delta_indices;
    result.operators = base.operators;

    return result;
}

Term WickOrderedTerm::to_term(const std::string& base) const {
    Term result(multiplicity, coefficients, sums, wick_expression.operators);
    result.delta_indices = delta_indices;
    result.delta_momenta = delta_momenta;

    for (const auto& expec : this->operators) {
        result.coefficients.push_back(expec.to_coefficient(base));
    }

    return result;
}

void WickOrderedTerm::include_template_result(const TemplateResult::SingleResult& result) {
    this->delta_indices.insert(this->delta_indices.begin(), result.index_deltas.begin(), result.index_deltas.end());
    this->operators.push_back(result.op);
    this->multiplicity *= result.factor;
}

bool WickOrderedTerm::resolve_deltas()
{
    if (!resolve_momentum_deltas())
        return false;
    if (!resolve_index_deltas())
        return false;

    return true;
}

void WickOrderedTerm::discard_zero_momenta() 
{
    AbstractTerm<WickOperator>::discard_zero_momenta();
    for (auto& op : wick_expression) {
        op.momentum.remove_zeros();
    }
}

void WickOrderedTerm::sort()
{
    std::sort(operators.begin(), operators.end());
    std::sort(sums.momenta.begin(), sums.momenta.end());
    std::sort(sums.spins.begin(), sums.spins.end());
}

bool WickOrderedTerm::is_identity() const noexcept {
    return wick_expression.empty();
}

bool WickOrderedTerm::is_bilinear() const noexcept {
    return wick_expression.size() == 2U;
}

bool WickOrderedTerm::is_quartic() const noexcept {
    return wick_expression.size() == 4U;
}

WickOrderedTerm& WickOrderedTerm::hermitian_conjugate_inplace() noexcept
{
    wick_expression.hermitian_conjugate_inplace();
    return *this;
}

WickOrderedTerm WickOrderedTerm::hermitian_conjugate() const noexcept
{
    WickOrderedTerm copy(*this);
    copy.hermitian_conjugate_inplace();
    return copy;
}


std::ostream& operator<<(std::ostream& os, const WickOrderedTerm& term)
{
    if (term.multiplicity > 0) {
        os << "+";
    }
    os << term.multiplicity << " ";
    os << term.sums;
    os << term.coefficients << " ";
    for (const auto& delta : term.delta_momenta) {
        os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
    }
    for (const auto& delta : term.delta_indices) {
        os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
    }
    for (const auto& op : term.operators) {
        os << op << " ";
    }
    os << term.wick_expression;

    return os;
}

std::strong_ordering WickOrderedTerm::operator<=>(const WickOrderedTerm& other) const {
    if (auto cmp = AbstractTerm<WickOperator>::operator<=>(other); cmp != 0)
        return cmp;
    return wick_expression.operators <=> other.wick_expression.operators;
}

bool WickOrderedTerm::operator==(const WickOrderedTerm& other) const {
    return (*this <=> other) == 0;
}

}