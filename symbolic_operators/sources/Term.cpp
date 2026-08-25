#include <mrock/symbolic_operators/AbstractTerm.hpp>
#include <mrock/symbolic_operators/Coefficient.hpp>
#include <mrock/symbolic_operators/Fractional.hpp>
#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <mrock/symbolic_operators/KroneckerDelta.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/MomentumList.hpp>
#include <mrock/symbolic_operators/MomentumSymbol.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/SumContainer.hpp>
#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/detail/container_helper.hpp>

#include <cassert>
#include <compare>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace mrock::symbolic_operators {
Term::Term(IntFractional _multiplicity,
           std::vector<Coefficient> _coefficients,
           const SumContainer& _sums,
           const std::vector<Operator>& _operators)
    : AbstractTerm<Operator>(_multiplicity, _coefficients, _sums, _operators) {}

Term::Term(IntFractional _multiplicity,
           Coefficient _coefficient,
           const SumContainer& _sums,
           const std::vector<Operator>& _operators)
    : AbstractTerm<Operator>(_multiplicity, _coefficient, _sums, _operators) {}

Term::Term(IntFractional _multiplicity,
           Coefficient _coefficient,
           const MomentumSum& _sum_momenta,
           const std::vector<Operator>& _operators)
    : AbstractTerm<Operator>(_multiplicity, _coefficient, SumContainer{_sum_momenta, {}}, _operators) {}

Term::Term(IntFractional _multiplicity,
           Coefficient _coefficient,
           const IndexSum& _sum_indizes,
           const std::vector<Operator>& _operators)
    : AbstractTerm<Operator>(_multiplicity, _coefficient, SumContainer{{}, _sum_indizes}, _operators) {}

Term::Term(IntFractional _multiplicity, Coefficient _coefficient, const std::vector<Operator>& _operators)
    : AbstractTerm<Operator>(_multiplicity, _coefficient, _operators) {}

Term::Term(IntFractional _multiplicity, const SumContainer& _sums, const std::vector<Operator>& _operators)
    : AbstractTerm<Operator>(_multiplicity, std::vector<Coefficient>(), _sums, _operators) {}

Term::Term(IntFractional _multiplicity, const MomentumSum& _sum_momenta, const std::vector<Operator>& _operators)
    : AbstractTerm<Operator>(_multiplicity, std::vector<Coefficient>(), SumContainer{_sum_momenta, {}}, _operators) {}

Term::Term(IntFractional _multiplicity, const IndexSum& _sum_indizes, const std::vector<Operator>& _operators)
    : AbstractTerm<Operator>(_multiplicity, std::vector<Coefficient>(), SumContainer{{}, _sum_indizes}, _operators) {}

Term::Term(IntFractional _multiplicity, const std::vector<Operator>& _operators)
    : AbstractTerm<Operator>(_multiplicity, _operators) {}

void Term::print() const {
    std::cout << *this << std::endl;
}

bool Term::resolve_deltas() {
    if (!resolve_momentum_deltas())
        return false;
    if (!resolve_index_deltas())
        return false;

    // Check for the Pauli principle
    for (auto it = operators.begin(); it != operators.end(); ++it) {
        for (auto jt = it + 1; jt != operators.end(); ++jt) {
            if (*it == *jt)
                return false;
            if (it->is_daggered != jt->is_daggered)
                break;
        }
    }

    return true;
}

void Term::sort() {
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
                flip_sign();
            }
        }
    }
    std::size_t new_n;
    std::size_t n = operators.size();
    while (n > 1U) {
        new_n = 0U;
        for (std::size_t i = 1U; i < n; ++i) {
            if (operators[i].is_daggered != operators[i - 1].is_daggered)
                continue;
            if (operators[i].is_fermion != operators[i - 1].is_fermion)
                continue;
            const Index l_idx = operators[i - 1].first_index();
            const Index r_idx = operators[i].first_index();
            if (operators[i].is_daggered) {
                // c^+ c^+
                if (r_idx == Index::SpinUp && l_idx != Index::SpinUp) {
                    perform_operator_swap(operators[i], operators[i - 1]);
                    new_n = i;
                } else if (l_idx == Index::SpinDown && r_idx != Index::SpinDown) {
                    perform_operator_swap(operators[i], operators[i - 1]);
                    new_n = i;
                } else if (l_idx > Index::SpinDown && r_idx > Index::SpinDown && l_idx > r_idx) {
                    perform_operator_swap(operators[i], operators[i - 1]);
                    new_n = i;
                }
            } else {
                // c c
                if (r_idx == Index::SpinDown && l_idx != Index::SpinDown) {
                    perform_operator_swap(operators[i], operators[i - 1]);
                    new_n = i;
                } else if (l_idx == Index::SpinUp && r_idx != Index::SpinUp) {
                    perform_operator_swap(operators[i], operators[i - 1]);
                    new_n = i;
                } else if (l_idx > Index::SpinDown && r_idx > Index::SpinDown && l_idx < r_idx) {
                    perform_operator_swap(operators[i], operators[i - 1]);
                    new_n = i;
                }
            }
        }
        n = new_n;
    }

    n = operators.size();
    while (n > 1U) {
        new_n = 0U;
        for (std::size_t i = 1U; i < n; ++i) {
            if (operators[i].is_daggered != operators[i - 1].is_daggered)
                continue;
            if (operators[i].first_index() != operators[i - 1].first_index())
                continue;

            if (momentum_order(operators[i - 1].momentum, operators[i].momentum)) {
                perform_operator_swap(operators[i], operators[i - 1]);
                new_n = i;
            }
        }
        n = new_n;
    }

    // Sort the occurring coefficients in alphabetical order
    std::sort(coefficients.begin(), coefficients.end(), [](const Coefficient& a, const Coefficient& b) {
        if (a.name == b.name) {
            return (a.is_daggered && (!b.is_daggered));
        }
        return a.name < b.name;
    });

    // check whether we can swap the sign of each momentum in the coefficients
    // 26.04.2024, I have no idea what I did here, nor do I know why I did what I did
    for (const auto& coeff : coefficients) {
        if (!(coeff.inversion_symmetry))
            return;
        if (std::any_of(coeff.momenta.begin(), coeff.momenta.end(),
                        [](Momentum const& momentum) { return momentum.momentum_list.size() > 1U; }))
            return;
    }

    for (const auto& sum_mom : sums.momenta) {
        bool first_occurance = true;
        for (auto& op : operators) {
            int i = op.momentum.is_used_at(sum_mom);
            if (i > -1) {
                if (first_occurance) {
                    if (op.momentum.momentum_list[i].factor < 0) {
                        first_occurance = false;
                    } else {
                        break;
                    }
                }
                op.momentum.momentum_list[i].factor *= -1;
            }
        }
    }
}

bool Term::is_equal(const Term& other) const {
    if (this->coefficients != other.coefficients)
        return false;
    if (this->sums != other.sums)
        return false;
    if (this->delta_indizes != other.delta_indizes)
        return false;
    if (this->delta_momenta != other.delta_momenta)
        return false;
    if (this->operators != other.operators)
        return false;
    return true;
}

bool Term::is_normal_ordered() const {
    for (std::size_t i = 1U; i < operators.size(); ++i) {
        if (operators[i - 1].is_fermion == operators[i].is_fermion) {
            if (!operators[i - 1].is_daggered && operators[i].is_daggered) {
                return false;
            }
        }
    }
    return true;
}

std::string Term::to_string_without_prefactor() const {
    std::ostringstream os;
    if (!this->sums.spins.empty()) {
        os << "\\sum_{ ";
        for (const auto& index : this->sums.spins) {
            os << index << " ";
        }
        os << "}";
    }
    if (!this->sums.momenta.empty()) {
        os << "\\sum_{ ";
        for (const auto& momentum : this->sums.momenta) {
            os << momentum << " ";
        }
        os << "}";
    }
    os << this->coefficients << " ";
    for (const auto& delta : delta_momenta) {
        os << delta;
    }
    for (const auto& delta : delta_indizes) {
        os << delta;
    }

    if (this->is_identity()) {
        os << " \\hat{1} ";
        return os.str();
    }
    for (const auto& op : this->operators) {
        os << op << " ";
    }
    return os.str();
}

Term& Term::hermitian_conjugate_inplace() {
    std::reverse(this->operators.begin(), this->operators.end());
    for (auto& op : this->operators) {
        op.hermitian_conjugate_inplace();
    }
    for (auto& coeff : this->coefficients) {
        coeff.hermitian_conjugate_inplace();
    }
    return *this;
}

Term Term::hermitian_conjugate() const {
    Term copy(*this);
    copy.hermitian_conjugate_inplace();
    return copy;
}

void Term::rename_indizes(Index what, Index to) {
    if (what == to)
        return;
    for (auto& index_sum : sums.spins) {
        if (index_sum == to) {
            throw std::invalid_argument("You are replacing an index sum with an index that already exists!");
        }
        if (index_sum == what) {
            index_sum = to;
        }
    }
    for (auto& coeff : coefficients) {
        for (auto& index : coeff.indizes) {
            if (index == what) {
                index = to;
            }
        }
    }
    for (auto& op : operators) {
        for (auto& index : op.indizes) {
            if (index == what) {
                index = to;
            }
        }
    }
    for (auto& delta : delta_indizes) {
        if (delta.first == what) {
            delta.first = to;
        }
        if (delta.second == what) {
            delta.second = to;
        }
    }
}

void Term::swap_momenta(const MomentumSymbol::name_type a, const MomentumSymbol::name_type b) {
    this->rename_momenta(a, PLACEHOLDER_SYMBOL);
    this->rename_momenta(b, a);
    this->rename_momenta(PLACEHOLDER_SYMBOL, b);
}

Term& Term::multiply_from_the_left(const Term& other)
{
    this->multiplicity *= other.multiplicity;

    this->rename_duplicate_sums(&other);
    this->sums.append(other.sums);

    append_vector(this->coefficients, other.coefficients);
    append_vector(this->delta_momenta, other.delta_momenta);
    append_vector(this->delta_indizes, other.delta_indizes);

    // IMPORTANT: prepend instead of append!
    prepend_vector(this->operators, other.operators);

    return *this;
}

Term& Term::multiply_from_the_right(const Term& other)
{
    return ((*this) *= other);
}

Term& Term::operator*=(const Term& rhs) {
    this->rename_duplicate_sums(&rhs);
    this->sums.append(rhs.sums);

    this->multiplicity *= rhs.multiplicity;
    append_vector(this->coefficients, rhs.coefficients);
    append_vector(this->delta_momenta, rhs.delta_momenta);
    append_vector(this->delta_indizes, rhs.delta_indizes);
    append_vector(this->operators, rhs.operators);

    return *this;
}

std::ostream& operator<<(std::ostream& os, const Term& term) {
    if (term.multiplicity > 0) {
        os << "+";
    }
    os << term.multiplicity << " ";
    os << term.sums;
    os << term.coefficients << " ";
    for (const auto& delta : term.delta_momenta) {
        os << delta;
    }
    for (const auto& delta : term.delta_indizes) {
        os << delta;
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
}  // namespace mrock::symbolic_operators