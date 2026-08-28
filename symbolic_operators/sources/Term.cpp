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

void Term::sort_operators_by_indizes() {
    auto general_swap_predicate = [this](const std::size_t& i, const std::size_t& index_pos) {
        if (operators[i - 1].is_daggered) {
            if (operators[i - 1].indizes.front() == Index::SpinUp)
                return false;
            if (operators[i].indizes.front() == Index::SpinDown)
                return false;

            if (operators[i - 1].indizes.front() == Index::SpinDown || operators[i].indizes.front() == Index::SpinUp) {
                return true;
            }

            if (operators[i - 1].indizes[index_pos] > operators[i].indizes[index_pos]) {
                return true;
            }
        } else {
            // The comparison operator and SpinUp <-> SpinDown are the only changes
            if (operators[i - 1].indizes.front() == Index::SpinDown)
                return false;
            if (operators[i].indizes.front() == Index::SpinUp)
                return false;

            if (operators[i - 1].indizes.front() == Index::SpinUp || operators[i].indizes.front() == Index::SpinDown) {
                return true;
            }

            if (operators[i - 1].indizes[index_pos] < operators[i].indizes[index_pos]) {
                return true;
            }
        }
        return false;
    };

    std::size_t n = operators.size();
    std::size_t new_n{};
    // First sort so that the spins indizes are always ordered the same way
    // Without destroying the previously achieved normal order
    while (n > 1U) {
        new_n = 0U;
        for (std::size_t i = 1U; i < operators.size(); ++i) {
            if (operators[i - 1].is_daggered != operators[i].is_daggered) {
                continue;
            }
            if (operators[i - 1].is_fermion != operators[i].is_fermion) {
                continue;
            }
            std::size_t j = 0U;
            while (j < operators[i - 1].indizes.size() && j < operators[i].indizes.size()) {
                if (general_swap_predicate(i, j)) {
                    perform_operator_swap(i - 1, i);
                    new_n = i;
                    break;
                }
                ++j;
            }
        }
        n = new_n;
    }
}

void Term::structure_momentum_dependencies() {
    // This function does not do anything to an identity term
    if (operators.empty())
        return;

    auto most_similar_daggered_operator = [this](const IndexWrapper& indizes) {
        int best_match_count = 0;
        std::vector<Operator>::const_iterator best_it = operators.end();

        for (std::vector<Operator>::const_iterator it = operators.begin(); it != operators.end() && it->is_daggered;
             ++it) {
            // Penalize size mismatch
            int current_match_count = static_cast<int>(indizes.size()) - static_cast<int>(it->indizes.size());
            for (std::size_t i = 0U; i < it->indizes.size() && i < indizes.size(); ++i) {
                if (it->indizes[i] == indizes[i]) {
                    ++current_match_count;
                }
            }
            if (current_match_count > best_match_count) {
                best_match_count = current_match_count;
                best_it = it;
            }
        }
        return best_it;
    };

    MomentumSum::const_iterator sum_it = sums.momenta.begin();
    std::vector<MomentumSymbol::name_type> do_not_touch;
    do_not_touch.reserve(sums.momenta.size());

    for (auto op_it = operators.begin(); op_it < operators.end() && sum_it != sums.momenta.end(); ++op_it) {
        if (!op_it->is_daggered)
            break;

        try {
            const MomentumSymbol::name_type target = *sum_it;
            redistribute_momenta(op_it->momentum, target, do_not_touch);
            do_not_touch.push_back(target);
            ++sum_it;
        } catch (redistribution_error& e) {
        }
    }

    // Iterate backwards
    for (auto op_it = operators.rbegin(); op_it < operators.rend() && sum_it != sums.momenta.end(); ++op_it) {
        if (op_it->is_daggered)
            break;

        auto best_it = most_similar_daggered_operator(op_it->indizes);
        if (best_it == operators.end())
            continue;
        try {
            const MomentumSymbol::name_type target = *sum_it;
            redistribute_momenta(op_it->momentum, target, do_not_touch);
            do_not_touch.push_back(target);
            ++sum_it;
            transform_momentum_sum(target, best_it->momentum + Momentum(PLACEHOLDER_SYMBOL), PLACEHOLDER_SYMBOL);
            rename_momenta(PLACEHOLDER_SYMBOL, target);
        } catch (redistribution_error& e) {
        }
    }

    for (auto coeff_it = coefficients.rbegin(); coeff_it != coefficients.rend(); ++coeff_it) {
        for (auto mom_it = coeff_it->momenta.rbegin(); mom_it != coeff_it->momenta.rend(); ++mom_it) {
            if (sum_it == sums.momenta.end())
                return;
            try {
                const MomentumSymbol::name_type target = *sum_it;
                redistribute_momenta(*mom_it, target, do_not_touch);
                do_not_touch.push_back(target);
                ++sum_it;
            } catch (redistribution_error& e) {
            }
        }
    }
}

void Term::structure() {
    sort_operators_by_indizes();

    std::size_t new_n;
    std::size_t n = operators.size();
    while (n > 1U) {
        new_n = 0U;
        for (std::size_t i = 1U; i < n; ++i) {
            if (operators[i].is_daggered != operators[i - 1].is_daggered)
                continue;
            if (operators[i].first_index() != operators[i - 1].first_index())
                continue;

            if (momentum_order(operators[i - 1].momentum, operators[i].momentum)) {
                perform_operator_swap(i, i - 1);
                new_n = i;
            }
        }
        n = new_n;
    }

    structure_momentum_dependencies();

    // Sort the occurring coefficients in alphabetical order
    std::sort(coefficients.begin(), coefficients.end(), [](const Coefficient& a, const Coefficient& b) {
        if (a.name == b.name) {
            return (a.is_daggered && (!b.is_daggered));
        }
        return a.name < b.name;
    });

    for (auto& coeff : coefficients) {
        if (coeff.momenta.size() == 3U) {
            if (coeff.momenta[0].empty() || coeff.momenta[1].empty())
                continue;
            if (coeff.momenta[0].front().name > coeff.momenta[1].front().name) {
                coeff.use_symmetric_interaction_exchange();
            }
        }
    }

    for (auto& coeff : coefficients) {
        if (coeff.momenta.empty())
            continue;
        if (coeff.momenta.back().size() != 1U || !coeff.momenta.back().first_momentum_is_negative())
            continue;

        const auto target_name = coeff.momenta.back().front().name;
        if (!sums.momenta.is_summed_over(target_name))
            continue;

        if (std::none_of(operators.begin(), operators.end(),
                         [&target_name](const Operator& op) { return op.momentum.uses(target_name); })) {
            // If none of the operators uses the momentum, we can change it without
            // disrupting the previously achieved operator structure
            invert_momentum_sum(target_name);
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
            if (coeff.Q_changes_sign && momentum.add_PI) {
                momentum.add_PI = false;
                flip_sign();
            }
        }
    }

    // check whether we can swap the sign of each momentum in the coefficients
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

Term& Term::multiply_from_the_left(const Term& other) {
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

Term& Term::multiply_from_the_right(const Term& other) {
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
}  // namespace mrock::symbolic_operators