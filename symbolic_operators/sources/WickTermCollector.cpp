#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/WickTermCollector.hpp>

#include <functional>
#include <memory>
#include <vector>

namespace mrock::symbolic_operators {

void WickTermCollector::clean_up(const int trivial_spin_summation_factor) {
    clean_up(std::vector<std::unique_ptr<WickSymmetry>>{}, trivial_spin_summation_factor);
}

void WickTermCollector::clear_etas() {
    for (auto it = terms.begin(); it != terms.end();) {
        bool isEta = false;
        for (const auto& op : it->operators) {
            if (op.type == OperatorType::Eta) {
                isEta = true;
                break;
            }
        }
        if (isEta) {
            it = terms.erase(it);
        } else {
            ++it;
        }
    }
}

void WickTermCollector::clean_up(
    const std::vector<std::unique_ptr<WickSymmetry>>& symmetries,
    const int trivial_spin_summation_factor) {
    for (auto& term : terms) {
        for (std::vector<Coefficient>::iterator it = term.coefficients.begin(); it != term.coefficients.end();) {
            if (it->name == "") {
                it = term.coefficients.erase(it);
            } else {
                ++it;
            }
        }
    }
    for (WickTermCollector::iterator it = terms.begin(); it != terms.end();) {
        if (!(it->resolve_deltas())) {
            it = terms.erase(it);
            continue;
        }
        it->discard_zero_momenta();
        it->rename_sums();
        it->sort();

        /*
        This function must be called before symmetries are applied!
        Sometimes <o^dagger> = <o> is a symmetry, which would transform <o^dagger> <o> into <o><o>.
        is_pauli_forbidden() transforms <o> back into the original operators and checks, whether its legal.
        That is, if <o> = <c_-k c_k>, then
        <o^dagger> <o> becomes <c_k^dagger c_-k^dagger c_-k c_k> which is finite.
        Applying the aforementioned symmetry gives
        <o><o> = <c_-k c_k c_-k c_k> which would be Pauli forbidden. */
        if (it->is_pauli_forbidden()) {
            it = terms.erase(it);
            continue;
        }

        for (const auto& symmetry : symmetries) {
            symmetry->apply_to(it->operators);
        }

        for (auto jt = it->sums.spins.begin(); jt != it->sums.spins.end();) {
            if (it->uses_index(*jt)) {
                ++jt;
            } else {
                it->multiplicity *= trivial_spin_summation_factor;
                jt = it->sums.spins.erase(jt);
            }
        }
        for (auto& coeff : it->coefficients) {
            coeff.use_custom_symmetry();
        }
        ++it;
    }

    // Setup so that we always have a structure like delta_(l,k+something)
    for (auto& term : terms) {
        for (auto& delta : term.delta_momenta) {
            assert(delta.first.momentum_list.size() == 1U);
            int l_is_at = delta.first.is_used_at('l');
            if (l_is_at == 0)
                continue;

            l_is_at = delta.second.is_used_at('l');
            if (l_is_at == -1) {
                // No l in the delta, skip the logic
                continue;
            }
            const Momentum l_mom('l', delta.second.momentum_list[l_is_at].factor);
            const Momentum remainder = delta.second - l_mom;
            delta -= remainder;
            std::swap(delta.first, delta.second);
            if (delta.first.add_PI) {
                delta.second.add_PI = !delta.second.add_PI;
                delta.first.add_PI = false;
            }
        }
    }

    for (auto& term : terms) {
        std::sort(term.operators.begin(), term.operators.end());
        std::sort(term.coefficients.begin(), term.coefficients.end());
        std::sort(term.delta_indices.begin(), term.delta_indices.end());
        std::sort(term.delta_momenta.begin(), term.delta_momenta.end());
        term.sums.sort();
    }

    combine_duplicates();

    // Sort terms
    std::sort(terms.begin(), terms.end(), std::greater<WickTerm>());
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

std::ostream& operator<<(std::ostream& os, const WickTermCollector& terms) {
    if (terms.empty()) {
        return (os << "0");
    }
    for (WickTermCollector::const_iterator it = terms.begin(); it != terms.end(); ++it) {
        if (terms.size() == 1U) {
            os << "\t";
        }
        else {
            os << "\t&";
        }
        os << *it;
        if (it != terms.end() - 1) {
            os << " \\\\";
        }
        os << "\n";
    }
    return os;
}

}  // namespace mrock::symbolic_operators
