#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/WickTermCollector.hpp>

#include <memory>
#include <vector>

namespace mrock::symbolic_operators {

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
    const std::vector<std::unique_ptr<WickSymmetry>>& symmetries /*= std::vector<std::unique_ptr<WickSymmetry>>{}*/) {
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
        Applyng the aforementioned symmetry gives
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
                // We are assuming there are only spin indizes here (spin 1/2)
                // If another kind of index arises I have to readress this section.
                it->multiplicity *= 2;
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

    combine_duplicates();

    auto predicate = [](const WickTerm& left, const WickTerm& right) -> bool {
        if (left.delta_momenta.empty() && right.delta_momenta.size() > 0) {
            return true;
        } else if (left.delta_momenta.size() > 0 && right.delta_momenta.size() > 0) {
            if (left.delta_momenta.size() < right.delta_momenta.size()) {
                return true;
            } else if (left.delta_momenta.size() == right.delta_momenta.size()) {
                if (left.delta_momenta[0].second.add_PI && !(right.delta_momenta[0].second.add_PI)) {
                    return true;
                } else if (!left.coefficients.empty() && right.coefficients[0].name < left.coefficients[0].name) {
                    return true;
                } else if ((!left.coefficients.empty() && right.coefficients[0].name == left.coefficients[0].name) ||
                           left.coefficients.empty()) {
                    if (!left.operators.empty() && right.operators.empty()) {
                        return true;
                    } else if ((!left.operators.empty() && !right.operators.empty()) &&
                               left.operators.front().type < right.operators.front().type) {
                        return true;
                    }
                }
            }
        } else if (left.delta_momenta.empty() && right.delta_momenta.empty()) {
            if (left.coefficients.size() < right.coefficients.size()) {
                return true;
            } else if (!left.coefficients.empty() && !right.coefficients.empty()) {
                if (right.coefficients[0].name < left.coefficients[0].name) {
                    return true;
                } else if (right.coefficients[0].name == left.coefficients[0].name) {
                    if (left.operators.size() > right.operators.size()) {
                        return true;
                    } else if ((!left.operators.empty() && !right.operators.empty()) &&
                               left.operators.front().type < right.operators.front().type) {
                        return true;
                    }
                }
            }
        }
        return false;
    };

    // Sort terms
    for (std::size_t i = 0U; i < terms.size(); i++) {
        for (std::size_t j = i + 1U; j < terms.size(); j++) {
            if (predicate(terms[i], terms[j]))
                std::swap(terms[i], terms[j]);
        }
    }
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
    for (WickTermCollector::const_iterator it = terms.begin(); it != terms.end(); ++it) {
        os << "\t&" << *it;
        if (it != terms.end() - 1) {
            os << " \\\\";
        }
        os << "\n";
    }
    return os;
}

}  // namespace mrock::symbolic_operators
