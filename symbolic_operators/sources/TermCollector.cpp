#include <mrock/symbolic_operators/TermCollector.hpp>
#include <mrock/symbolic_operators/detail/container_helper.hpp>

#include <cstddef>
#include <string>

namespace mrock::symbolic_operators {

TermCollector commutator(const Term& left, const TermCollector& right) {
    const TermCollector buffer = {left};
    return commutator(buffer, right);
}

TermCollector commutator(const TermCollector& left, const Term& right) {
    const TermCollector buffer = {right};
    return commutator(left, buffer);
}

void TermCollector::hermitian_conjugate() {
    for (auto& t : terms) {
        t.hermitian_conjugate_inplace();
    }
}

TermCollector commutator(const Term& left, const Term& right) {
    TermCollector reciever(2);
    reciever[0] = left;
    reciever[0].multiply_from_the_right(right);
    reciever[1] = left;
    reciever[1].multiply_from_the_left(right);
    reciever[1].flip_sign();

    reciever.normal_order();
    return reciever;
}

TermCollector commutator(const TermCollector& left, const TermCollector& right) {
    TermCollector reciever = -left;
    reciever.multiply_from_the_left(right);
    append_vector(reciever, left * right);
    reciever.normal_order();
    return reciever;
}

void TermCollector::normal_order() {
    for (std::size_t t = 0U; t < terms.size();) {
    normal_order_outerLoop:
        if (t >= terms.size())
            break;
        std::size_t n = terms[t].operators.size();
        std::size_t new_n{};
        // First sort so that the bosons are upfront
        while (n > 1U) {
            new_n = 0U;
            for (std::size_t i = 1U; i < terms[t].operators.size(); ++i) {
                if (terms[t].operators[i - 1].is_fermion && !terms[t].operators[i].is_fermion) {
                    new_n = i;
                    std::swap(terms[t].operators[i - 1], terms[t].operators[i]);
                }
            }
            n = new_n;
        }

        n = terms[t].operators.size();
        new_n = 0U;
        while (n > 1U) {
            new_n = 0U;
            for (std::size_t i = 1U; i < terms[t].operators.size(); ++i) {
                if (!terms[t].operators[i - 1].is_fermion && terms[t].operators[i].is_fermion)
                    continue;
                if (!(terms[t].operators[i - 1].is_daggered) && (terms[t].operators[i].is_daggered)) {
                    bool other_deltas = false;
                    new_n = i;
                    // Swap cc^+
                    std::swap(terms[t].operators[i - 1], terms[t].operators[i]);

                    // Add a new term where cc^+ is replaced by the appropriate delta
                    Term new_term(terms[t]);
                    // flip the signs if we have fermions
                    if (terms[t].operators[i - 1].is_fermion && terms[t].operators[i].is_fermion) {
                        terms[t].flip_sign();
                    }
                    if (new_term.operators[i - 1].indices.size() != new_term.operators[i].indices.size()) {
                        throw std::invalid_argument("Operators do not have the same index count.");
                    }

                    if ((new_term.operators[i - 1].first_index() == Index::SpinUp &&
                         new_term.operators[i].first_index() == Index::SpinDown) ||
                        (new_term.operators[i - 1].first_index() == Index::SpinDown &&
                         new_term.operators[i].first_index() == Index::SpinUp)) {
                        continue;
                    } else if (new_term.operators[i - 1].first_index() != new_term.operators[i].first_index()) {
                        new_term.delta_indices.push_back(
                            make_delta(new_term.operators[i - 1].first_index(), new_term.operators[i].first_index()));
                    }
                    for (std::size_t c = 1; c < new_term.operators[i - 1].indices.size(); c++) {
                        // if the indices are not the same we emplace a delta
                        // otherwise no action is required
                        if (new_term.operators[i - 1].indices[c] != new_term.operators[i].indices[c]) {
                            other_deltas = true;
                            new_term.delta_indices.push_back(
                                make_delta(new_term.operators[i - 1].indices[c], new_term.operators[i].indices[c]));
                        }
                    }
                    if (new_term.operators[i - 1].momentum != new_term.operators[i].momentum) {
                        other_deltas = true;
                        new_term.delta_momenta.push_back(
                            make_delta(new_term.operators[i - 1].momentum, new_term.operators[i].momentum));
                    } else {
                        other_deltas = true;
                    }

                    new_term.operators.erase(new_term.operators.begin() + i - 1, new_term.operators.begin() + i + 1);

                    if (other_deltas)
                        terms.push_back(new_term);
                } else if (terms[t].operators[i - 1] == terms[t].operators[i]) {
                    if (terms[t].operators[i - 1].is_fermion) {
                        // two identical fermion operators = 0
                        terms.erase(terms.begin() + t);
                        goto normal_order_outerLoop;
                    }
                }
            }
            n = new_n;
        }
        ++t;
    }

    for (auto& term : terms) {
        term.sort_operators_by_indices();
    }

    combine_duplicates();
}

TermCollector operator-(TermCollector terms) {
    for (auto& term : terms) {
        term.multiplicity *= -1;
    }
    return terms;
}

TermCollector& TermCollector::operator+=(const TermCollector& other) {
    if (other.empty()) {
        return (*this);
    }
    append_vector(this->terms, other.terms);
    return (*this);
}

TermCollector& TermCollector::operator-=(const TermCollector& other) {
    if (other.empty()) {
        return (*this);
    }
    // we effectively compute -(-(*this) + other)
    // Doing so forgoes creating the temporary of -other
    for (auto& term : this->terms) {
        term.multiplicity *= -1;
    }
    append_vector(this->terms, other.terms);
    for (auto& term : this->terms) {
        term.multiplicity *= -1;
    }
    return (*this);
}

TermCollector& TermCollector::operator*=(const TermCollector& other) {
    if (other.empty()) {
        this->clear();
        return (*this);
    };
    const std::size_t n = this->size();
    duplicate_n_inplace(this->terms, other.size() - 1U);
    for (std::size_t l = 0U; l < n; ++l) {
        for (std::size_t r = 0U; r < other.size(); ++r) {
            this->terms[l + r * n] *= other.terms[r];
        }
    }
    return (*this);
}

std::ostream& operator<<(std::ostream& os, const TermCollector& terms) {
    for (TermCollector::const_iterator it = terms.begin(); it != terms.end(); ++it) {
        os << "\t&" << *it;
        if (it != terms.end() - 1) {
            os << " \\\\";
        }
        os << "\n";
    }
    return os;
}

void TermCollector::clean_up() {
    for (auto it = terms.begin(); it != terms.end();) {
        if (!(it->resolve_deltas())) {
            it = terms.erase(it);
            continue;
        }
        it->discard_zero_momenta();
        it->rename_sums();
        it->structure();
        ++it;
    }

    // Setup so that we always have a structure like delta_(l,k+something)
    for (auto& term : terms) {
        for (auto& delta : term.delta_momenta) {
            assert(delta.first.momentum_list.size() == 1U);
            int l_is = delta.first.is_used_at('l');
            if (l_is == 0)
                continue;

            l_is = delta.second.is_used_at('l');
            if (l_is == -1) {
                std::cout << term << std::endl;
                throw;
            }
            const Momentum l_mom('l', delta.second.momentum_list[l_is].factor);
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

    // Sort terms
    for (std::size_t i = 0; i < terms.size(); i++) {
        for (std::size_t j = i + 1; j < terms.size(); j++) {
            if (terms[i].sums.momenta.empty() && terms[j].sums.momenta.size() > 0) {
                std::swap(terms[i], terms[j]);
            }
            if (terms[i].sums.momenta.size() > 0 && terms[j].sums.momenta.size() > 0) {
                if (terms[i].sums.momenta.size() < terms[j].sums.momenta.size()) {
                    std::swap(terms[i], terms[j]);
                } else if (terms[i].sums.momenta.size() == terms[j].sums.momenta.size()) {
                    if (terms[i].coefficients.size() > 0) {
                        if (terms[j].coefficients[0].name < terms[i].coefficients[0].name) {
                            std::swap(terms[i], terms[j]);
                        }
                    }
                }
            } else if (terms[i].sums.momenta.empty() && terms[j].sums.momenta.empty()) {
                if (terms[i].coefficients.size() > 0) {
                    if (terms[j].coefficients[0].name < terms[i].coefficients[0].name) {
                        std::swap(terms[i], terms[j]);
                    }
                }
            }
        }
    }

#ifndef NDEBUG
    for (const auto& term : terms) {
        assert(term.is_normal_ordered());
    }
#endif
}

TermCollector& TermCollector::multiply_from_the_left(const TermCollector& other) {
    if (other.empty()) {
        this->clear();
        return (*this);
    };
    const std::size_t n = this->size();
    duplicate_n_inplace(this->terms, other.size() - 1U);
    for (std::size_t l = 0U; l < n; ++l) {
        for (std::size_t r = 0U; r < other.size(); ++r) {
            this->terms[l + r * n].multiply_from_the_left(other.terms[r]);
        }
    }
    return (*this);
}

TermCollector& TermCollector::multiply_from_the_right(const TermCollector& other) {
    return ((*this) *= other);
}

std::string TermCollector::to_string_without_prefactor() const {
    std::string ret = "";
    for (std::size_t i = 0U; i < terms.size(); i++) {
        if (terms[i].multiplicity < 0) {
            ret += "-";
        } else if (i > 0) {
            ret += "+";
        }
        ret += terms[i].to_string_without_prefactor();
    }
    return ret;
}

}  // namespace mrock::symbolic_operators
