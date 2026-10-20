#include <mrock/symbolic_operators/TermCollector.hpp>
#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/detail/OperatorOrder.hpp>
#include <mrock/symbolic_operators/detail/container_helper.hpp>
#include <mrock/symbolic_operators/experimental/WickOrderedCollector.hpp>

#include <algorithm>
#include <list>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

namespace mrock::symbolic_operators::experimental {

const std::string WICK_OP_TO_STRING_TEMP_BASE = "WICK_OP_TO_STRING_TEMP_BASE";

// Computes n!!
template <typename IntegerType, typename RealType = std::size_t>
constexpr RealType double_factorial(const IntegerType n) {
    RealType result{1};
    for (std::make_signed_t<IntegerType> i = n; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

using IndexVector = std::vector<std::size_t>;
using Pair = std::pair<std::size_t, std::size_t>;
using Pairing = std::vector<Pair>;
using PairingVector = std::vector<Pairing>;

void generate_pairings(const IndexVector& remaining, std::size_t n, Pairing& current, PairingVector& result) {
    // We have generated enough pairs
    if (current.size() == n) {
        result.push_back(current);
        return;
    }

    // Not enough elements left to make another pair
    if (remaining.size() < 2) {
        return;
    }

    // Always consider the smallest remaining index first
    std::size_t first = remaining[0];

    // Option 1: leave it unused
    IndexVector skipped(remaining.begin() + 1, remaining.end());
    generate_pairings(skipped, n, current, result);

    // Option 2: pair it with another remaining index
    for (std::size_t i = 1; i < remaining.size(); ++i) {
        std::size_t second = remaining[i];

        IndexVector newRemaining;
        newRemaining.reserve(remaining.size() - 2);

        for (std::size_t j = 1; j < remaining.size(); ++j) {
            if (j != i) {
                newRemaining.push_back(remaining[j]);
            }
        }

        current.push_back({first, second});

        generate_pairings(newRemaining, n, current, result);

        current.pop_back();
    }
}

PairingVector get_pairings(std::size_t n_indices, std::size_t n_pairings) {
    PairingVector result;
    Pairing current;

    IndexVector indices(n_indices);
    std::iota(indices.begin(), indices.end(), std::size_t{});

    generate_pairings(indices, n_pairings, current, result);

    return result;
}

int pair_parity(Pairing pairs) {
    int parity = 0;

    // Normalize ordering
    for (auto& p : pairs) {
        if (p.first > p.second)
            std::swap(p.first, p.second);

        // internal swaps needed to make this pair adjacent
        parity ^= ((p.second - p.first - 1) & 1);
    }

    // crossings
    for (size_t i = 0; i < pairs.size(); ++i) {
        auto [a, b] = pairs[i];

        for (size_t j = i + 1; j < pairs.size(); ++j) {
            auto [c, d] = pairs[j];

            if ((a < c && c < b && b < d) || (c < a && a < d && d < b)) {
                parity ^= 1;
            }
        }
    }

    return parity ? -1 : +1;
}

template <typename T>
std::vector<T> remove_elements(const std::vector<T>& vec, const Pairing& pairing) {
    std::vector<bool> remove(vec.size(), false);

    for (auto [a, b] : pairing) {
        remove[a] = true;
        remove[b] = true;
    }

    std::vector<T> result;
    result.reserve(vec.size());

    for (size_t i = 0; i < vec.size(); i++) {
        if (!remove[i]) {
            result.push_back(vec[i]);
        }
    }

    return result;
}

void wick_decompose_internal(const std::vector<Operator>& operators,
                             const std::vector<WickOperatorTemplate>& templates,
                             WickOrderedTerm&& attach_to,
                             WickOrderedCollector& result) {
    WickOrderedTerm contraction_base = attach_to;
    result.emplace_back(std::move(attach_to));
    result.back().wick_expression.operators = operators;

    struct ContractionOption {
        TemplateResult::SingleResult result;
        KroneckerDelta<Momentum> momentum_delta;
    };

    for (std::size_t n_contractions = 1U; n_contractions <= operators.size() / 2; ++n_contractions) {
        const PairingVector pairings = get_pairings(operators.size(), n_contractions);

        for (const auto& pairing : pairings) {
            std::vector<std::vector<ContractionOption>> options_by_pair;
            options_by_pair.reserve(pairing.size());
            bool pairing_is_valid = true;

            for (const auto& pair : pairing) {
                std::vector<ContractionOption> options;
                for (const auto& wick_template : templates) {
                    TemplateResult template_result =
                        wick_template.create_from_operators(operators[pair.first], operators[pair.second]);
                    template_result.clean_up();
                    for (const auto& single_result : template_result.results) {
                        options.push_back({single_result, template_result.momentum_delta});
                    }
                }
                if (options.empty()) {
                    pairing_is_valid = false;
                    break;
                }
                options_by_pair.push_back(std::move(options));
            }
            if (!pairing_is_valid)
                continue;

            WickOrderedTerm base = WickOrderedTerm::from_wick_term_skip_wick_expression(contraction_base);
            base.multiplicity *= pair_parity(pairing);
            auto append_contraction_options = [&](auto&& self, std::size_t pair_index,
                                                  WickOrderedTerm partial) -> void {
                if (pair_index == options_by_pair.size()) {
                    partial.wick_expression.operators = remove_elements(operators, pairing);
                    result.emplace_back(std::move(partial));
                    return;
                }

                for (const auto& option : options_by_pair[pair_index]) {
                    WickOrderedTerm branch = partial;
                    branch.include_template_result(option.result);
                    if (!option.momentum_delta.isOne())
                        branch.delta_momenta.push_back(option.momentum_delta);
                    self(self, pair_index + 1U, std::move(branch));
                }
            };
            append_contraction_options(append_contraction_options, 0U, std::move(base));
        }
    }
}

WickOrderedCollector wick_decompose(const Term& term,
                                    const std::vector<WickOperatorTemplate>& templates,
                                    const bool clean_up) {
    const std::size_t estimated_size = double_factorial(term.get_operators().size()) * templates.size();
    WickOrderedCollector result;
    result.reserve(estimated_size);

    WickOrderedTerm initial(term);
    wick_decompose_internal(term.get_operators(), templates, std::move(initial), result);

    if (clean_up) {
        clean_wick_ordered_terms(result);
    }

    return result;
}

WickOrderedCollector wick_decompose(const TermCollector& terms,
                                    const std::vector<WickOperatorTemplate>& templates,
                                    const bool clean_up) {
    WickOrderedCollector result;
    std::size_t estimated_size =
        std::accumulate(terms.begin(), terms.end(), std::size_t{}, [&templates](std::size_t current, const auto& term) {
            return current + double_factorial(term.get_operators().size()) * templates.size();
        });
    result.reserve(estimated_size);

    for (const auto& term : terms) {
        WickOrderedTerm initial(term);
        wick_decompose_internal(term.get_operators(), templates, std::move(initial), result);
    }

    if (clean_up) {
        clean_wick_ordered_terms(result);
    }
    return result;
}

WickOrderedCollector multiply_impl(const WickOrderedTerm& lhs,
                                   const WickOrderedTerm& rhs,
                                   const std::vector<WickOperatorTemplate>& templates) {
    WickOrderedCollector result;

    Term term_l = lhs.to_term(WICK_OP_TO_STRING_TEMP_BASE);
    Term term_r = rhs.to_term(WICK_OP_TO_STRING_TEMP_BASE);

    for (auto& op : term_l.operators) {
        op.indices.push_back(Index::WickA);
    }
    for (auto& op : term_r.operators) {
        op.indices.push_back(Index::WickB);
    }
    term_l *= term_r;

    result = wick_decompose(term_l, templates);

    const auto is_tmp_index = [](Index index) { return (index == Index::WickA || index == Index::WickB); };

    for (auto& term : result) {
        for (auto& op : term.operators) {
            std::erase_if(op.indices.indices, is_tmp_index);
        }
        for (auto& op : term.wick_expression) {
            std::erase_if(op.indices.indices, is_tmp_index);
        }
        for (auto it = term.coefficients.begin(); it != term.coefficients.end();) {
            if (it->momenta.size() != 1U) {
                ++it;
                continue;
            }
            WickOperator transform = WickOperator::from_coefficient(*it, WICK_OP_TO_STRING_TEMP_BASE);
            if (transform.type == OperatorType::Undefined) {
                ++it;
            } else {
                term.operators.push_back(std::move(transform));
                it = term.coefficients.erase(it);
            }
        }
    }
    clean_wick_ordered_terms(result);
    return result;
}

void include_cross_contraction_in_templates(std::vector<WickOperatorTemplate>& templates) {
    const std::size_t old_size = templates.size();
    duplicate_n_inplace(templates, 1);
    for (std::size_t i = 0U; i < old_size; ++i) {
        templates[i].indexComparison.push_back(IndexComparison{false, Index::WickA, Index::WickB});
        templates[i + old_size].indexComparison.push_back(IndexComparison{false, Index::WickB, Index::WickA});
    }
}

WickOrderedCollector multiply(const WickOrderedTerm& lhs,
                              const WickOrderedTerm& rhs,
                              std::vector<WickOperatorTemplate> templates) {
    include_cross_contraction_in_templates(templates);
    return multiply_impl(lhs, rhs, templates);
}

WickOrderedCollector commutator_impl(const WickOrderedTerm& lhs,
                                     const WickOrderedTerm& rhs,
                                     const std::vector<WickOperatorTemplate>& templates) {
    WickOrderedCollector result = multiply_impl(lhs, rhs, templates) - multiply_impl(rhs, lhs, templates);
    clean_wick_ordered_terms(result);
    return result;
}

WickOrderedCollector commutator(const WickOrderedTerm& lhs,
                                const WickOrderedTerm& rhs,
                                std::vector<WickOperatorTemplate> templates) {
    include_cross_contraction_in_templates(templates);
    return commutator_impl(lhs, rhs, templates);
}

WickOrderedCollector commutator(const WickOrderedCollector& lhs,
                                const WickOrderedCollector& rhs,
                                std::vector<WickOperatorTemplate> templates) {
    include_cross_contraction_in_templates(templates);
    WickOrderedCollector result;
    for (const auto& lhs_term : lhs) {
        for (const auto& rhs_term : rhs) {
            result += commutator_impl(lhs_term, rhs_term, templates);
        }
    }
    return result;
}

void clean_wick_ordered_terms(
    WickOrderedCollector& terms,
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

    for (WickOrderedCollector::iterator it = terms.begin(); it != terms.end();) {
        if (!(it->resolve_deltas())) {
            it = terms.erase(it);
            continue;
        }

        it->discard_zero_momenta();
        it->rename_sums();

        for (const auto& symmetry : symmetries) {
            symmetry->apply_to(it->operators);
        }
        for (auto& coeff : it->coefficients) {
            coeff.use_custom_symmetry();
            coeff.use_inversion_symmetry();
        }

        if (!it->creation_operators_to_the_left()) {
            it = terms.erase(it);
            continue;
        }
        if (sort_lexicographically(it->wick_expression.operators)) {
            it->multiplicity *= -1;
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
        if (sort_operators_by_indices(term.wick_expression.operators)) {
            term.multiplicity *= -1;
        }
        detail::structure_momentum_dependencies_impl(term, term.wick_expression.operators);
        for (auto& coeff : term.coefficients) {
            coeff.use_custom_symmetry();
            coeff.use_inversion_symmetry();
        }
        std::sort(term.coefficients.begin(), term.coefficients.end());
    }

    terms.combine_duplicates();

    terms.sort();
}

std::ostream& operator<<(std::ostream& os, const WickOrderedCollector& terms) {
    for (WickOrderedCollector::const_iterator it = terms.begin(); it != terms.end(); ++it) {
        os << "\t&" << *it;
        if (it != terms.end() - 1) {
            os << " \\\\";
        }
        os << "\n";
    }
    return os;
}

WickOrderedCollector& WickOrderedCollector::hermitian_conjugate_inplace() noexcept {
    for (auto& term : terms) {
        term.hermitian_conjugate_inplace();
    }
    return *this;
}

WickOrderedCollector WickOrderedCollector::hermitian_conjugate() const noexcept {
    WickOrderedCollector copy(*this);
    copy.hermitian_conjugate_inplace();
    return copy;
}

void WickOrderedCollector::sort() {
    for (auto& term : terms) {
        term.sort();
    }

    std::sort(terms.begin(), terms.end());
}

WickOrderedCollector& WickOrderedCollector::operator+=(const WickOrderedCollector& other) {
    append_vector(this->terms, other.terms);
    return *this;
}

WickOrderedCollector& WickOrderedCollector::operator-=(const WickOrderedCollector& other) {
    const std::size_t old_size = size();
    append_vector(this->terms, other.terms);

    for (std::size_t i = old_size; i < size(); ++i) {
        terms[i].multiplicity *= -1;
    }
    return *this;
}

WickOrderedCollector operator+(WickOrderedCollector lhs, const WickOrderedCollector& rhs) {
    return (lhs += rhs);
}

WickOrderedCollector operator-(WickOrderedCollector lhs, const WickOrderedCollector& rhs) {
    return (lhs -= rhs);
}

}  // namespace mrock::symbolic_operators::experimental