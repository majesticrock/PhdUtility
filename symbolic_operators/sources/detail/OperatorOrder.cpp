#include <mrock/symbolic_operators/detail/OperatorOrder.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/Term.hpp>

#include <cstddef>
#include <vector>

namespace mrock::symbolic_operators {

// Returns true if a fermionic sign change was induced by the swap, false otherwise
bool _perform_operator_swap(std::vector<Operator>& operators, std::size_t i, std::size_t j) {
    std::size_t min = std::min(i, j);
    const std::size_t max = std::max(i, j);
    bool parity = operators[min].is_fermion && operators[max].is_fermion;
    while (++min < max) {
        if (operators[min].is_fermion && operators[max].is_fermion) {
            parity = !parity;
        }
    }
    std::swap(operators[i], operators[j]);
    return parity;
}

bool sort_operators_by_indices(std::vector<Operator>& operators)
{
    auto general_swap_predicate = [&operators](const std::size_t& i, const std::size_t& index_pos) {
        if (operators[i - 1].is_daggered) {
            if (operators[i - 1].indices.front() == Index::SpinUp)
                return false;
            if (operators[i].indices.front() == Index::SpinDown)
                return false;

            if (operators[i - 1].indices.front() == Index::SpinDown || operators[i].indices.front() == Index::SpinUp) {
                return true;
            }

            if (operators[i - 1].indices[index_pos] > operators[i].indices[index_pos]) {
                return true;
            }
        } else {
            // The comparison operator and SpinUp <-> SpinDown are the only changes
            if (operators[i - 1].indices.front() == Index::SpinDown)
                return false;
            if (operators[i].indices.front() == Index::SpinUp)
                return false;

            if (operators[i - 1].indices.front() == Index::SpinUp || operators[i].indices.front() == Index::SpinDown) {
                return true;
            }

            if (operators[i - 1].indices[index_pos] < operators[i].indices[index_pos]) {
                return true;
            }
        }
        return false;
    };

    bool parity = false;

    std::size_t n = operators.size();
    std::size_t new_n{};
    // First sort so that the spins indices are always ordered the same way
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
            while (j < operators[i - 1].indices.size() && j < operators[i].indices.size()) {
                if (general_swap_predicate(i, j)) {
                    if (_perform_operator_swap(operators, i - 1, i)) {
                        parity = !parity;
                    }
                    new_n = i;
                    break;
                }
                ++j;
            }
        }
        n = new_n;
    }
    return parity;
}

void structure_momentum_dependencies(Term& term)
{
    detail::structure_momentum_dependencies_impl(term, term.operators);
}

} // namespace mrock::symbolic_operators