#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_ABSTRACTTERM_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_ABSTRACTTERM_HPP
/**
 * @file AbstractTerm.hpp
 * @brief Defines the AbstractTerm structure, which serves as a parent to both \c Term and \c WickTerm.
 */

#include "detail/container_helper.hpp"
#include "Coefficient.hpp"
#include "Fractional.hpp"
#include "KroneckerDelta.hpp"
#include "KroneckerDeltaUtility.hpp"
#include "SumContainer.hpp"

#include <vector>
#include <functional>

namespace mrock::symbolic_operators {
/**
 * @class AbstractTerm
 * @brief Serves as a parent to \c Term and \c WickTerm.
 * Defines and implements certain methods that are used by both classes as to avoid code duplication
 *
 * @sa Term, WickTerm
 */
template <class tOperatorType>
class AbstractTerm {
protected:
    // Here, we define buffers to be used later on.
    // Later, we want to bring all terms to the same notation.
    // This will entail renaming sums, so that the first sum is always q, the second p, the third r, etc.
    // To avoid name clashes, we first rename all sums to a buffer name, e.g., the first sum is renamed to :, the second
    // to ;, the third to |, etc. As of now, this is limited to 11 sums (which I doubt will become a problem in the
    // future). If more is needed, the buffer_list and name_list need to be extended, and the N_BUFFER constant needs to
    // be changed accordingly.
    constexpr static int N_BUFFER = 11;
    constexpr static MomentumSymbol::name_type name_list[N_BUFFER] = {'q', 'p', 'r', 's', 't', 'u',
                                                                      'v', 'w', 'x', 'y', 'z'};
    constexpr static MomentumSymbol::name_type buffer_list[N_BUFFER] = {':', ';', '|', '?', '!', '.',
                                                                        '-', '_', '+', '/', '='};

    virtual void for_each_momentum_except_deltas(const std::function<void(Momentum&)>& f) 
    {
        for (auto& op : operators) {
            f(op.momentum);
        }
        for (auto& coeff : coefficients) {
            for (auto& momentum : coeff.momenta) {
                f(momentum);
            }
        }
    }

    virtual void for_each_index_except_deltas(const std::function<void(IndexWrapper&)>& f) 
    {
        for (auto& op : operators) {
            f(op.indizes);
        }
        for (auto& coeff : coefficients) {
            f(coeff.indizes);
        }
    }

public:
    IntFractional multiplicity;             ///< Multiplicity of the term.
    std::vector<Coefficient> coefficients;  ///< Coefficients of the term.
    SumContainer sums;                      ///< Sum container for the term. Contains e.g. \sum_{k,l} \sum_{sigma}
    std::vector<KroneckerDelta<Momentum>> delta_momenta;  ///< Kronecker delta for momenta.
    std::vector<KroneckerDelta<Index>> delta_indizes;     ///< Kronecker delta for indices.
    std::vector<tOperatorType>
        operators;  ///< Operators in the term, if empty the term is considered to contain the identiy operator

    /**
     * @brief Default constructor.
     */
    AbstractTerm() = default;

    /**
     * @brief Constructs a Term with a summation over momenta and spins and multiple coefficients and Kronecker deltas
     *
     * @param _multiplicity The multiplicity of the term
     * @param _sums The sums
     * @param _coefficients The coefficients
     * @param _operators The operators
     * @param _delta_momenta The Kronecker deltas for the momenta
     * @param _delta_indizes The Kronecker deltas for the indizes
     */
    AbstractTerm(const IntFractional& _multiplicity,
                 const std::vector<Coefficient>& _coefficients,
                 const SumContainer& _sums,
                 const std::vector<KroneckerDelta<Momentum>>& _delta_momenta,
                 const std::vector<KroneckerDelta<Index>>& _delta_indizes,
                 const std::vector<tOperatorType>& _operators)
        : multiplicity{_multiplicity},
          coefficients{_coefficients},
          sums{_sums},
          delta_momenta{_delta_momenta},
          delta_indizes{_delta_indizes},
          operators{_operators} {};

    /**
     * @brief Constructs a Term with a summation over momenta and spins and multiple coefficients
     *
     * @param _multiplicity The multiplicity of the term
     * @param _sums The sums
     * @param _coefficients The coefficients
     * @param _operators The operators
     */
    AbstractTerm(const IntFractional& _multiplicity,
                 const std::vector<Coefficient>& _coefficients,
                 const SumContainer& _sums,
                 const std::vector<tOperatorType>& _operators = std::vector<tOperatorType>())
        : multiplicity{_multiplicity}, coefficients{_coefficients}, sums{_sums}, operators{_operators} {};

    /**
     * @brief Constructs a Term with a summation over momenta and spins and one coefficient
     *
     * @param _multiplicity The multiplicity of the term
     * @param _sums The sums
     * @param _coefficient The coefficient
     * @param _operators The operators
     */
    AbstractTerm(const IntFractional& _multiplicity,
                 const Coefficient& _coefficient,
                 const SumContainer& _sums,
                 const std::vector<tOperatorType>& _operators = std::vector<tOperatorType>())
        : multiplicity{_multiplicity}, coefficients(1, _coefficient), sums{_sums}, operators{_operators} {};

    /**
     * @brief Constructs a Term with no summations
     *
     * @param _multiplicity The multiplicity of the term
     * @param _coefficient The coefficient
     * @param _operators The operators
     */
    AbstractTerm(const IntFractional& _multiplicity,
                 const Coefficient& _coefficient,
                 const std::vector<tOperatorType>& _operators = std::vector<tOperatorType>())
        : multiplicity{_multiplicity}, coefficients(1, _coefficient), operators{_operators} {};

    /**
     * @brief Constructs a Term with no summations and no coefficient
     *
     * @param _multiplicity The multiplicity of the term
     * @param _coefficient The coefficient
     * @param _operators The operators
     */
    explicit AbstractTerm(const IntFractional& _multiplicity,
                          const std::vector<tOperatorType>& _operators = std::vector<tOperatorType>())
        : multiplicity{_multiplicity}, operators{_operators} {};

    /**
     * @brief Virtual destructor.
     */
    virtual ~AbstractTerm() = default;

    /**
     * @brief Replaces a momentum symbol everywhere in the term.
     *
     * @param replaceWhat The momentum symbol to replace.
     * @param replaceWith The momentum expression that replaces the symbol.
     * @param skip A predicate that skips specific momentum Kronecker deltas during replacement.
     */
    void replace_each_momentum(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith, 
        std::function<bool(std::vector<KroneckerDelta<Momentum>>::iterator)> skip = [](auto) { return false; });

    /**
     * @brief Replaces an index everywhere in the term.
     *
     * @param target The index to replace.
     * @param replace_with The index that replaces the target.
     * @param skip A predicate that skips specific index Kronecker deltas during replacement.
     */
    void replace_each_index(Index target, Index replace_with, 
        std::function<bool(std::vector<KroneckerDelta<Index>>::iterator)> skip = [](auto) { return false; });

    /**
     * @brief Discard momenta that are zero (k+0=k)
     */
    void discard_zero_momenta();

    /**
     * @brief Resolves the Kronecker deltas of the momenta in the term.
     * @return True if successful, false otherwise.
     */
    bool resolve_momentum_deltas();

    /**
     * @brief Resolves the Kronecker deltas of the indizes in the term.
     * @return True if successful, false otherwise.
     */
    bool resolve_index_deltas();

    /**
     * @brief Renames the sums in the term.
     */
    void rename_sums();

    /**
     * @brief Checks if the term is an identity term.
     *
     * @return true if the term is an identity term.
     * @return false otherwise.
     */
    bool is_identity() const noexcept;

    /**
     * @brief Flips the sign of the term.
     */
    void flip_sign();

    /**
     * @brief Gets the operators in the term.
     * @return The operators.
     */
    const std::vector<tOperatorType>& get_operators() const;

    /**
     * @brief Inverts a momentum in the term.
     *
     * @param what The momentum to invert.
     */
    void invert_momentum(const MomentumSymbol::name_type what);

    /**
     * @brief Inverts a momentum sum in the term.
     *
     * @param what The momentum sum to invert.
     */
    void invert_momentum_sum(const MomentumSymbol::name_type what);

    /**
     * @brief Removes a momentum contribution from the term.
     *
     * @param value The momentum value to remove.
     */
    void remove_momentum_contribution(const MomentumSymbol::name_type value);

    /**
     * @brief Renames momenta in the term.
     * @param what The momentum to rename.
     * @param to The new momentum.
     */
    void rename_momenta(const MomentumSymbol::name_type what, const MomentumSymbol::name_type to);

    /**
     * @brief Transforms a momentum sum in the term.
     * @param what The momentum to transform.
     * @param to The new momentum.
     * @param new_sum_index The new sum index.
     */
    void transform_momentum_sum(const MomentumSymbol::name_type what,
                                const Momentum to,
                                const MomentumSymbol::name_type new_sum_index);

    /**
     * @brief Restructures the momentum summations so that a given \c current Momentum is turned into \c should_be
     * 
     * @param current The current state of the Momentum to target
     * @param should_be The desired final state; constructed via \c Momentum(should_be)
     * @param do_not_use Optionally provide a vector of symbols that should not be used for the transformation
     * The idea is, if you already ordered 'q', you probably do not want to destroy what you already achieved.
     * So you can pass 'q' as \c do_not_use and the algorithm will skip it.
     */
    void redistribute_momenta(const Momentum& current, const MomentumSymbol::name_type& should_be,
        const std::vector<MomentumSymbol::name_type>& do_not_use = {});

    /**
     * @brief Renames all occuring \c what to \c to
     * 
     * @param what The Index to look for.
     * @param to The Index that \c what should be changed to.
     */
    void rename_indizes(const Index what, const Index to);

    /**
     * @brief Restructures the index summations so that a given \c current Index is turned into \c should_be
     * 
     * @param current The current state of the Index to target
     * @param should_be The desired final state
     * @param do_not_use Optionally provide a vector of symbols that should not be used for the transformation
     * The idea is, if you already ordered 'sigma', you probably do not want to destroy what you already achieved.
     * So you can pass 'sigma' as \c do_not_use and the algorithm will skip it.
     */
    void redistribute_indizes(const Index current, const Index should_be,
        const std::vector<Index>& do_not_use = {});
};

// Implementations
template <class tOperatorType>
void AbstractTerm<tOperatorType>::replace_each_momentum(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith, 
    std::function<bool(std::vector<KroneckerDelta<Momentum>>::iterator)> skip)
{
    for_each_momentum_except_deltas([&replaceWhat, &replaceWith](Momentum& momentum){
        momentum.replace_occurances(replaceWhat, replaceWith);
    });
    for (std::vector<KroneckerDelta<Momentum>>::iterator it = delta_momenta.begin(); it != delta_momenta.end(); ++it) {
        if (skip(it)) {
            continue;
        }
        it->first.replace_occurances(replaceWhat, replaceWith);
        it->second.replace_occurances(replaceWhat, replaceWith);
    }
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::replace_each_index(Index target, Index replace_with, 
    std::function<bool(std::vector<KroneckerDelta<Index>>::iterator)> skip) 
{
    for_each_index_except_deltas([&target, &replace_with](IndexWrapper& indizes) {
        indizes.replace_index(target, replace_with);
    });
    for (std::vector<KroneckerDelta<Index>>::iterator it = delta_indizes.begin(); it != delta_indizes.end(); ++it) {
        if (skip(it)) {
            continue;
        }
        if (it->first == target) {
            it->first = replace_with;
        }
        if (it->second == target) {
            it->second = replace_with;
        }
    }
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::discard_zero_momenta() {
    for (auto& op : operators) {
        op.momentum.remove_zeros();
    }
    for (auto& coeff : coefficients) {
        coeff.momenta.remove_zeros();
    }
}

template <class tOperatorType>
bool AbstractTerm<tOperatorType>::resolve_momentum_deltas() {
    if (is_always_zero(delta_momenta))
        return false;

    // Remove delta^2
    remove_delta_squared(this->delta_momenta);
    // Erase delta_k,k etc
    remove_delta_is_one(this->delta_momenta);

    for (auto delta_it = delta_momenta.begin(); delta_it != delta_momenta.end();) {
        delta_it->first -= delta_it->second;
        delta_it->second = Momentum();

        if (delta_it->first.momentum_list.empty() && delta_it->second.momentum_list.empty()) {
            // 0 = Q can never be achieved
            if (delta_it->first.add_Q != delta_it->second.add_Q)
                return false;
            // delta_(0,0) = 1
            delta_it = delta_momenta.erase(delta_it);
            continue;
        }

        MomentumSymbol resolve_to{*(delta_it->first.begin())};
        bool found_sum{};

        // Try to find a momentum in the sums that is also in the delta
        for (auto sum_it = sums.momenta.begin(); sum_it != sums.momenta.end(); ++sum_it) {
            const auto found_it =
                std::find_if(delta_it->first.begin(), delta_it->first.end(),
                             [&sum_it](const MomentumSymbol& symbol) { return symbol.name == *sum_it; });
            if (found_it != delta_it->first.end()) {
                resolve_to = *found_it;
                sums.momenta.erase(sum_it);
                found_sum = true;
                break;
            }
        }

        delta_it->second = Momentum(resolve_to);
        delta_it->second.flip_momentum();
        delta_it->first += delta_it->second;
        if (delta_it->second.front().factor < 0) {
            delta_it->first.flip_momentum();
            delta_it->second.flip_momentum();
        }

        // Fractional momenta were never implemented, since they were never needed.
        // If they are ever needed, you unfortunately have to implement them yourself
        for (MomentumSymbol& symbol : delta_it->first) {
            assert(symbol.factor % delta_it->second.front().factor == 0);
            symbol.factor /= delta_it->second.front().factor;
        }

        // Replace set the delta everywhere, e.g., delta_{k,l+q} would replace each k with l+q
        replace_each_momentum(delta_it->second.front().name, delta_it->first, [&delta_it](auto it) { return it == delta_it; });

        if (found_sum) {
            delta_it = delta_momenta.erase(delta_it);
        } else {
            ++delta_it;
        }
    }

    // Make sure that delta.first has always exactly one momentum symbol (or delta_{0,0})
    for (auto& delta : delta_momenta) {
        if (delta.first.empty() && !delta.second.empty()) {
            std::swap(delta.first, delta.second);
        }
        if (delta.first.size() > 1U) {
            const Momentum shift = delta.first - Momentum(delta.first.front());
            delta.first -= shift;
            delta.second -= shift;
        }
        if (delta.first.front().factor < 0) {
            delta.first.flip_momentum();
            delta.second.flip_momentum();
        }
    }

    // Remove delta^2
    remove_delta_squared(this->delta_momenta);
    // Erase delta_k,k etc
    remove_delta_is_one(this->delta_momenta);

    return true;
}

template <class tOperatorType>
bool AbstractTerm<tOperatorType>::resolve_index_deltas() {
    if (is_always_zero(delta_indizes))
        return false;

    // Remove delta^2
    remove_delta_squared(this->delta_indizes);
    // Erase delta_k,k etc
    remove_delta_is_one(this->delta_indizes);

    for (auto delta_it = delta_indizes.begin(); delta_it != delta_indizes.end();) {
        Index to_resolve{Index::UndefinedIndex};
        Index change_to{Index::UndefinedIndex};
        bool found_sum{};

        // try to find a spin in the sums that is also in the delta
        auto sum_it = std::find_if(sums.spins.begin(), sums.spins.end(),
                                   [&delta_it](const Index& idx) { return idx == delta_it->first; });
        if (sum_it != sums.spins.end()) {
            to_resolve = delta_it->first;
            change_to = delta_it->second;
            found_sum = true;
        } else {
            sum_it = std::find_if(sums.spins.begin(), sums.spins.end(),
                                  [&delta_it](const Index& idx) { return idx == delta_it->second; });
            if (sum_it != sums.spins.end()) {
                to_resolve = delta_it->second;
                change_to = delta_it->first;
                found_sum = true;
            }
        }

        if (to_resolve == Index::UndefinedIndex) {
            if (is_mutable(delta_it->first)) {
                to_resolve = delta_it->first;
                change_to = delta_it->second;
            } else if (is_mutable(delta_it->second)) {
                to_resolve = delta_it->second;
                change_to = delta_it->first;
            } else if (delta_it->first != delta_it->second) {
                // Two differing, immutable indizes can never be equal
                return false;
            }
        }

        if(remove_delta_is_one(delta_indizes)) {
            delta_it = delta_indizes.begin();
            continue;
        }
        if(remove_delta_squared(delta_indizes)) {
            delta_it = delta_indizes.begin();
            continue;
        }
        replace_each_index(to_resolve, change_to, [&delta_it](auto it) { return it == delta_it; });

        if (found_sum) {
            sums.spins.erase(sum_it);
            delta_it = delta_indizes.erase(delta_it);
        } else {
            ++delta_it;
        }
    }

    // Remove delta^2
    remove_delta_squared(this->delta_indizes);
    // Erase delta_k,k etc
    remove_delta_is_one(this->delta_indizes);

    return true;
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::rename_sums() {
    for (std::size_t i = 0U; i < sums.momenta.size(); ++i) {
        if (i >= N_BUFFER) {
            std::cerr << "More than " << N_BUFFER << "momenta, time to implement this..." << std::endl;
            break;
        }
        if (sums.momenta[i] == name_list[i])
            continue;

        replace_each_momentum(sums.momenta[i], Momentum(buffer_list[i]), [](auto) { return true; });
        sums.momenta[i] = name_list[i];
    }

    for (std::size_t i = 0U; i < sums.momenta.size(); ++i) {
        replace_each_momentum(buffer_list[i], Momentum(name_list[i]), [](auto) { return true; });
    }

    std::sort(sums.spins.begin(), sums.spins.end());
    for (std::size_t i=0U; i < sums.spins.size(); ++i){
        const Index should_be = static_cast<Index>(static_cast<std::size_t>(Index::Sigma) + i);
        if (static_cast<std::size_t>(sums.spins[i]) != static_cast<std::size_t>(should_be)) {
            const Index old_sum = sums.spins[i];
            sums.spins[i] = should_be;
            replace_each_index(old_sum, should_be, [](auto) { return true; });
        }
    }
}

template <class tOperatorType>
bool AbstractTerm<tOperatorType>::is_identity() const noexcept {
    return this->operators.empty();
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::flip_sign() {
    this->multiplicity *= -1;
}

template <class tOperatorType>
const std::vector<tOperatorType>& AbstractTerm<tOperatorType>::get_operators() const {
    return this->operators;
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::invert_momentum(const MomentumSymbol::name_type what) {
    for_each_momentum_except_deltas([what](Momentum& momentum){
        momentum.flip_single(what);
    });
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::invert_momentum_sum(const MomentumSymbol::name_type what) {
    if (std::find(sums.momenta.begin(), sums.momenta.end(), what) == sums.momenta.end()) {
        throw std::invalid_argument(
            "You are trying to perform a sum transformation on a momentum that is not being summed over!");
    }
    invert_momentum(what);
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::remove_momentum_contribution(const MomentumSymbol::name_type value) {
    for (auto& coeff : coefficients) {
        coeff.remove_momentum_contribution(value);
    }
    for (auto& op : operators) {
        op.remove_momentum_contribution(value);
    }
    for (auto& delta : delta_momenta) {
        delta.first.remove_contribution(value);
        delta.second.remove_contribution(value);
    }
    std::erase_if(sums.momenta.summations, [&](const MomentumSymbol::name_type sum_idx) { return sum_idx == value; });
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::transform_momentum_sum(const MomentumSymbol::name_type what,
                                                        const Momentum to,
                                                        const MomentumSymbol::name_type new_sum_index) {
    auto pos = std::find(sums.momenta.begin(), sums.momenta.end(), what);
    if (pos == sums.momenta.end()) {
        throw std::invalid_argument(
            "You are trying to perform a sum transformation on a momentum that is not being summed over!");
    } else {
        *pos = new_sum_index;
    }

    replace_each_momentum(what, to, [](auto) { return true; });
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::rename_momenta(const MomentumSymbol::name_type what, const MomentumSymbol::name_type to) {
    if (what == to)
        return;
    for (auto& mom_sum : sums.momenta) {
        if (mom_sum == to) {
            throw std::invalid_argument("You are replacing a momentum sum with an index that already exists!");
        }
        if (mom_sum == what) {
            mom_sum = to;
        }
    }

    replace_each_momentum(what, Momentum(to));
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::redistribute_momenta(const Momentum& current, const MomentumSymbol::name_type& should_be,
    const std::vector<MomentumSymbol::name_type>& do_not_use)
{
    if (current == Momentum(should_be) || current.empty()) return;

    std::size_t i=0;
    MomentumSymbol::name_type transformer = current[i].name;

    while (!sums.momenta.is_summed_over(transformer) || exists_in(do_not_use, transformer)) {
        if (++i >= current.size()) {
            throw std::invalid_argument("There is no summation that would allow the desired momentum transformation!");
        }
        transformer = current[i].name;
        
    }

    if (current[i].factor < 0) {
        invert_momentum_sum(transformer);
    }
    Momentum target = -current;
    target += Momentum(transformer);
    target += Momentum('?');

    transform_momentum_sum(transformer, target, '?');

    rename_momenta(should_be, transformer);
    rename_momenta('?', should_be);
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::rename_indizes(const Index what, const Index to)
{
    if (what == to) return;

    for (auto& index_sum : sums.spins) {
        if (index_sum == to) {
            throw std::invalid_argument("You are replacing an index sum with an index that already exists!");
        }
        if (index_sum == what) {
            index_sum = to;
        }
    }

    replace_each_index(what, to);
}

template <class tOperatorType>
void AbstractTerm<tOperatorType>::redistribute_indizes(const Index current, const Index should_be,
        const std::vector<Index>& do_not_use)
{
    if (current == should_be) return;
    if (exists_in(do_not_use, current)) return;

    if (sums.spins.is_summed_over(current)) {
        rename_indizes(should_be, Index::PlaceHolderIndex);
        rename_indizes(current, should_be);
        rename_indizes(Index::PlaceHolderIndex, current);
    }
    else {
        throw std::invalid_argument("There is no summation that would allow the desired index transformation!");
    }
}

}  // namespace mrock::symbolic_operators
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_ABSTRACTTERM_HPP
