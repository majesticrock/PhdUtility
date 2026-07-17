#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_KRONECKERDELTAUTILITY_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_KRONECKERDELTAUTILITY_HPP
/**
 * @file KroneckerDeltaUtility.hpp
 * @brief Utility functions for manipulating KroneckerDelta objects.
 */

#include "IndexWrapper.hpp"
#include "KroneckerDelta.hpp"
#include "Momentum.hpp"
#include "detail/defines_arithmetic_operators.hpp"

#include <concepts>
#include <cstddef>

namespace mrock::symbolic_operators {

/**
 * @brief Removes squared KroneckerDelta objects from the vector. Note that delta_{a,b}^N = delta_{a,b}.
 * Specialized for types that define + and - operators, i.e., that are linearly combinable.
 *
 * @tparam T The LinearlyCombinable (defines + and -) type of the elements.
 * @param deltas The vector of KroneckerDelta objects.
 */
template <LinearlyCombinable T>
    requires std::default_initializable<T>
void remove_delta_squared(std::vector<KroneckerDelta<T>>& deltas) {
    for (std::size_t i = 0U; i < deltas.size(); i++) {
        for (std::size_t j = i + 1; j < deltas.size(); j++) {
            const T buffer = (deltas[i].first - deltas[i].second) - (deltas[j].first - deltas[j].second);
            if (buffer == T()) {
                deltas.erase(deltas.begin() + j);
                break;
            }
        }
    }
}

template <class T>
void remove_delta_squared(std::vector<KroneckerDelta<T>>& deltas) {
    for (std::size_t i = 0U; i < deltas.size(); i++) {
        for (std::size_t j = i + 1; j < deltas.size(); j++) {
            if (deltas[i] == deltas[j]) {
                deltas.erase(deltas.begin() + j);
                --i;
                break;
            }
        }
    }
}

/**
 * @brief Removes KroneckerDelta objects that are one from the vector. Note that delta_{a,a} = 1.
 *
 * @tparam T The type of the elements.
 * @param deltas The vector of KroneckerDelta objects.
 */
template <class T>
void remove_delta_is_one(std::vector<KroneckerDelta<T>>& deltas) {
    auto new_end =
        std::remove_if(deltas.begin(), deltas.end(), [](const KroneckerDelta<T>& delta) { return delta.isOne(); });
    deltas.erase(new_end, deltas.end());
}

/**
 * @brief Checks if the vector of KroneckerDelta<Index> objects is always zero.
 *
 * @param deltas The vector of KroneckerDelta<Index> objects.
 * @return true if the vector is always zero.
 * @return false otherwise.
 */
inline bool is_always_zero(const std::vector<KroneckerDelta<Index>>& deltas) {
    return std::any_of(deltas.begin(), deltas.end(), [](const KroneckerDelta<Index>& delta) {
        return (delta.first != delta.second && (!is_mutable(delta.first) && !is_mutable(delta.second)));
    });
}

/**
 * @brief Checks if the vector of KroneckerDelta<Momentum> objects is always zero.
 *
 * @param deltas The vector of KroneckerDelta<Momentum> objects.
 * @return true if the vector is always zero.
 * @return false otherwise.
 */
inline bool is_always_zero(const std::vector<KroneckerDelta<Momentum>>& deltas) {
    return std::any_of(deltas.begin(), deltas.end(), [](const KroneckerDelta<Momentum>& delta) {
        return delta.first.differs_only_in_Q(delta.second);
    });
}

/**
 * @brief Removes double occurrences in a KroneckerDelta<Momentum> object.
 *
 * @param delta The KroneckerDelta<Momentum> object.
 */
inline void remove_double_occurances(KroneckerDelta<Momentum>& delta) {
    if (delta.first.add_Q) {
        delta.first.add_Q = false;
        delta.second.add_Q = !(delta.second.add_Q);
    }
    for (auto it = delta.first.momentum_list.begin(); it != delta.first.momentum_list.end();) {
        const int idx = delta.second.is_used_at(it->name);
        if (idx < 0) {
            ++it;
            continue;
        }

        delta.second.momentum_list[idx].factor -= it->factor;
        it = delta.first.momentum_list.erase(it);
        if (delta.second.momentum_list[idx].factor == 0) {
            delta.second.momentum_list.erase(delta.second.momentum_list.begin() + idx);
        }
    }
}
}  // namespace mrock::symbolic_operators
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_KRONECKERDELTAUTILITY_HPP
