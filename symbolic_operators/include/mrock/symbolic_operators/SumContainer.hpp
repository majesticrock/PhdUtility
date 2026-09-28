#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_SUMCONTAINER_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_SUMCONTAINER_HPP
/**
 * @file SumContainer.hpp
 * @brief Defines the SumContainer structure and related operators for symbolic operations.
 */

#include "IndexWrapper.hpp"
#include "MomentumSymbol.hpp"
#include "SymbolicSum.hpp"

#include <compare>
#include <ostream>
#include <vector>

namespace mrock::symbolic_operators {
/**
 * @typedef IndexSum
 * @brief Typedef for SymbolicSum with Index type.
 */
typedef SymbolicSum<Index> IndexSum;

/**
 * @typedef MomentumSum
 * @brief Typedef for SymbolicSum with MomentumSymbol::name_type type.
 */
typedef SymbolicSum<MomentumSymbol::name_type> MomentumSum;

/**
 * @struct SumContainer
 * @brief A container for holding symbolic sums of momenta and spins.
 *
 * Sums are contained within the \c SumContainer class.
 * It hosts both sums of momenta and sums of spins, each one is accessible via the appropriate class member
 * and its \c operator[], e.g., \c container.momenta[i].
 *
 * @sa Index, MomentumSymbol, MomentumSymbol::name_type
 */
struct SumContainer {
    MomentumSum momenta;  ///< Container for momentum sums.
    IndexSum spins;       ///< Container for spin sums.

    /**
     * @brief Serializes the SumContainer object.
     * @tparam Archive The type of the archive.
     * @param ar The archive to serialize to.
     * @param version The version of the serialization.
     */
    template <class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar& this->momenta;
        ar& this->spins;
    }

    /**
     * @brief Appends another SumContainer to this one.
     * @param other The other SumContainer to append.
     * @return Reference to this SumContainer.
     */
    SumContainer& append(const SumContainer& other);

    /**
     * @brief Appends a MomentumSum to this SumContainer.
     * @param other The MomentumSum to append.
     * @return Reference to this SumContainer.
     */
    SumContainer& append(const MomentumSum& other);

    /**
     * @brief Appends an IndexSum to this SumContainer.
     * @param other The IndexSum to append.
     * @return Reference to this SumContainer.
     */
    SumContainer& append(const IndexSum& other);

    /**
     * @brief Prepends another SumContainer to this one.
     * @param other The other SumContainer to prepend.
     * @return Reference to this SumContainer.
     */
    SumContainer& prepend(const SumContainer& other);

    /**
     * @brief Prepends a MomentumSum to this SumContainer.
     * @param other The MomentumSum to prepend.
     * @return Reference to this SumContainer.
     */
    SumContainer& prepend(const MomentumSum& other);

    /**
     * @brief Prepends an IndexSum to this SumContainer.
     * @param other The IndexSum to prepend.
     * @return Reference to this SumContainer.
     */
    SumContainer& prepend(const IndexSum& other);

    /**
     * @brief Sorts the summations vectors of both \c momenta and \c spins using the default comparison operator
     */
    void sort();

    /**
     * @brief Pushes back a momentum into the momenta container.
     * @param momentum The momentum to push back.
     */
    inline void push_back(const MomentumSymbol::name_type momentum);

    /**
     * @brief Pushes back a spin into the spins container.
     * @param spin The spin to push back.
     */
    inline void push_back(const Index spin);

    /**
     * @brief Checks if the container has any momenta.
     * @return True if the container has momenta, false otherwise.
     */
    inline bool has_momentum() const noexcept;

    /**
     * @brief Checks if the container has any spins.
     * @return True if the container has spins, false otherwise.
     */
    inline bool has_spins() const noexcept;

    inline std::strong_ordering operator<=>(const SumContainer& other) const {
        if (auto cmp = momenta.size() <=> other.momenta.size(); cmp != 0)
            return cmp;
        if (auto cmp = spins.size() <=> other.spins.size(); cmp != 0)
            return cmp;
        
        if (auto cmp = momenta <=> other.momenta; cmp != 0)
            return cmp;
        return spins <=> other.spins;
    };

    bool operator==(const SumContainer& other) const = default;
};

/**
 * @brief Stream insertion operator for SumContainer.
 * @param os The output stream.
 * @param sums The SumContainer to insert into the stream.
 * @return Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& os, const SumContainer& sums) {
    os << sums.momenta << sums.spins;
    return os;
}

// Inline definitions
void SumContainer::push_back(const MomentumSymbol::name_type momentum) {
    this->momenta.push_back(momentum);
}
void SumContainer::push_back(const Index spin) {
    this->spins.push_back(spin);
}
bool SumContainer::has_momentum() const noexcept {
    return !momenta.empty();
}
bool SumContainer::has_spins() const noexcept {
    return !spins.empty();
}
}  // namespace mrock::symbolic_operators
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_SUMCONTAINER_HPP
