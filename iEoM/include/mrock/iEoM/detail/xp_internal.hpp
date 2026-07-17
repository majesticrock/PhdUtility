#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_XP_INTERNAL_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_XP_INTERNAL_HPP
#include "../XPStartingState.hpp"

#include <type_traits>

namespace mrock::iEoM::detail {
struct state_iterator_tag {};

template <typename It>
concept StateIterator = requires { typename It::state_iterator_tag; } &&
                        std::same_as<typename It::state_iterator_tag, mrock::iEoM::detail::state_iterator_tag>;

template <typename It>
concept ConstStateIterator = StateIterator<It> && It::is_const_iterator;

/**
 * @brief Iterator implementation for filtering phase or amplitude starting states.
 *
 * This iterator traverses a vector of StartingState objects and skips entries
 * that do not contain the requested state type (phase or amplitude).
 *
 * @tparam RealType Numeric type used by the contained StartingState vectors.
 * @tparam iteratesPhase When true, the iterator traverses phase states; otherwise amplitude states.
 * @tparam isConst When true, the iterator is a const iterator.
 */
template <class RealType, bool iteratesPhase, bool isConst>
class StateIteratorImpl {
    using container_type = std::
        conditional_t<isConst, const std::vector<XPStartingState<RealType>>, std::vector<XPStartingState<RealType>>>;

    using value_type_impl = std::conditional_t<isConst, const XPStartingState<RealType>, XPStartingState<RealType>>;

public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = XPStartingState<RealType>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type_impl*;
    using reference = value_type_impl&;

    using Vector = typename XPStartingState<RealType>::Vector;
    static constexpr bool is_const_iterator = isConst;
    using state_iterator_tag = mrock::iEoM::detail::state_iterator_tag;

    /**
     * @brief Construct an iterator for a vector of starting states.
     *
     * @param vec Reference to the starting state container.
     * @param pos Initial position within the container.
     */
    StateIteratorImpl(container_type& vec, std::size_t pos = 0U) : m_vec(vec), m_pos(pos) { skip_invalid(); }
    /**
     * @brief Create an iterator referring to the first valid state.
     *
     * @param vec Reference to the starting state container.
     * @return Iterator at the beginning of the valid state sequence.
     */
    static StateIteratorImpl begin(container_type& vec) { return StateIteratorImpl(vec, 0U); }
    /**
     * @brief Create an iterator representing the end of the container.
     *
     * @param vec Reference to the starting state container.
     * @return Iterator one past the last element.
     */
    static StateIteratorImpl end(container_type& vec) { return StateIteratorImpl(vec, vec.size()); }

    /**
     * @brief Advance the iterator to the next valid state.
     *
     * @return Reference to the incremented iterator.
     */
    StateIteratorImpl& operator++() {
        ++m_pos;
        skip_invalid();
        return *this;
    }
    /**
     * @brief Post-increment the iterator.
     *
     * @return Previous iterator position.
     */
    StateIteratorImpl operator++(int) {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    bool operator==(const StateIteratorImpl& other) const { return m_pos == other.m_pos; }
    bool operator!=(const StateIteratorImpl& other) const { return !(*this == other); }
    reference operator*() { return m_vec[m_pos]; }
    pointer operator->() { return &m_vec[m_pos]; }

    /**
     * @brief Access the current phase or amplitude state vector.
     *
     * @return Reference to the currently selected state vector.
     */
    auto& state() const {
        if constexpr (iteratesPhase)
            return m_vec[m_pos].phase_state;
        else
            return m_vec[m_pos].amplitude_state;
    }

    /**
     * @brief Returns the number of the current state, while skipping invalid states
     * Example: Let there be 3 types of starting states, but only the first and the last have
     * a phase state. Then, we have a phase iterator object that points to the last state,
     * and ask it \c get_state_number()
     * It will count: 0 ( \c begin() ), [1 skip because there is no phase state], 1 (pointing to the last element)
     * and return 1.
     * Note that m_pos will contain 2 and is thus different than the result of this function.
     *
     * @return The number of the current state, while skipping invalid states
     */
    std::size_t get_state_number() const {
        std::size_t n{};
        StateIteratorImpl init = StateIteratorImpl::begin(m_vec);
        while (init != *this) {
            ++init;
            ++n;
        }
        return n;
    }

private:
    container_type& m_vec;
    std::size_t m_pos;
    /**
     * @brief Advance the position until a valid state is found.
     *
     * Skips starting state entries that do not contain the selected state type.
     */
    void skip_invalid() {
        if constexpr (iteratesPhase) {
            while (m_pos < m_vec.size() && !m_vec[m_pos].contains_phase_state()) {
                ++m_pos;
            }
        } else {
            while (m_pos < m_vec.size() && !m_vec[m_pos].contains_amplitude_state()) {
                ++m_pos;
            }
        }
    }
};

template <class RealType>
using PhaseIterator = StateIteratorImpl<RealType, true, false>;
template <class RealType>
using ConstPhaseIterator = StateIteratorImpl<RealType, true, true>;

template <class RealType>
using AmplitudeIterator = StateIteratorImpl<RealType, false, false>;
template <class RealType>
using ConstAmplitudeIterator = StateIteratorImpl<RealType, false, true>;
}  // namespace mrock::iEoM::detail
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_XP_INTERNAL_HPP
