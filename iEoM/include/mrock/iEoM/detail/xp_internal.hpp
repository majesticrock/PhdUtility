#pragma once
#include <type_traits>
#include "../XPStartingState.hpp"

namespace mrock::iEoM::detail { 
    struct state_iterator_tag {};

    template<typename It>
    concept StateIterator = requires { 
            typename It::state_iterator_tag; 
        } && std::same_as<typename It::state_iterator_tag, mrock::iEoM::detail::state_iterator_tag>;

    template<typename It>
    concept ConstStateIterator = StateIterator<It> && It::is_const_iterator;

    template<class RealType, bool iteratesPhase, bool isConst>
    class StateIteratorImpl {
        using container_type = std::conditional_t<isConst,
            const std::vector<StartingState<RealType>>,
            std::vector<StartingState<RealType>>>;
        
        using value_type_impl = std::conditional_t<isConst,
            const StartingState<RealType>,
            StartingState<RealType>>;
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = StartingState<RealType>;
        using difference_type   = std::ptrdiff_t;
        using pointer           = value_type_impl*;
        using reference         = value_type_impl&;

        using Vector = typename StartingState<RealType>::Vector;
        static constexpr bool is_const_iterator = isConst;
        using state_iterator_tag = mrock::iEoM::detail::state_iterator_tag;

        StateIteratorImpl(container_type& vec, size_t pos = 0U)
            : m_vec(vec), m_pos(pos) { skip_invalid(); }
        static StateIteratorImpl begin(container_type& vec) {
            return StateIteratorImpl(vec, 0U);
        }
        static StateIteratorImpl end(container_type& vec) {
            return StateIteratorImpl(vec, vec.size());
        }

        StateIteratorImpl& operator++() {
            ++m_pos;
            skip_invalid();
            return *this;
        }
        StateIteratorImpl operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const StateIteratorImpl& other) const { return m_pos == other.m_pos; }
        bool operator!=(const StateIteratorImpl& other) const { return !(*this == other); }
        reference operator*() { return m_vec[m_pos]; }
        pointer operator->() { return &m_vec[m_pos]; }

        auto& state() const 
        { 
            if constexpr (iteratesPhase)
                return m_vec[m_pos].phase_state; 
            else 
                return m_vec[m_pos].amplitude_state;
        }
    private:
        container_type& m_vec;
        size_t m_pos;
        void skip_invalid() {
            if constexpr (iteratesPhase) {
                while (m_pos < m_vec.size() && !m_vec[m_pos].contains_phase_state()) {
                    ++m_pos;
                }
            }
            else {
                while (m_pos < m_vec.size() && !m_vec[m_pos].contains_amplitude_state()) {
                    ++m_pos;
                }
            }
        }
    };

    template<class RealType>
    using PhaseIterator = StateIteratorImpl<RealType, true, false>;
    template<class RealType>
    using ConstPhaseIterator = StateIteratorImpl<RealType, true, true>;

    template<class RealType>
    using AmplitudeIterator = StateIteratorImpl<RealType, false, false>;
    template<class RealType>
    using ConstAmplitudeIterator = StateIteratorImpl<RealType, false, true>;
}