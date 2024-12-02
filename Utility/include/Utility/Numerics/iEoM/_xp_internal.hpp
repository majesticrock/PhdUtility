#pragma once
#include "../../VectorWrapper.hpp"
#include <string>
#include <numeric>
#include <algorithm>

namespace Utility::Numerics::iEoM { 
    template<class RealType>
	struct StartingState {
		using Vector = Eigen::Vector<RealType, Eigen::Dynamic>;
		Vector phase_state;
		Vector amplitude_state;
		std::string name;

		inline bool contains_amplitude_state() const noexcept {
			return amplitude_state.size() != 0;
		}
		inline bool contains_phase_state() const noexcept {
			return phase_state.size() != 0;
		}
		inline size_t size() const noexcept {
			return static_cast<size_t>(contains_amplitude_state()) + static_cast<size_t>(contains_phase_state());
		}
		inline Vector& operator[](int i) {
			if(i == 0) return phase_state;
			return amplitude_state;
		}
		inline const Vector& operator[](int i) const {
			if(i == 0) return phase_state;
			return amplitude_state;
		}
	};

    template<class RealType>
    auto phase_size(std::vector<StartingState<RealType>> const& states) {
        return std::count_if(states.begin(), states.end(), [](const StartingState<RealType>& state) { 
            return state.contains_phase_state();
            });
    }

    template<class RealType>
    auto amplitude_size(std::vector<StartingState<RealType>> const& states) {
        return std::count_if(states.begin(), states.end(), [](const StartingState<RealType>& state) { 
            return state.contains_amplitude_state();
            });
    }

    template<class RealType>
    int total_size(std::vector<StartingState<RealType>> const& states) {
        return std::accumulate(states.begin(), states.end(), int{}, [](int current, const StartingState<RealType>& state) { 
            return std::move(current) + state.size(); 
            });
    }

    template<class RealType>
    class PhaseIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = StartingState<RealType>;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type*;
        using reference = value_type&;

        PhaseIterator(std::vector<StartingState<RealType>>& vec, size_t pos = 0U)
            : m_vec(vec), m_pos(pos) { skip_invalid(); }
        static PhaseIterator begin(std::vector<StartingState<RealType>>& vec) {
            return PhaseIterator<RealType>(vec, 0U);
        }
        static PhaseIterator end(std::vector<StartingState<RealType>>& vec) {
            return PhaseIterator<RealType>(vec, vec.size());
        }

        PhaseIterator& operator++() {
            ++m_pos;
            skip_invalid();
            return *this;
        }
        PhaseIterator operator++(int) {
            PhaseIterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const PhaseIterator& other) const { return m_pos == other.m_pos; }
        bool operator!=(const PhaseIterator& other) const { return !(*this == other); }
        reference operator*() { return m_vec[m_pos]; }
        pointer operator->() { return &m_vec[m_pos]; }

    private:
        std::vector<StartingState<RealType>>& m_vec;
        size_t m_pos;
        void skip_invalid() {
            while (m_pos < m_vec.size() && !m_vec[m_pos].contains_phase_state()) {
                ++m_pos;
            }
        }
    };

    template<class RealType>
    class AmplitudeIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = StartingState<RealType>;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type*;
        using reference = value_type&;

        AmplitudeIterator(std::vector<StartingState<RealType>>& vec, size_t pos = 0U)
            : m_vec(vec), m_pos(pos) { skip_invalid(); }
        static AmplitudeIterator begin(std::vector<StartingState<RealType>>& vec) {
            return AmplitudeIterator<RealType>(vec, 0U);
        }
        static AmplitudeIterator end(std::vector<StartingState<RealType>>& vec) {
            return AmplitudeIterator<RealType>(vec, vec.size());
        }

        AmplitudeIterator& operator++() {
            ++m_pos;
            skip_invalid();
            return *this;
        }
        AmplitudeIterator operator++(int) {
            AmplitudeIterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const AmplitudeIterator& other) const { return m_pos == other.m_pos; }
        bool operator!=(const AmplitudeIterator& other) const { return !(*this == other); }
        reference operator*() { return m_vec[m_pos]; }
        pointer operator->() { return &m_vec[m_pos]; }

    private:
        std::vector<StartingState<RealType>>& m_vec;
        size_t m_pos;
        void skip_invalid() {
            while (m_pos < m_vec.size() && !m_vec[m_pos].contains_amplitude_state()) {
                ++m_pos;
            }
        }
    };
}