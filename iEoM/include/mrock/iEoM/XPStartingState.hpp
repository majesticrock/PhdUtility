#pragma once
#include <Eigen/Dense>
#include <string>
#include <algorithm>

namespace mrock::iEoM {
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
}