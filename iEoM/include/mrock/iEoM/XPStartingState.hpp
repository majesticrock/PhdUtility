#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_XPSTARTINGSTATE_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_XPSTARTINGSTATE_HPP
#include <Eigen/Dense>
#include <string>
#include <numeric>
#include <utility>

namespace mrock::iEoM {
	/**
	 * @brief Represents an initial state for phase and amplitude variables.
	 *
	 * The starting state may contain a phase state vector, an amplitude state
	 * vector, or both. The name field can be used to identify the state.
	 *
	 * @tparam RealType Numeric type for the state vectors.
	 */
	template<class RealType>
	struct XPStartingState {
	    using Vector = Eigen::Vector<RealType, Eigen::Dynamic>;

	    /** Phase state vector. Empty when no phase state is present. */
	    Vector phase_state;

	    /** Amplitude state vector. Empty when no amplitude state is present. */
	    Vector amplitude_state;

	    /** Optional name for the starting state. */
	    std::string name;

	    /**
	     * @brief Returns true when an amplitude state vector is present.
	     * @return true if amplitude_state is non-empty.
	     */
	    inline bool contains_amplitude_state() const noexcept {
	        return amplitude_state.size() != 0;
	    }

	    /**
	     * @brief Returns true when a phase state vector is present.
	     * @return true if phase_state is non-empty.
	     */
	    inline bool contains_phase_state() const noexcept {
	        return phase_state.size() != 0;
	    }

	    /**
	     * @brief Returns the number of state vectors that are present.
	     * @return 0, 1, or 2 depending on which state vectors are non-empty.
	     */
	    inline size_t size() const noexcept {
	        return static_cast<size_t>(contains_amplitude_state()) + static_cast<size_t>(contains_phase_state());
	    }

	    /**
	     * @brief Accesses one of the stored state vectors by index.
	     * @param i Index of the state vector. 0 returns phase_state; any other value returns amplitude_state.
	     * @return Reference to the selected state vector.
	     */
	    inline Vector& operator[](int i) {
	        if(i == 0) return phase_state;
	        return amplitude_state;
	    }

	    /**
	     * @brief Accesses one of the stored state vectors by index.
	     * @param i Index of the state vector. 0 returns phase_state; any other value returns amplitude_state.
	     * @return Const reference to the selected state vector.
	     */
	    inline const Vector& operator[](int i) const {
	        if(i == 0) return phase_state;
	        return amplitude_state;
	    }
	};

	/**
	 * @brief Count how many starting states contain a phase state vector.
	 * @tparam RealType Numeric type for the state vectors.
	 * @param states Vector of starting states to inspect.
	 * @return Number of states with a non-empty phase_state.
	 */
	template<class RealType>
	auto phase_size(std::vector<XPStartingState<RealType>> const& states) {
	    return std::count_if(states.begin(), states.end(), [](const XPStartingState<RealType>& state) {
	        return state.contains_phase_state();
	    });
	}

	/**
	 * @brief Count how many starting states contain an amplitude state vector.
	 * @tparam RealType Numeric type for the state vectors.
	 * @param states Vector of starting states to inspect.
	 * @return Number of states with a non-empty amplitude_state.
	 */
	template<class RealType>
	auto amplitude_size(std::vector<XPStartingState<RealType>> const& states) {
	    return std::count_if(states.begin(), states.end(), [](const XPStartingState<RealType>& state) {
	        return state.contains_amplitude_state();
	    });
	}

	/**
	 * @brief Compute the total number of non-empty state vectors across all states.
	 * @tparam RealType Numeric type for the state vectors.
	 * @param states Vector of starting states to inspect.
	 * @return Total count of phase and amplitude vectors contained in all states.
	 */
	template<class RealType>
	int total_size(std::vector<XPStartingState<RealType>> const& states) {
	    return std::accumulate(states.begin(), states.end(), int{}, [](int current, const XPStartingState<RealType>& state) {
	        return std::move(current) + state.size();
	    });
	}
} // namespace mrock::iEoM
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_XPSTARTINGSTATE_HPP
