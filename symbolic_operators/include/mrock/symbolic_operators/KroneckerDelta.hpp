/**
 * @file KroneckerDelta.hpp
 * @brief Defines the KroneckerDelta structure used in symbolic operators.
 */

#pragma once
#include <utility>
#include <iostream>
#include <mrock/utility/defines_arithmetic_operators.hpp>

namespace mrock::symbolic_operators {

	/**
	 * @class KroneckerDelta
	 * @brief A structure representing the Kronecker Delta.
	 * 
	 * @tparam T The type of the elements.
	 */
	template<typename T>
	struct KroneckerDelta {
		T first{};
		T second{};

		/**
		 * @brief Serializes the KroneckerDelta object.
		 * 
		 * @tparam Archive The type of the archive.
		 * @param ar The archive object.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& first;
			ar& second;
		}

		/**
		 * @brief Checks if the Kronecker Delta is one.
		 * 
		 * @return true if first equals second.
		 * @return false otherwise.
		 */
		constexpr bool isOne() const {
			return first == second;
		};
	};

	/**
	 * @brief Creates a KroneckerDelta object.
	 * 
	 * @tparam T The type of the elements.
	 * @param first The first element.
	 * @param second The second element.
	 * @return constexpr auto A KroneckerDelta object.
	 */
	template<typename T>
	constexpr auto make_delta(const T& first, const T& second) {
		return KroneckerDelta<T>{first, second};
	}

	/**
	 * @brief Creates a KroneckerDelta object with rvalue references.
	 * 
	 * @tparam T The type of the elements.
	 * @param first The first element.
	 * @param second The second element.
	 * @return constexpr auto A KroneckerDelta object.
	 */
	template<typename T>
	constexpr auto make_delta(std::decay_t<T>&& first, std::decay_t<T>&& second) {
		return KroneckerDelta<std::decay_t<T>>{std::move(first), std::move(second)};
	}

	/**
	 * @brief Equality operator for KroneckerDelta.
	 * 
	 * @tparam T The type of the elements.
	 * @param lhs The left-hand side KroneckerDelta.
	 * @param rhs The right-hand side KroneckerDelta.
	 * @return true if the two KroneckerDelta objects are equal.
	 * @return false otherwise.
	 */
	template<typename T>
	bool operator==(const KroneckerDelta<T>& lhs, const KroneckerDelta<T>& rhs) {
		if (lhs.first == rhs.first && lhs.second == rhs.second) return true;
		if (lhs.first == rhs.second && lhs.second == rhs.first) return true;
		return false;
	}

	/**
	 * @brief Inequality operator for KroneckerDelta.
	 * 
	 * @tparam T The type of the elements.
	 * @param lhs The left-hand side KroneckerDelta.
	 * @param rhs The right-hand side KroneckerDelta.
	 * @return true if the two KroneckerDelta objects are not equal.
	 * @return false otherwise.
	 */
	template<typename T>
	bool operator!=(const KroneckerDelta<T>& lhs, const KroneckerDelta<T>& rhs) {
		return !(lhs == rhs);
	}

	/**
	 * @brief Addition assignment operator for KroneckerDelta.
	 * 
	 * @tparam T The type of the elements.
	 * @param lhs The left-hand side KroneckerDelta.
	 * @param rhs The right-hand side element.
	 * @return KroneckerDelta<T>& The updated KroneckerDelta.
	 */
	template<typename T> requires mrock::utility::defines_plus<T>::value
	inline KroneckerDelta<T>& operator+=(KroneckerDelta<T>& lhs, T& rhs) {
		lhs.first += rhs;
		lhs.second += rhs;
		return lhs;
	}

	/**
	 * @brief Subtraction assignment operator for KroneckerDelta.
	 * 
	 * @tparam T The type of the elements.
	 * @param lhs The left-hand side KroneckerDelta.
	 * @param rhs The right-hand side element.
	 * @return KroneckerDelta<T>& The updated KroneckerDelta.
	 */
	template<typename T> requires mrock::utility::defines_minus<T>::value
	inline KroneckerDelta<T>& operator-=(KroneckerDelta<T>& lhs, const T& rhs) {
		lhs.first -= rhs;
		lhs.second -= rhs;
		return lhs;
	}

	/**
	 * @brief Addition operator for KroneckerDelta.
	 * 
	 * @tparam T The type of the elements.
	 * @param lhs The left-hand side KroneckerDelta.
	 * @param rhs The right-hand side element.
	 * @return KroneckerDelta<T> The resulting KroneckerDelta.
	 */
	template<typename T> requires mrock::utility::defines_plus<T>::value
	inline KroneckerDelta<T> operator+(KroneckerDelta<T> lhs, T const& rhs) {
		return (lhs += rhs);
	}

	/**
	 * @brief Subtraction operator for KroneckerDelta.
	 * 
	 * @tparam T The type of the elements.
	 * @param lhs The left-hand side KroneckerDelta.
	 * @param rhs The right-hand side element.
	 * @return KroneckerDelta<T> The resulting KroneckerDelta.
	 */
	template<typename T> requires mrock::utility::defines_minus<T>::value
	inline KroneckerDelta<T> operator-(KroneckerDelta<T> lhs, T const& rhs) {
		return (lhs -= rhs);
	}

	/**
	 * @brief Stream insertion operator for KroneckerDelta.
	 * 
	 * @tparam T The type of the elements.
	 * @param os The output stream.
	 * @param delta The KroneckerDelta object.
	 * @return std::ostream& The updated output stream.
	 */
	template<typename T>
	inline std::ostream& operator<<(std::ostream& os, const KroneckerDelta<T>& delta) {
		os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		return os;
	};
} // namespace mrock::symbolic_operators