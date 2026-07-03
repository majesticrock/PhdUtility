/**
 * @file OperatorType.hpp
 * @brief Defines the OperatorType enum and related functions for symbolic operators.
 */

#pragma once
#include <iostream>
#include <map>

namespace mrock::symbolic_operators {
	/**
	 * @enum Index
	 * @brief Enumeration representing various symbolic indices.
	 * 
	 * The indices include:
	 * - Number: <n> (0).
	 * - CDW: <c_{k+Q}^+ c_k> (1).
	 * - SC: <c_{-k down} c_{k up}> (2).
	 * - Eta: <c_{-k-Q down} c_{k up}> (3).
	 * - Undefined: Represents an undefined type (4).
	 */
	enum class OperatorType { Number = 0, CDW, SC, Eta, Undefined = 255 };
	
	/**
	 * @brief Checks whether a given \c OperatorType is a pair creation or annihilation operators
	 * 
	 * @param type The \c OperatorType to check
	 * @return \c true if \c type is either \c SC or \c Eta and \c false otherwise
	 */
	constexpr bool is_bcs_like(OperatorType type) {
		return (type == OperatorType::SC || type == OperatorType::Eta);
	}

	/**
	 * @var string_to_wick
	 * @brief A map that associates string representations with their corresponding Wick operators values.
	 */
	inline const std::map<std::string, OperatorType> string_to_wick = {
		{"n", OperatorType::Number}, {"g", OperatorType::CDW}, {"f", OperatorType::SC}, {"\\eta", OperatorType::Eta}
	};

	/**
	 * @brief Overloads the stream insertion operator for the Index enum.
	 * 
	 * @param os The output stream.
	 * @param index The OperatorType to insert into the stream.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const OperatorType op);
} // namespace mrock::symbolic_operators