#pragma once
#include <iostream>
#include <map>

namespace mrock::symbolic_operators {
	/**
	 * @enum Index
	 * @brief Enumeration representing various symbolic indices.
	 * 
	 * The indices include:
	 * - Number_Type: <n> (0).
	 * - CDW_Type: <c_{k+Q}^+ c_k> (1).
	 * - SC_Type: <c_{-k down} c_{k up}> (2).
	 * - Eta_Type: <c_{-k-Q down} c_{k up}> (3).
	 * - Undefined_Type: Represents an undefined type (4).
	 */
	enum OperatorType { Number_Type = 0, CDW_Type, SC_Type, Eta_Type, Undefined_Type = 255 };
	
	/**
	 * @var string_to_wick
	 * @brief A map that associates string representations with their corresponding Wick operators values.
	 */
	inline const std::map<std::string, OperatorType> string_to_wick = {
		{"n", Number_Type}, {"g", CDW_Type}, {"f", SC_Type}, {"\\eta", Eta_Type}
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