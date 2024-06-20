#pragma once
#include <iostream>
#include <map>

namespace SymbolicOperators {
    enum OperatorType { Number_Type = 0, CDW_Type, SC_Type, Eta_Type, Undefined_Type };
	inline const std::map<std::string, OperatorType> string_to_wick = {
		{"n", Number_Type}, {"g", CDW_Type}, {"f", SC_Type}, {"\\eta", Eta_Type}
	};
	std::ostream& operator<<(std::ostream& os, const OperatorType op);
}