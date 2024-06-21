#ifndef _WICK
#define _WICK

#include "WickTerm.hpp"
#include "WickSymmetry.hpp"
#include <vector>
#include <memory>

namespace SymbolicOperators
{
	WickTermCollector identifyWickOperators(const WickTerm& source, const std::vector<WickOperatorTemplate>& operator_templates);

	void wicks_theorem(const std::vector<Term>& terms, const std::vector<WickOperatorTemplate>& operator_templates, WickTermCollector& reciever);

    // Call this function if <eta> = 0
	void clearEtas(WickTermCollector& terms);
    
	void cleanWicks(WickTermCollector& terms, const std::vector<std::unique_ptr<WickSymmetry>>& symmetries = std::vector<std::unique_ptr<WickSymmetry>>{});
}
#endif