/**
 * @file Wick.hpp
 * @brief Functions for applying Wick's theorem and manipulating Wick terms.
 */

#ifndef _WICK
#define _WICK

#include "WickTerm.hpp"
#include "WickSymmetry.hpp"
#include <vector>
#include <memory>

namespace mrock::symbolic_operators {

	/**
	 * @brief Identifies Wick operators in a given Wick term.
	 * 
	 * @param source The source Wick term.
	 * @param operator_templates The vector of Wick operator templates.
	 * @return WickTermCollector The collected Wick terms.
	 */
	WickTermCollector identify_wick_operators(const WickTerm& source, const std::vector<WickOperatorTemplate>& operator_templates);

	/**
	 * @brief Applies Wick's theorem to a set of terms.
	 * 
	 * @param terms The vector of terms.
	 * @param operator_templates The vector of Wick operator templates.
	 * @param reciever The WickTermCollector to receive the results.
	 */
	void wicks_theorem(const std::vector<Term>& terms, const std::vector<WickOperatorTemplate>& operator_templates, WickTermCollector& reciever);

	/**
	 * @brief Clears eta terms from the WickTermCollector. Intended for use if <eta>=0
	 * 
	 * @param terms The WickTermCollector containing the terms.
	 */
	void clear_etas(WickTermCollector& terms);
    
	/**
	 * @brief Cleans Wick terms using the provided symmetries.
	 * 
	 * @param terms The WickTermCollector containing the terms.
	 * @param symmetries The vector of unique pointers to WickSymmetry objects.
	 */
	void clean_wicks(WickTermCollector& terms, const std::vector<std::unique_ptr<WickSymmetry>>& symmetries = std::vector<std::unique_ptr<WickSymmetry>>{});
} // namespace mrock::symbolic_operators

#endif