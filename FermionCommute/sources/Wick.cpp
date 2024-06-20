#include "Wick.hpp"
#include <Utility/RangeUtility.hpp>
#include "KroneckerDeltaUtility.hpp"
#include <Utility/Numerics/MathFunctions.hpp>
#include <variant>
#include <numeric>

namespace SymbolicOperators {
	void wick_processor(const std::vector<Operator>& remaining, WickTermCollector& reciever_list, std::variant<WickTerm, Term> buffer)
	{
		if (remaining.empty()) {
			reciever_list.push_back(std::get<WickTerm>(buffer));
			return;
		}
		for (size_t i = 1U; i < remaining.size(); ++i)
		{
			if (std::holds_alternative<Term>(buffer)) {
				WickTerm temp(std::get<Term>(buffer));
				buffer = temp;
			}
			if ((i % 2) == 0) {
				std::get<WickTerm>(buffer).multiplicity *= -1;
			}
			std::get<WickTerm>(buffer).temporary_operators.reserve(std::get<WickTerm>(buffer).temporary_operators.size() + 2);
			std::get<WickTerm>(buffer).temporary_operators.push_back(remaining[0]);
			std::get<WickTerm>(buffer).temporary_operators.push_back(remaining[i]);

			std::vector<Operator> copy_operators = remaining;
			copy_operators.erase(copy_operators.begin() + i);
			copy_operators.erase(copy_operators.begin());
			wick_processor(copy_operators, reciever_list, buffer);

			// delete last two elements, as they are to be updated in the next iteration
			std::get<WickTerm>(buffer).temporary_operators.pop_back();
			std::get<WickTerm>(buffer).temporary_operators.pop_back();
			if ((i % 2) == 0) {
				std::get<WickTerm>(buffer).multiplicity *= -1;
			}
		}
	}
	
	// Sorts the operators in 'term' according to Wick's theorem within 'temporary_operators'
	// Afterwards, these can be rewritten in terms of 'WickOperator's.
	WickTermCollector prepare_wick(const std::vector<Term>& terms)
	{
		WickTermCollector prepared_wick;
		const size_t estimated_size = std::accumulate(terms.begin(), terms.end(), size_t{}, [](size_t current, const Term& term) {
			return current + Utility::Numerics::double_factorial(term.getOperators().size());
			});

		prepared_wick.reserve(estimated_size);
		for (const auto& term : terms) {
			if (term.isIdentity()) {
				prepared_wick.push_back(WickTerm(term));
			}
			else {
				wick_processor(term.getOperators(), prepared_wick, term);
			}
		}

		return prepared_wick;
	}

	WickTermCollector identifyWickOperators(const WickTerm& source, const std::vector<WickOperatorTemplate>& operator_templates)
	{
		WickTermCollector ret;
		ret.push_back(source);
		ret.front().temporary_operators.clear();

		for (size_t i = 0U; i < source.temporary_operators.size(); i += 2U)
		{
			std::vector<TemplateResult> template_results;
			for (const auto& operator_template : operator_templates) {
				auto template_result = operator_template.createFromOperators(source.temporary_operators[i], source.temporary_operators[i + 1U]);
				if (template_result)
					template_results.push_back(std::move(template_result));
			}

			const size_t current_size = ret.size();
			const size_t number_additional_elements = std::accumulate(template_results.begin(), template_results.end(), size_t{}, [](size_t current, const TemplateResult& tr) {
				return current + tr.results.size();
				});
			if (number_additional_elements > 1U) {
				Utility::duplicate_n_inplace(ret, number_additional_elements - 1U);
			}

			size_t template_result_it{};
			size_t old_it{};
			for (const auto& tr : template_results) {
				old_it = template_result_it;
				for (const auto& tr_result : tr.results) {
					for (auto it = ret.begin(); it != ret.begin() + current_size; ++it) {
						(it + template_result_it * current_size)->includeTemplateResult(tr_result);
					}
					++template_result_it;
				}
				std::for_each(ret.begin() + old_it * current_size, ret.begin() + current_size * template_result_it, [&tr](WickTerm& ret_element) {
					if (!tr.momentum_delta.isOne())
						ret_element.delta_momenta.push_back(tr.momentum_delta);
					});
			}
		}

		return ret;
	}

	void wicks_theorem(const std::vector<Term>& terms, const std::vector<WickOperatorTemplate>& operator_templates, WickTermCollector& reciever)
	{
		WickTermCollector prepared_wick = prepare_wick(terms);

		for (auto& w_term : prepared_wick) {
			Utility::append_if(reciever, identifyWickOperators(w_term, operator_templates), [](const WickTerm& wick) {
				return !(is_always_zero(wick.delta_indizes) || is_always_zero(wick.delta_momenta));
				});
		}
	}

	void clearEtas(WickTermCollector& terms)
	{
		for (auto it = terms.begin(); it != terms.end();) {
			bool isEta = false;
			for (const auto& op : it->operators) {
				if (op.type == Eta_Type) {
					isEta = true;
					break;
				}
			}
			if (isEta) {
				it = terms.erase(it);
			}
			else {
				++it;
			}
		}
	}

	void cleanWicks(WickTermCollector& terms, const std::vector<std::unique_ptr<WickSymmetry>>& symmetries /*= std::vector<std::unique_ptr<WickSymmetry>>{}*/)
	{
		for (auto& term : terms) {
			for (std::vector<Coefficient>::iterator it = term.coefficients.begin(); it != term.coefficients.end();) {
				if (it->name == "") {
					it = term.coefficients.erase(it);
				}
				else {
					++it;
				}
			}
		}
		for (WickTermCollector::iterator it = terms.begin(); it != terms.end();) {
			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discardZeroMomenta();
			if (!(it->computeSums())) {
				it = terms.erase(it);
				continue;
			}
			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discardZeroMomenta();
			it->renameSums();
			it->sort();

			for(const auto& symmetry : symmetries) {
				symmetry->apply_to(*it);
			}

			//it->applyPhaseSymmetry();
			//it->applySpinSymmetry();
			//it->applyTranslationalSymmetry();

			for (auto jt = it->sums.spins.begin(); jt != it->sums.spins.end();)
			{
				if (it->usesIndex(*jt)) {
					++jt;
				}
				else {
					// We are assuming there are only spin indizes here (spin 1/2)
					// If another kind of index arises I have to readress this section.
					it->multiplicity *= 2;
					jt = it->sums.spins.erase(jt);
				}
			}
			// sort momentum lists in coefficients
			for (auto& coeff : it->coefficients) {
				coeff.momenta.sort();
			}
			++it;
		}

		// remove duplicates
		for (int i = 0; i < terms.size(); i++)
		{
			for (int j = i + 1; j < terms.size(); j++)
			{
				if (terms[i] == terms[j]) {
					terms[i].multiplicity += terms[j].multiplicity;
					terms.erase(terms.begin() + j);
					--i;
					break;
				}
			}
		}
		// removes any terms that have a 0 prefactor
		for (auto it = terms.begin(); it != terms.end();)
		{
			if (it->multiplicity == 0) {
				it = terms.erase(it);
			}
			else {
				++it;
			}
		}

		// Sort terms
		for (size_t i = 0; i < terms.size(); i++)
		{
			for (size_t j = i + 1; j < terms.size(); j++)
			{
				if (terms[i].delta_momenta.empty() && terms[j].delta_momenta.size() > 0) {
					std::swap(terms[i], terms[j]);
				}
				if (terms[i].delta_momenta.size() > 0 && terms[j].delta_momenta.size() > 0) {
					if (terms[i].delta_momenta.size() < terms[j].delta_momenta.size()) {
						std::swap(terms[i], terms[j]);
					}
					else if (terms[i].delta_momenta.size() == terms[j].delta_momenta.size()) {
						if (terms[i].delta_momenta[0].second.add_Q && !(terms[j].delta_momenta[0].second.add_Q)) {
							std::swap(terms[i], terms[j]);
						}
						else if (terms[i].coefficients.size() > 0) {
							if (terms[j].coefficients[0].name <terms[i].coefficients[0].name) {
								std::swap(terms[i], terms[j]);
							}
						}
					}
				}
				else if (terms[i].delta_momenta.empty() && terms[j].delta_momenta.empty()) {
					if (terms[i].coefficients.size() > 0) {
						if (terms[j].coefficients[0].name <terms[i].coefficients[0].name) {
							std::swap(terms[i], terms[j]);
						}
					}
				}
			}
		}
	}
}