#include "WickOperatorTemplate.hpp"
#include <algorithm>
#include <iterator>
#include "Momentum.hpp"
#include "KroneckerDelta.hpp"
#include "KroneckerDeltaUtility.hpp"

namespace SymbolicOperators {
	TemplateResult::TemplateResult(size_t initial_size, OperatorType operator_type, const Momentum& base_momentum)
		: results(initial_size)
	{
		for (auto& result : results)
		{
			result.op.type = operator_type;
			result.op.momentum = base_momentum;
			result.factor = 1;
		}
	}

	TemplateResult WickOperatorTemplate::_handle_sc_type(const Operator& left, const Operator& right) const {
		// c_{-k-q} c_{k} or c_{k}^+ c_{-k-q}^+
		const Operator& base{ left.isDaggered ? right : left };
		const Operator& other{ left.isDaggered ? left : right };
		// q
		const Momentum momentum_diff{ -(base.momentum + other.momentum) };

		TemplateResult result(1U, this->type, base.momentum);
		result.results.front().op.isDaggered = left.isDaggered;
		result.momentum_delta = make_delta(this->momentum_difference, momentum_diff);

		for (size_t i = 0U; i < indexComparison.size(); ++i)
		{
			if (indexComparison[i].any_identical) {
				result.add_index_delta(make_delta(base.indizes[i], other.indizes[i]));
				result.operation_on_each([&base, &i](TemplateResult::SingleResult& res) {
					res.op.indizes.push_back(base.indizes[i]);
					});
				// c c can be swapped for the cost of a sign
				const size_t previous_size{ result.create_branch() };
				result.operation_on_range([](TemplateResult::SingleResult& res) { res.factor *= -1; }, previous_size, previous_size);
				result.operation_on_range([&other](TemplateResult::SingleResult& res) { res.op.momentum = other.momentum; }, previous_size, previous_size);
			}
			else {
				const size_t previous_size{ result.create_branch() };
				result.add_index_delta_range(make_delta(base.indizes[i], indexComparison[i].base), 0U, previous_size);
				result.add_index_delta_range(make_delta(other.indizes[i], indexComparison[i].other), 0U, previous_size);

				// c c can be swapped for the cost of a sign
				result.add_index_delta_range(make_delta(base.indizes[i], indexComparison[i].other), previous_size, previous_size);
				result.add_index_delta_range(make_delta(other.indizes[i], indexComparison[i].base), previous_size, previous_size);
				result.operation_on_range([](TemplateResult::SingleResult& res) { res.factor *= -1; }, previous_size, previous_size);
				result.operation_on_range([&other](TemplateResult::SingleResult& res) { res.op.momentum = other.momentum; }, previous_size, previous_size);
			}
		}
		for (auto& res : result.results) {
			remove_delta_is_one(res.index_deltas);
			remove_delta_squared(res.index_deltas);
		}
		return result;
	}
	TemplateResult WickOperatorTemplate::_handle_num_type(const Operator& left, const Operator& right) const {
		// c_{k}^+ c_{k+q}
		// q
		const Momentum momentum_diff = right.momentum - left.momentum;
		KroneckerDelta<Momentum> momentum_delta{ this->momentum_difference, momentum_diff };

		std::vector<KroneckerDelta<Index>> index_delta;

		TemplateResult result(1U, this->type, left.momentum);
		result.results.front().op.isDaggered = false;
		result.momentum_delta = make_delta(this->momentum_difference, momentum_diff);

		for (size_t i = 0U; i < indexComparison.size(); ++i)
		{
			if (indexComparison[i].any_identical) {
				result.add_index_delta(make_delta(left.indizes[i], right.indizes[i]));
				result.operation_on_each([&left, &i](TemplateResult::SingleResult& res) {
					res.op.indizes.push_back(left.indizes[i]);
					});
			}
			else {
				const size_t previous_size{ result.create_branch() };
				result.add_index_delta_range(make_delta(left.indizes[i], indexComparison[i].base), 0U, previous_size);
				result.add_index_delta_range(make_delta(right.indizes[i], indexComparison[i].other), 0U, previous_size);
			}
		}
		for (auto& res : result.results) {
			remove_delta_is_one(res.index_deltas);
			remove_delta_squared(res.index_deltas);
		}
		return result;
	}

	TemplateResult WickOperatorTemplate::createFromOperators(const Operator& left, const Operator& right) const {
		if (this->is_sc_type) {
			if (left.isDaggered != right.isDaggered)
				return TemplateResult::null_result();
			return this->_handle_sc_type(left, right);
		}

		if (left.isDaggered == right.isDaggered)
			return TemplateResult::null_result();
		// The input needs to be normal ordered.
		// This means that the left input must be daggered here
		assert(left.isDaggered);
		return this->_handle_num_type(left, right);
	}
}