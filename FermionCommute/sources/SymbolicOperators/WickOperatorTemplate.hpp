#pragma once

#include "Operator.hpp"
#include "WickOperator.hpp"
#include "KroneckerDelta.hpp"
#include <optional>
#include <algorithm>

namespace SymbolicOperators {
	struct IndexComparison {
		bool any_identical;
		Index base{ UndefinedIndex };
		Index other{ UndefinedIndex };
	};

	struct TemplateResult {
		struct SingleResult {
			int factor{};
			WickOperator op;
			std::vector<KroneckerDelta<Index>> index_deltas;

			inline void clear_delta_equals_one() {
				auto new_end = std::remove_if(this->index_deltas.begin(), this->index_deltas.end(), [](const KroneckerDelta<Index>& delta) {
					return delta.first == delta.second;
					});
				this->index_deltas.erase(new_end, this->index_deltas.end());
			};
			inline bool contains_impossible_delta() const {
				return std::any_of(this->index_deltas.begin(), this->index_deltas.end(), [](const KroneckerDelta<Index>& delta) {
					return (!is_mutable(delta.first) && !is_mutable(delta.second) && delta.first != delta.second);
					});
			};
		};
		std::vector<SingleResult> results;
		KroneckerDelta<Momentum> momentum_delta;

		TemplateResult() = default;
		TemplateResult(size_t initial_size, OperatorType operator_type, const Momentum& base_momentum);

		inline static TemplateResult null_result() {
			return {};
		};

		template<class UnaryOperation>
		void operation_on_range(const UnaryOperation& operation, size_t begin, size_t n) {
			for (size_t i = begin; i < begin + n; ++i)
			{
				operation(results[i]);
			}
		}
		template<class UnaryOperation>
		void operation_on_each(const UnaryOperation& operation) {
			for (auto& res : results)
			{
				operation(res);
			}
		}
		inline void add_index_delta_range(const KroneckerDelta<Index>& index, size_t begin, size_t n) {
			operation_on_range([&index](SingleResult& res) { res.index_deltas.push_back(index); }, begin, n);
		}
		inline void add_index_delta(const KroneckerDelta<Index>& index) {
			operation_on_each([&index](SingleResult& res) { res.index_deltas.push_back(index); });
		}
		inline size_t create_branch() {
			const size_t current_size{ results.size() };
			results.reserve(2 * current_size);
			std::copy_n(results.begin(), current_size, std::back_inserter(results));
			return current_size;
		}

		inline void clear_impossible() {
			auto new_end = std::remove_if(this->results.begin(), this->results.end(), [](const SingleResult& result) {
				return result.contains_impossible_delta();
				});
			this->results.erase(new_end, this->results.end());
		};
		inline void clean_up() {
			for (auto& result : results) {
				result.clear_delta_equals_one();
			}
			this->clear_impossible();
		}

		inline explicit operator bool() const {
			return !this->results.empty();
		}
	};

	struct WickOperatorTemplate {
		std::vector<IndexComparison> indexComparison;
		Momentum momentum_difference;
		OperatorType type;
		bool is_sc_type{};

		// Returns the corresponding WickOperator if construction is possible
		// Otherwise, it returns an empty optional
		TemplateResult createFromOperators(const Operator& left, const Operator& right) const;

	private:
		TemplateResult _handle_sc_type(const Operator& left, const Operator& right) const;
		TemplateResult _handle_num_type(const Operator& left, const Operator& right) const;
	};
}