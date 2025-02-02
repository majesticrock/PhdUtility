/**
 * @file WickOperatorTemplate.hpp
 * @brief Defines templates for creating  Wick operators from a set of normal operators.
 */

#pragma once

#include "Operator.hpp"
#include "WickOperator.hpp"
#include "KroneckerDelta.hpp"
#include <optional>
#include <algorithm>

namespace mrock::symbolic_operators {

	/**
	 * @struct IndexComparison
	 * @brief A structure for comparing indices. E.g. <n_k> merely requires that the spin indizes of the composing operators are identical,
	 * 		but <f_k> requires the first index to be spin down.
	 */
	struct IndexComparison {
		bool any_identical; ///< Indicates if any identical indices are identical.
		Index base{ Index::UndefinedIndex }; ///< The base index.
		Index other{ Index::UndefinedIndex }; ///< The other index.
	};

	/**
	 * @struct TemplateResult
	 * @brief A structure for storing the result of a template operation.
	 */
	struct TemplateResult {

		/**
		 * @struct SingleResult
		 * @brief A structure for storing a single result.
		 */
		struct SingleResult {
			int factor{}; ///< The factor of the result.
			WickOperator op; ///< The Wick operator.
			std::vector<KroneckerDelta<Index>> index_deltas; ///< The index deltas.

			/**
			 * @brief Clears KroneckerDelta objects that are one, i.e., delta_{a,a} = 1
			 */
			inline void clear_delta_equals_one();

			/**
			 * @brief Checks if the result contains an impossible delta, e.g., delta_{down,up}
			 * 
			 * @return true if the result contains an impossible delta and false otherwise.
			 */
			inline bool contains_impossible_delta() const;
		};

		std::vector<SingleResult> results; ///< The vector of single results.
		KroneckerDelta<Momentum> momentum_delta; ///< The momentum delta.

		/**
		 * @brief Default constructor for TemplateResult.
		 */
		TemplateResult() = default;

		/**
		 * @brief Constructs a TemplateResult object.
		 * 
		 * @param initial_size The initial size of the results vector.
		 * @param operator_type The type of the operator.
		 * @param base_momentum The base momentum.
		 */
		TemplateResult(size_t initial_size, OperatorType operator_type, const Momentum& base_momentum);

		/**
		 * @brief Creates a null TemplateResult.
		 * 
		 * @return TemplateResult A null TemplateResult.
		 */
		inline static TemplateResult null_result() { return {}; }

		/**
		 * @brief Applies an operation on a range of results.
		 * 
		 * @tparam UnaryOperation The type of the operation.
		 * @param operation The operation to apply.
		 * @param begin The beginning of the range.
		 * @param n The number of elements in the range.
		 */
		template<class UnaryOperation>
		void operation_on_range(const UnaryOperation& operation, size_t begin, size_t n) {
			for (size_t i = begin; i < begin + n; ++i)
			{
				operation(results[i]);
			}
		}

		/**
		 * @brief Applies an operation on each result.
		 * 
		 * @tparam UnaryOperation The type of the operation.
		 * @param operation The operation to apply.
		 */
		template<class UnaryOperation>
		void operation_on_each(const UnaryOperation& operation) {
			for (auto& res : results)
			{
				operation(res);
			}
		}

		/**
		 * @brief Adds an index delta to a range of results.
		 * 
		 * @param index The index delta to add.
		 * @param begin The beginning of the range.
		 * @param n The number of elements in the range.
		 */
		inline void add_index_delta_range(const KroneckerDelta<Index>& index, size_t begin, size_t n);

		/**
		 * @brief Adds an index delta to each result.
		 * 
		 * @param index The index delta to add.
		 */
		inline void add_index_delta(const KroneckerDelta<Index>& index);

		/**
		 * @brief Creates a branch in the results vector.
		 * 
		 * @return size_t The size of the current results vector.
		 */
		size_t create_branch();

		/**
		 * @brief Clears impossible results.
		 */
		void clear_impossible();

		/**
		 * @brief Cleans up the results by clearing deltas that are one and removing impossible results.
		 */
		void clean_up();

		/**
		 * @brief Checks if the TemplateResult is valid.
		 * 
		 * @return true if the TemplateResult is valid and false otherwise.
		 */
		inline explicit operator bool() const { return !this->results.empty(); }
	};

	/**
	 * @class WickOperatorTemplate
	 * @brief A template for creating Wick operators.
	 */
	struct WickOperatorTemplate {
		std::vector<IndexComparison> indexComparison; ///< The vector of index comparisons.
		Momentum momentum_difference; ///< The momentum difference.
		OperatorType type; ///< The type of the operator.
		bool is_sc_type{}; ///< Indicates if the operator is of SC type.

		/**
		 * @brief Creates a WickOperator from two operators if possible.
		 * 
		 * @param left The left operator.
		 * @param right The right operator.
		 * @return TemplateResult The result of the creation.
		 */
		TemplateResult create_from_operators(const Operator& left, const Operator& right) const;

	private:
		/**
		 * @brief Handles the creation of SC type operators.
		 * 
		 * @param left The left operator.
		 * @param right The right operator.
		 * @return TemplateResult The result of the creation.
		 */
		TemplateResult _handle_sc_type(const Operator& left, const Operator& right) const;

		/**
		 * @brief Handles the creation of NUM type operators.
		 * 
		 * @param left The left operator.
		 * @param right The right operator.
		 * @return TemplateResult The result of the creation.
		 */
		TemplateResult _handle_num_type(const Operator& left, const Operator& right) const;
	};

	// Inline definitions
	void TemplateResult::SingleResult::clear_delta_equals_one() {
		auto new_end = std::remove_if(this->index_deltas.begin(), this->index_deltas.end(), [](const KroneckerDelta<Index>& delta) {
			return delta.first == delta.second;
			});
		this->index_deltas.erase(new_end, this->index_deltas.end());
	}
	bool TemplateResult::SingleResult::contains_impossible_delta() const {
		return std::any_of(this->index_deltas.begin(), this->index_deltas.end(), [](const KroneckerDelta<Index>& delta) {
			return (!is_mutable(delta.first) && !is_mutable(delta.second) && delta.first != delta.second);
			});
	}

	void TemplateResult::add_index_delta_range(const KroneckerDelta<Index>& index, size_t begin, size_t n) {
		operation_on_range([&index](SingleResult& res) { res.index_deltas.push_back(index); }, begin, n);
	}
	void TemplateResult::add_index_delta(const KroneckerDelta<Index>& index) {
		operation_on_each([&index](SingleResult& res) { res.index_deltas.push_back(index); });
	}
} // namespace mrock::symbolic_operators