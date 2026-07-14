#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_DETAIL_CONTAINER_HELPER_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_DETAIL_CONTAINER_HELPER_HPP
/**
 * @file container_helper.hpp
 * @brief Utility helpers for sequence-like container operations.
 *
 * This header provides generic helper functions for duplicating container
 * contents and appending elements from one container to another, including
 * overloads for move-aware append operations and predicate-based filtering.
 */

#include <iterator>
#include <algorithm>

namespace mrock::symbolic_operators {
/**
 * @brief Duplicates the existing contents of a container in place.
 *
 * The contents of `target` are copied `n` times to the end of the container.
 * The function reserves the final size before insertion to minimize
 * reallocations.
 *
 * @tparam Vector Container type supporting `size()`, `begin()`, `end()`, and back insertion.
 * @param target Container whose contents are duplicated.
 * @param n Number of additional copies to append.
 */
	template<class Vector>
	void duplicate_n_inplace(Vector& target, std::size_t n) {
		const std::size_t original_size = target.size();
		target.reserve(original_size * (n + 1U));
		for (std::size_t i = 0U; i < n; ++i) {
			std::copy_n(target.begin(), original_size, std::back_inserter(target));
		}
	}

/**
 * @brief Appends the contents of a source container to a target container.
 *
 * @tparam Vector Container type supporting `insert()` and iterators.
 * @param target Destination container.
 * @param source Source container whose elements are copied.
 */
	template<class Vector>
	void append_vector(Vector& target, const Vector& source) {
		target.insert(target.end(), source.begin(), source.end());
	}

/**
 * @brief Appends the contents of a source container to a target container using move semantics.
 *
 * @tparam Vector Container type supporting `insert()` and move iterators.
 * @param target Destination container.
 * @param source Source container whose elements are moved.
 */
	template<class Vector>
	void append_vector(Vector& target, Vector&& source) {
		target.insert(target.end(), std::make_move_iterator(source.begin()), std::make_move_iterator(source.end()));
	}

/**
 * @brief Appends elements from the source container that satisfy a predicate.
 *
 * @tparam Vector Container type.
 * @tparam UnaryPred Predicate callable returning `bool` for each element.
 * @param target Destination container.
 * @param source Source container whose elements are copied.
 * @param predicate Predicate used to filter appended elements.
 */
	template<class Vector, class UnaryPred>
	void append_if(Vector& target, const Vector& source, const UnaryPred& predicate) {
		std::copy_if(source.begin(), source.end(), std::back_inserter(target), predicate);
	}

/**
 * @brief Appends elements from the source container that satisfy a predicate using move semantics.
 *
 * @tparam Vector Container type.
 * @tparam UnaryPred Predicate callable returning `bool` for each element.
 * @param target Destination container.
 * @param source Source container whose elements are moved.
 * @param predicate Predicate used to filter appended elements.
 */
	template<class Vector, class UnaryPred>
	void append_if(Vector& target, Vector&& source, const UnaryPred& predicate) {
		std::copy_if(std::make_move_iterator(source.begin()), std::make_move_iterator(source.end()), std::back_inserter(target), predicate);
	}
}
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_DETAIL_CONTAINER_HELPER_HPP
