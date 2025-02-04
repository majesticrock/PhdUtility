#pragma once
#include <iterator>
#include <algorithm>

namespace mrock::utility {
	template<class Vector>
	void duplicate_n_inplace(Vector& target, size_t n) {
		const size_t original_size = target.size();
		target.reserve(original_size * (n + 1U));
		for (size_t i = 0U; i < n; ++i) {
			std::copy_n(target.begin(), original_size, std::back_inserter(target));
		}
	}

	template<class Vector>
	void append_vector(Vector& target, const Vector& source) {
		target.insert(target.end(), source.begin(), source.end());
	}
	template<class Vector>
	void append_vector(Vector& target, Vector&& source) {
		target.insert(target.end(), std::make_move_iterator(source.begin()), std::make_move_iterator(source.end()));
	}
	template<class Vector, class UnaryPred>
	void append_if(Vector& target, const Vector& source, const UnaryPred& predicate) {
		std::copy_if(source.begin(), source.end(), std::back_inserter(target), predicate);
	}
	template<class Vector, class UnaryPred>
	void append_if(Vector& target, Vector&& source, const UnaryPred& predicate) {
		std::copy_if(std::make_move_iterator(source.begin()), std::make_move_iterator(source.end()), std::back_inserter(target), predicate);
	}
}