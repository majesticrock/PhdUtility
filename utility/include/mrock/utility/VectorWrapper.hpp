#pragma once
#include <vector>
#include <utility>
#include "detail/vector_macro.h"
/* Provides a wrapper for std::vector that can savely be derived from, i.e., provides a virtual deconstructor
* Forwards function calls to the contained std::vector, but only implements those
* that I needed so far. If others are required they can easily be added in the same manner
*/

namespace mrock::utility {
	template<class T, class Allocator = std::allocator<T>>
	struct VectorWrapper {
	public:
		std::vector<T, Allocator> _vector;

		virtual ~VectorWrapper() = default;
		VectorWrapper() = default;
		VectorWrapper(const std::vector<T, Allocator>& vector) : _vector(vector) {};
		VectorWrapper(std::vector<T, Allocator>&& vector) : _vector(std::move(vector)) {};
		explicit VectorWrapper(size_t size) : _vector(size) {};
		VectorWrapper(size_t size, const T& value) : _vector(size, value) {};
		VectorWrapper(std::initializer_list<T> init, const Allocator& alloc = Allocator()) : _vector(init, alloc) {};

		MROCK_VECTOR_WRAPPER_FILL_MEMBERS(T, _vector)

		inline auto operator<=>(const VectorWrapper<T, Allocator>& rhs) const = default;
	};
}