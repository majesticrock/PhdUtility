#pragma once
#include <vector>
#include <utility>
/* Provides a wrapper for std::vector that can savely be derived from, i.e., provides a virtual deconstructor
* Forwards function calls to the contained std::vector, but only implements those
* that I needed so far. If others are required they can easily be added in the same manner
*/

namespace Utility {
	template<class T, class Allocator = std::allocator<T>>
	class VectorWrapper {
	protected:
		std::vector<T, Allocator> _vector;

	public:
		using value_type = typename std::vector<T>::value_type;
		using allocator_type = typename std::vector<T>::allocator_type;
		using size_type = typename std::vector<T>::size_type;
		using difference_type = typename std::vector<T>::difference_type;
		using reference = typename std::vector<T>::reference;
		using const_reference = typename std::vector<T>::const_reference;
		using pointer = typename std::vector<T>::pointer;
		using const_pointer = typename std::vector<T>::const_pointer;
		using iterator = typename std::vector<T>::iterator;
		using const_iterator = typename std::vector<T>::const_iterator;
		using reverse_iterator = typename std::vector<T>::reverse_iterator;
		using constreverse_iterator = typename std::vector<T>::const_reverse_iterator;

		virtual ~VectorWrapper() = default;
		VectorWrapper() = default;
		VectorWrapper(const std::vector<T, Allocator>& vector) : _vector(vector) {};
		VectorWrapper(std::vector<T, Allocator>&& vector) : _vector(std::move(vector)) {};
		explicit VectorWrapper(size_t size) : _vector(size) {};
		VectorWrapper(size_t size, const T& value) : _vector(size, value) {};
		VectorWrapper(std::initializer_list<T> init, const Allocator& alloc = Allocator()) : _vector(init, alloc) {};

		inline reference operator[](size_t i) {
			return _vector[i];
		};
		inline const_reference operator[](size_t i) const {
			return _vector[i];
		};

		inline auto begin() noexcept {
			return _vector.begin();
		}
		inline auto begin() const noexcept {
			return _vector.begin();
		}
		inline auto end() noexcept {
			return _vector.end();
		}
		inline auto end() const noexcept {
			return _vector.end();
		}
		inline auto rbegin() noexcept {
			return _vector.rbegin();
		}
		inline auto rbegin() const noexcept {
			return _vector.rbegin();
		}
		inline auto rend() noexcept {
			return _vector.rend();
		}
		inline auto rend() const noexcept {
			return _vector.rend();
		}

		inline bool empty() const noexcept {
			return _vector.empty();
		}
		inline size_t size() const noexcept {
			return _vector.size();
		}

		inline const_reference front() const {
			return _vector.front();
		}
		inline reference front() {
			return _vector.front();
		}
		inline const_reference back() const {
			return _vector.back();
		}
		inline reference back() {
			return _vector.back();
		}

		inline void push_back(const T& element) {
			_vector.push_back(element);
		}
		inline void push_back(T&& element) {
			_vector.push_back(std::move(element));
		}
		inline void reserve(size_t new_capacity) {
			_vector.reserve(new_capacity);
		}
		inline void resize(size_t new_size) {
			_vector.resize(new_size);
		}
		inline void clear() noexcept {
			_vector.clear();
		}
		inline void pop_back() {
			_vector.pop_back();
		}

		template <class iterator>
		inline auto erase(iterator pos) {
			return _vector.erase(pos);
		}
		template <class iterator>
		inline auto erase(iterator first, iterator last) {
			return _vector.erase(first, last);
		}
		template <class iterator>
		inline auto insert(iterator pos, const T& value) {
			return _vector.insert(pos, value);
		}
		template <class iterator, class input_iterator>
		inline auto insert(iterator pos, input_iterator first, input_iterator last) {
			return _vector.insert(pos, first, last);
		}

		inline auto operator<=>(const VectorWrapper<T, Allocator>& rhs) const = default;
	};
}