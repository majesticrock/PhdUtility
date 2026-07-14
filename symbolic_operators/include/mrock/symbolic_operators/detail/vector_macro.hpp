#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_DETAIL_VECTOR_MACRO_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_DETAIL_VECTOR_MACRO_HPP
/**
 * @file vector_macro.hpp
 * @brief This file contains macros for defining wrapper types and member functions for vector-like containers.
 *
 * The macros provided in this file help in creating wrapper classes for vector-like containers by defining
 * common type aliases and member functions. These macros can be used to reduce boilerplate code when creating
 * custom vector wrappers.
 *
 * Macros:
 * - \c MROCK_VECTOR_WRAPPER_FILL_MEMBERS(_T,_vector_name): Defines type aliases and common member functions for a \c std::vector<_T>.
 *
 * Example usage:
 * @code
 * class MyVectorWrapper {
 *     std::vector<int> my_vector;
 * 	   MROCK_VECTOR_WRAPPER_FILL_MEMBERS(int, my_vector)
 * };
 * @endcode
 */

#include <vector>
#include <cstddef>

#define MROCK_VECTOR_WRAPPER_FILL_MEMBERS(_T, _vector_name) using value_type = typename std::vector<_T>::value_type; \
		using allocator_type = typename std::vector<_T>::allocator_type; \
		using size_type = typename std::vector<_T>::size_type; \
		using difference_type = typename std::vector<_T>::difference_type; \
		using reference = typename std::vector<_T>::reference; \
		using const_reference = typename std::vector<_T>::const_reference; \
		using pointer = typename std::vector<_T>::pointer; \
		using const_pointer = typename std::vector<_T>::const_pointer; \
		using iterator = typename std::vector<_T>::iterator; \
		using const_iterator = typename std::vector<_T>::const_iterator; \
		using reverse_iterator = typename std::vector<_T>::reverse_iterator; \
		using constreverse_iterator = typename std::vector<_T>::const_reverse_iterator; \
		\
		inline reference operator[](size_type i) { \
			return _vector_name[i]; \
		}; \
		inline const_reference operator[](size_type i) const { \
			return _vector_name[i]; \
		}; \
 		\
		inline auto begin() noexcept { \
			return _vector_name.begin(); \
		} \
		inline auto begin() const noexcept { \
			return _vector_name.begin(); \
		} \
		inline auto end() noexcept { \
			return _vector_name.end(); \
		} \
		inline auto end() const noexcept { \
			return _vector_name.end(); \
		} \
		inline auto rbegin() noexcept { \
			return _vector_name.rbegin(); \
		} \
		inline auto rbegin() const noexcept { \
			return _vector_name.rbegin(); \
		} \
		inline auto rend() noexcept { \
			return _vector_name.rend(); \
		} \
		inline auto rend() const noexcept { \
			return _vector_name.rend(); \
		} \
 		\
		inline bool empty() const noexcept { \
			return _vector_name.empty(); \
		} \
		inline size_type size() const noexcept { \
			return _vector_name.size(); \
		} \
 		\
		inline const_reference front() const { \
			return _vector_name.front(); \
		} \
		inline reference front() { \
			return _vector_name.front(); \
		} \
		inline const_reference back() const { \
			return _vector_name.back(); \
		} \
		inline reference back() { \
			return _vector_name.back(); \
		} \
 		\
		inline void push_back(const value_type& element) { \
			_vector_name.push_back(element); \
		} \
		inline void push_back(value_type&& element) { \
			_vector_name.push_back(std::move(element)); \
		} \
		inline void reserve(size_type new_capacity) { \
			_vector_name.reserve(new_capacity); \
		} \
		inline void resize(size_type new_size) { \
			_vector_name.resize(new_size); \
		} \
		inline void resize(size_type new_size, value_type new_value) { \
			_vector_name.resize(new_size, new_value); \
		} \
		inline void clear() noexcept { \
			_vector_name.clear(); \
		} \
		inline void pop_back() { \
			_vector_name.pop_back(); \
		} \
 		\
		template <class iterator> \
		inline auto erase(iterator pos) { \
			return _vector_name.erase(pos); \
		} \
		template <class iterator> \
		inline auto erase(iterator first, iterator last) { \
			return _vector_name.erase(first, last); \
		} \
		template <class iterator> \
		inline auto insert(iterator pos, const value_type& value) { \
			return _vector_name.insert(pos, value); \
		} \
		template <class iterator, class input_iterator> \
		inline auto insert(iterator pos, input_iterator first, input_iterator last) { \
			return _vector_name.insert(pos, first, last); \
		} \

#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_DETAIL_VECTOR_MACRO_HPP