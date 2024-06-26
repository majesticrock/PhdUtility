#pragma once
#include <array>
#include <cassert>
#include <utility>
#include <type_traits>

namespace Utility::Numerics {
	namespace __detail {
		template<class _container>
		using element_t = std::decay_t<decltype(std::declval<_container>()[0])>;
	}

	template <class x_type, class y_type>
	constexpr y_type linearly_interpolate(x_type const& x, x_type const& x0, x_type const& x1, y_type const& y0, y_type const& y1) 
	{
		return (y0 * (x - x1) - y1 * (x - x0)) / (x0 - x1);
	}

	template <unsigned int n, class x_type, class y_type>
	constexpr y_type interpolate_lagrange(x_type const& x, std::array<x_type, n> const& x_data, std::array<y_type, n> const& y_data) 
	{
		static_assert(n > 1, "You need to provide at least 2 points!");
		y_type value{};
		for (unsigned int i = 0U; i < n; ++i) {
			x_type buffer{ 1 };
			for (unsigned int j = 0U; j < n; ++j) {
				if (i == j) continue;
				buffer *= (x - x_data[j]) / (x_data[i] - x_data[j]);
			}
			value += y_data[i] * buffer;
		}

		return value;
	}

	template <unsigned int n, class x_vector, class y_vector>
	inline __detail::element_t<y_vector> interpolate_from_vector(__detail::element_t<x_vector> const& x, 
		x_vector const& x_data, y_vector const& y_data, unsigned int start_index_x, int index_offset_y = 0)
	{
		assert(x_data.size() >= n && y_data.size() >= n + index_offset_y);

		if(start_index_x + n >= x_data.size()) {
			start_index_x = x_data.size() - n;
		}

		std::array<__detail::element_t<x_vector>, n> x_arr;
		std::array<__detail::element_t<y_vector>, n> y_arr;
		for(unsigned int i = 0U; i < n; ++i) {
			x_arr[i] = x_data[start_index_x + i];
			y_arr[i] = y_data[start_index_x + i + index_offset_y];
		}

		return interpolate_lagrange<n>(x, x_arr, y_arr);
	}
}