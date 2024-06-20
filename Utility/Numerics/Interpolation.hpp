#pragma once
#include <array>

namespace Utility::Numerics {
	template <class x_type, class y_type>
	constexpr y_type linearly_interpolate(x_type const& x, x_type const& x0, x_type const& x1, y_type const& y0, y_type const& y1) {
		return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
	}

	template <class x_type, class y_type, unsigned int n>
	constexpr y_type interpolate_lagrange(x_type const& x, std::array<x_type, n> const& x_i, std::array<y_type, n> const& y_i) {
		static_assert(n > 1, "You need to provide at least 2 points!");
		y_type value{};
		for (unsigned int i = 0; i < n; ++i) {
			x_type buffer{};
			if (i == 0) {
				buffer = (x - x_i[1]) / (x_i[i] - x_i[1]);
				for (unsigned int j = 2; j < n; ++j) {
					buffer *= (x - x_i[j]) / (x_i[i] - x_i[j]);
				}
			}
			else {
				buffer = (x - x_i[0]) / (x_i[i] - x_i[0]);
				for (unsigned int j = 0; j < n; ++j) {
					if (i == j) continue;
					buffer *= (x - x_i[j]) / (x_i[i] - x_i[j]);
				}
			}
			value += y_i[i] * buffer;
		}
	}
}