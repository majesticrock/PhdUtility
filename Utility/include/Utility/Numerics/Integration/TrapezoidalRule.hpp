#pragma once
#include <utility>

namespace Utility::Numerics::Integration {
	template <class UnaryFunction, class RealType, class... Args>
	auto trapezoidal_rule(const UnaryFunction& function, const RealType begin, const RealType end, unsigned long num_steps, Args&&... args) {
		const RealType step = (end - begin) / num_steps;
		decltype(function(begin)) value = 0.5 * (function(begin, std::forward<Args>(args)...) + function(end, std::forward<Args>(args)...));
		for (int n = 1U; n < num_steps; ++n) {
			value += function(begin + n * step, std::forward<Args>(args)...);
		}
		value *= step;
		return value;
	}

	template <class UnaryFunction, class RealType>
	auto trapezoidal_rule(const UnaryFunction& function, const RealType begin, const RealType end, unsigned long num_steps) {
		const RealType step = (end - begin) / num_steps;
		decltype(function(begin)) value{ 0.5 * (function(begin) + function(end)) };
		for (unsigned long n = 1U; n < num_steps; ++n) {
			value += function(begin + n * step);
		}

		value *= step;
		return value;
	}

	template <class UnaryFunction, class RealType>
	auto trapezoidal_rule_kahan(const UnaryFunction& function, const RealType begin, const RealType end, unsigned long num_steps) {
		const RealType step = (end - begin) / num_steps;
		decltype(function(begin)) y{ 0.5 * (function(begin) + function(end)) };
		decltype(function(begin)) t{ y };
		decltype(function(begin)) c{ t - y };

		decltype(function(begin)) value{ t };
		for (unsigned long n = 1U; n < num_steps; ++n) {
			y = function(begin + n * step) - c;
			t = value + y;
			c = (t - value) - y;
			value = t;
		}

		value *= step;
		return value;
	}

	template <class Vector, class RealType>
	auto trapezoidal_rule(const Vector& fx, const RealType step) {
		auto value = 0.5 * (fx.front() + fx.back());
		for (unsigned long n = 1U; n + 1U < fx.size(); ++n) {
			value += fx[n];
		}
		value *= step;
		return value;
	}
}