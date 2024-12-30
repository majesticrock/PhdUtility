#pragma once
#include <utility>
#include <iostream>

namespace mrock::utility::Numerics::Minimization {
	template <class Function, class RealType>
	RealType bisection(const Function& function, RealType begin, RealType end, RealType tol, int maxiter) {
		if (begin > end) {
			std::swap(begin, end);
		}
		RealType middle = 0.5 * (begin + end);
		RealType d{};
		while (end - begin > tol && --maxiter >= 0) {
			if (middle - begin > end - middle) {
				d = (begin + middle) * 0.5;
				if (function(d) < function(middle)) {
					end = middle;
					middle = d;
				}
				else {
					begin = d;
				}
			}
			else {
				d = (middle + end) * 0.5;
				if (function(d) < function(middle)) {
					begin = middle;
					middle = d;
				}
				else {
					end = d;
				}
			}
		}
		if (maxiter < 0) {
			std::cerr << "Bisection terminated by maxiter-constraint!" << std::endl;
		}
		return RealType{ 0.5 } *(begin + end);
	}
}