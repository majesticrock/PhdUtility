#pragma once
#include <utility>
#include <iostream>
#include <cmath>
#include <limits>

namespace Utility::Numerics::Roots {
	template <class Function, class RealType>
	RealType bisection(const Function& function, RealType begin, RealType end, RealType tol, int maxiter) {
        const auto is_zero = [](RealType val) {
            return std::abs(val) <= std::numeric_limits<RealType>::epsilon();
            };

        const RealType f_upper{ function(end) };
        if(is_zero(f_upper)) return end;
        const RealType f_lower{ function(begin) };
        if(is_zero(f_lower)) return begin;

        if(f_lower * f_upper > 0){
            throw std::invalid_argument("There is no root in the given interval!");
        }
        if(f_lower > 0) {
            // Ensure that the function is rising
            std::swap(begin, end);
        }
        RealType f_middle{ };
        RealType middle{ };
        do {
            middle = 0.5 * (end + begin);
            f_middle = function(middle);
            if(is_zero(f_middle)) return middle;
            if(f_middle < 0){
                begin = middle;
            } else {
                end = middle;
            }
        } while(std::abs(end - begin) > tol && --maxiter >= 0);
        if (maxiter < 0) {
			std::cerr << "Bisection terminated by maxiter-constraint! Function type:" << typeid(Function).name() << std::endl;
		}
        return middle;
    }
}