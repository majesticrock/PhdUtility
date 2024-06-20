/* This file serves as a testing environment for the various
* algorithms provided in this library
*/
#define _USE_MATH_DEFINES

#include "Numerics/Integration/MidpointRule.hpp"
#include "Numerics/Integration/TrapezoidalRule.hpp"
#include "Numerics/Minimization/Bisection.hpp"

#include <cmath>
#include <iostream>

using namespace Utility;
using namespace Utility::Numerics;

int main(int argc, char** argv) {
	auto test_func = [](double t, double f) {
		return cos(M_PI * f * t);
		};
	auto test_func_2 = [](double t) {
		return cos(M_PI * t);
		};

	constexpr int step_num = 1000;

	std::cout << "Int Test 1: f=1:\n";
	std::cout << Integration::midpoint_rule(test_func, 1., 5., step_num, 1) << "\t"
		<< Integration::trapezoidal_rule(test_func, 1., 5., step_num, 1) << std::endl;

	std::cout << "Int Test 1: f=0.1:\n";
	std::cout << Integration::midpoint_rule(test_func, 1., 2., step_num, 0.1) << "\t"
		<< Integration::trapezoidal_rule(test_func, 1., 2., step_num, 0.1) << std::endl;

	std::cout << "Int Test 2: f=0.1:\n";
	std::cout << Integration::midpoint_rule(test_func_2, 0.3, 0.3, step_num) << "\t"
		<< Integration::trapezoidal_rule(test_func_2, 0.3, 0.3, step_num) << std::endl;

	auto min_test = [](double x) {return x * x - x * std::exp(-3 * x); };
	std::cout << "\n#####################################\n\nMin Test:\n";
	std::cout << Minimization::bisection(min_test, -3., 6., 1e-9, 100) << "\tAnalytical result: 0.16036674296 -> error = "
		<< Minimization::bisection(min_test, -3., 6., 1e-9, 100) - 0.16036674296 << std::endl;

	std::cout << "Min on boundary:\n";
	std::cout << Minimization::bisection([](double x) {return x * (x - 2); }, 0., 1., 1e-10, 100) << "\texpected: 1" << std::endl;
	std::cout << Minimization::bisection([](double x) {return -x * (x - 2); }, 0., 1., 1e-10, 100) << "\texpected: 0 or 1" << std::endl;
	return 0;
}