#include "../include/Utility/Numerics/Minimization/Bisection.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

using namespace Utility::Numerics::Minimization;

int main(){
    auto f1 = [](double x){
        return std::exp(-x*x) + 1e-3 * x * x;
    }; // x_0 = +/- sqrt(ln(1000))
    const double x1 = std::sqrt(std::log(1e3));
    const double x1_bi = bisection(f1, 0., 3., sqrt(std::numeric_limits<double>::epsilon()), 200);
    if(std::abs(x1 - x1_bi) < sqrt(std::numeric_limits<double>::epsilon())){
        std::cout << "Test 1 passed! ";
    } else {
        return 1;
    }
    std::cout << std::scientific << std::setprecision(16) << "x0 = " << x1 << "   we found   x0 = " << x1_bi << std::endl;
    return 0;
}