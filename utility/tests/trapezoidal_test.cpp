#include <iostream>
#include <cmath>
#include <iomanip>
#include "../include/mrock/utility/Numerics/Integration/AdaptiveTrapezoidalRule.hpp"
#include <Eigen/Dense>

using namespace mrock::utility::Numerics::Integration;

int main() {
    std::cout << std::scientific << std::setprecision(16) << std::endl;

    auto cos_func = [](double x) { return std::cos(x); };
    auto inv_one_plus_cos_func = [](double x) { return 1 / (1 + std::cos(x)); };
    auto x_squared_exp_x_func = [](double x) { return x * x * std::exp(x); };

    double begin = 0.0;
    double end = 1.0;
    unsigned int num_steps = 1000;
    double max_error = 1e-8;

    double cos_analytical = std::sin(end) - std::sin(begin);
    double inv_one_plus_cos_func_analytical = std::tan(0.5 * end) - std::tan(0.5 * begin); 
    double x_squared_exp_x_analytical = (end*end - 2*end + 2) * std::exp(end) - (begin*begin - 2*begin + 2) * std::exp(begin);

    adapative_trapezoidal_rule<double, false> adaptive;
    num_steps = 10;
    double cos_adaptive = adaptive.integrate(cos_func, begin, end, num_steps, max_error);
    double inv_one_plus_cos_func_adaptive = adaptive.integrate(inv_one_plus_cos_func, begin, end, num_steps, max_error);
    double x_squared_exp_x_adaptive = adaptive.integrate(x_squared_exp_x_func, begin, end, num_steps, max_error);

    std::cout << "Integral of cos(x) from 0 to 1: " << cos_adaptive << " (Analytical: " << cos_analytical << ")" << std::endl;
    std::cout << "Integral of sin(x)/x from 0 to 1: " << inv_one_plus_cos_func_adaptive << " (Analytical: " << inv_one_plus_cos_func_analytical << ")" << std::endl;
    std::cout << "Integral of x^2 * e^x from 0 to 1: " << x_squared_exp_x_adaptive << " (Analytical: " << x_squared_exp_x_analytical << ")" << std::endl;

    relative_error<double> error_func;
    if (error_func(cos_adaptive, cos_analytical) > max_error) {
        return 1;
    }
    if (error_func(inv_one_plus_cos_func_adaptive, inv_one_plus_cos_func_analytical) > max_error) {
        return 1;
    }
    if (error_func(x_squared_exp_x_adaptive, x_squared_exp_x_analytical) > max_error) {
        return 1;
    }

    auto vector_func = [](double x) { 
        return Eigen::Vector3d{std::tan(x), x * std::cos(x), (std::sin(x) + std::cos(x)) * std::exp(-x)}; 
    };
    Eigen::Vector3d vector_analytical = { 
        std::log(std::abs(std::cos(begin))) - std::log(std::abs(std::cos(end))),
        end * std::sin(end) + std::cos(end) - begin * std::sin(begin) - std::cos(begin), 
        std::cos(begin) * std::exp(-begin) - std::cos(end) * std::exp(-end) 
    };

    vector_norm_error<Eigen::Vector3d> vector_error;
    vector_elementwise_error<Eigen::Vector3d> alt_vector_error;

    Eigen::Vector3d vector_adaptive = adaptive.integrate(vector_func, begin, end, num_steps, max_error, vector_error, Eigen::Vector3d::Zero());
    
    std::cout << "Integrating a vector valued function:\n" << vector_adaptive << "\nAnalytical:\n" << vector_analytical 
        << "\nError = " << vector_error(vector_analytical, vector_adaptive) << std::endl;
    std::cout << "Elementwise error = " << alt_vector_error(vector_analytical, vector_adaptive) << std::endl;

    if (vector_error(vector_analytical, vector_adaptive) > max_error) {
        return 1;
    }

    return 0;
}